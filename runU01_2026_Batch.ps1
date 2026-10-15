<#
runU01_2026_Batch.ps1 - batch driver for the 2026 U01 reprocessing.

Walks every session under -DataRoot (including the nested Repeatability\
folder), decides whether it is a brain or a calf session, and hands it to the
appropriate existing driver:

    brain -> .\processU01_Brain.ps1   (anatomic prep + HD-BET + brainseg.sh + 1H + 31P)
    calf  -> .\processU01_Calf.ps1    (1H + 31P only; no segmentation)

This script only dispatches -- it does not duplicate any processing logic.

SKIP RULES (each skipped session is reported, with its reason, in the summary
and in the manifest CSV):
  * phantom-only sessions (no subject)
  * sessions that lack EITHER a 1H (*svssel*.dat) or a 31P (*31p*.dat) raw file
    (with -Skip1H / -Skip31P only the modality being processed is required)
  * brain sessions with no MPRAGE* series -- brainseg.sh cannot run without a
    T1w, and processU01_Brain.ps1 asserts on it at stage 1 (with -Skip1H only
    -ExtRef sessions need it, for the HD-BET mask)

EXTERNAL REFERENCE:
  -ExtRef is passed through to the driver only for sessions whose folder name
  marks them as having an external phantom (external / d2o / pbs / extref).
  Enabling it elsewhere would add a spurious ~8 ppm singlet to the 31P basis.
  Sessions so marked but lacking a B1 map still get RefExt FITTED; they just
  cannot be converted to concentrations (extRefQuant31P reports ok=false).

  NOTE: -RefConc defaults to 500 mM.  External concentrations scale LINEARLY
  with it, and the D2O / PBS / phosphate phantoms are not necessarily all
  500 mM.  Set it per phantom before trusting the external numbers.

Usage:
    .\runU01_2026_Batch.ps1 -DryRun            # show the plan, process nothing
    .\runU01_2026_Batch.ps1                    # process everything
    .\runU01_2026_Batch.ps1 -Only 'S095*'      # subset by session-name wildcard
    .\runU01_2026_Batch.ps1 -Anatomy calf      # calf sessions only
    .\runU01_2026_Batch.ps1 -Skip31P           # 1H only, every session
    .\runU01_2026_Batch.ps1 -Skip1H            # 31P only, every session
    .\runU01_2026_Batch.ps1 -Skip31P -ReuseAnat  # 1H only, brain anatomy from processed\

Note: brain sessions need the VirtualBox VM running (brainseg.sh over ssh),
except with -Skip1H.  The script preflights this and, if the VM is
unreachable, skips brain sessions rather than failing 7 times in a row.
#>
param(
    [string]$DataRoot = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\2026',
    [string]$LogRoot  = 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_2026\_batchlogs',
    [string]$Only     = '*',                      # wildcard filter on session name
    [ValidateSet('all','brain','calf')]
    [string]$Anatomy  = 'all',
    [double]$RefConc  = 500,                      # external phantom [mM] -- see note above
    [switch]$DryRun,                              # classify + report, run nothing
    [switch]$SkipExisting,                        # skip sessions that already have processed\*_abs.mat
    [switch]$Skip1H,                              # passed to the drivers: 31P only
    [switch]$Skip31P,                             # passed to the drivers: 1H only
    [switch]$ReuseAnat,                           # passed to the brain driver: reuse processed\ anatomy
    [switch]$Plot
)

if ($Skip1H -and $Skip31P) { throw "-Skip1H and -Skip31P together leave nothing to process." }

$ErrorActionPreference = 'Stop'
$scriptDir = $PSScriptRoot
$batchStart = Get-Date
$stamp = $batchStart.ToString('yyyyMMdd_HHmmss')

if (-not (Test-Path $DataRoot)) { throw "DataRoot not found: $DataRoot" }
if (-not $DryRun) { New-Item -ItemType Directory -Force -Path $LogRoot | Out-Null }

# ---- Session discovery -------------------------------------------------
# Top-level session folders, plus one level down inside Repeatability\.
function Get-SessionDirs([string]$root) {
    $dirs = @(Get-ChildItem -Path $root -Directory)
    $out  = @()
    foreach ($d in $dirs) {
        if ($d.Name -eq 'Repeatability') {
            $out += @(Get-ChildItem -Path $d.FullName -Directory)
        } else {
            $out += $d
        }
    }
    return $out
}

# ---- Classification ----------------------------------------------------
function Get-ExamDirs([string]$sessionDir) {
    return @(Get-ChildItem -Path $sessionDir -Directory -Filter 'E*' -ErrorAction SilentlyContinue)
}

function Test-SeriesExists([string]$sessionDir, [string]$pattern) {
    foreach ($e in Get-ExamDirs $sessionDir) {
        if (@(Get-ChildItem -Path $e.FullName -Directory -Filter $pattern -ErrorAction SilentlyContinue).Count -gt 0) {
            return $true
        }
    }
    return $false
}

function Get-DatCount([string]$sessionDir, [string]$pattern) {
    $dd = Join-Path $sessionDir 'datfiles'
    if (-not (Test-Path $dd)) { return 0 }
    return @(Get-ChildItem -Path $dd -File -Filter $pattern -ErrorAction SilentlyContinue).Count
}

function Test-Mixed31P([string]$sessionDir) {
    # True when the session holds BOTH real 31P spectroscopy (isis_svs*) and
    # some other *31p*.dat (prep_moco_b1map_31P, fid_test1D_31P, ...).
    #
    # pick_latest_31p inside the MATLAB drivers takes the highest-MID
    # *31p*.dat, and the non-spectroscopy files often carry a HIGHER MID than
    # the real scan -- so it selects a file with no spectroscopy dwell time and
    # dat2mat_svs31p dies with 'unknown dwell time'.  Rather than guess which
    # file was intended, such sessions are skipped for manual sorting.
    $dd = Join-Path $sessionDir 'datfiles'
    if (-not (Test-Path $dd)) { return $false }
    $all = @(Get-ChildItem -Path $dd -File -Filter '*31p*.dat' -ErrorAction SilentlyContinue)
    if ($all.Count -lt 2) { return $false }
    $isis  = @($all | Where-Object { $_.Name -match '(?i)isis' })
    $other = @($all | Where-Object { $_.Name -notmatch '(?i)isis' })
    return ($isis.Count -gt 0 -and $other.Count -gt 0)
}

function Classify-Session($sess) {
    $name = $sess.Name
    $dir  = $sess.FullName

    $info = [ordered]@{
        Session  = $name
        Path     = $dir
        Anatomy  = ''
        N1H      = 0
        N31P     = 0
        HasMprage= $false
        HasB1    = $false
        ExtRef   = $false
        Action   = ''
        Reason   = ''
    }

    # phantom-only sessions carry no subject
    if ($name -match 'PHANTOM') {
        $info.Action = 'SKIP'; $info.Reason = 'phantom-only session (no subject)'
        return [pscustomobject]$info
    }

    $info.N1H  = Get-DatCount $dir '*svssel*.dat'
    $info.N31P = Get-DatCount $dir '*31p*.dat'

    # anatomy: folder name is authoritative, DICOM series names corroborate
    $isCalf = ($name -match '(?i)calf') -or
              (Test-SeriesExists $dir '*CALF*') -or
              (Test-SeriesExists $dir '*DIXON*')
    $info.Anatomy = if ($isCalf) { 'calf' } else { 'brain' }

    $info.HasMprage = Test-SeriesExists $dir 'MPRAGE*'
    $info.HasB1     = Test-SeriesExists $dir '*B1MAP*'
    # external phantom present?  folder-name marker is the only reliable signal
    $info.ExtRef    = ($name -match '(?i)external|d2o|pbs|extref')

    if ($info.N1H -eq 0 -and $info.N31P -eq 0) {
        $info.Action = 'SKIP'; $info.Reason = 'no 1H and no 31P raw data'
    } elseif (-not $Skip1H -and $info.N1H -eq 0) {
        $info.Action = 'SKIP'; $info.Reason = 'no 1H (*svssel*.dat) -- 31P only (use -Skip1H)'
    } elseif (-not $Skip31P -and $info.N31P -eq 0) {
        $info.Action = 'SKIP'; $info.Reason = 'no 31P (*31p*.dat) -- 1H only (use -Skip31P)'
    } elseif (-not $Skip31P -and (Test-Mixed31P $dir)) {
        $info.Action = 'SKIP'
        $info.Reason = 'mixed 31P: isis_svs + non-spectroscopy *31p*.dat -- pick one manually'
    } elseif (-not $isCalf -and -not $info.HasMprage -and -not $Skip1H) {
        $info.Action = 'SKIP'
        $info.Reason = 'brain session with no MPRAGE* -- segmentation impossible'
    } elseif (-not $isCalf -and -not $info.HasMprage -and $info.ExtRef) {
        $info.Action = 'SKIP'
        $info.Reason = 'brain session with no MPRAGE* -- no HD-BET mask for -ExtRef'
    } else {
        $info.Action = 'RUN'
    }

    return [pscustomobject]$info
}

# ---- Build the plan ----------------------------------------------------
Write-Host ("=" * 78) -ForegroundColor Cyan
Write-Host "U01 2026 batch reprocessing" -ForegroundColor Cyan
Write-Host "  DataRoot : $DataRoot"
Write-Host "  Spectra  : $(if ($Skip1H) { '31P only' } elseif ($Skip31P) { '1H only' } else { '1H + 31P' })"
Write-Host "  Started  : $($batchStart.ToString('yyyy-MM-dd HH:mm:ss'))"
Write-Host ("=" * 78) -ForegroundColor Cyan

$plan = @()
foreach ($s in Get-SessionDirs $DataRoot) {
    if ($s.Name -notlike $Only) { continue }
    $plan += Classify-Session $s
}

if ($Anatomy -ne 'all') {
    $plan = @($plan | Where-Object { $_.Anatomy -eq $Anatomy -or $_.Action -eq 'SKIP' })
}

if ($SkipExisting) {
    foreach ($p in $plan) {
        if ($p.Action -ne 'RUN') { continue }
        $pd   = Join-Path $p.Path 'processed'
        $want = @()
        if (-not $Skip1H)  { $want += '{0}_{1}_1H_abs.mat'  -f $p.Session, $p.Anatomy }
        if (-not $Skip31P) { $want += '{0}_{1}_31P_abs.mat' -f $p.Session, $p.Anatomy }
        $have = @($want | Where-Object { Test-Path (Join-Path $pd $_) })
        if ($have.Count -eq $want.Count) {
            $p.Action = 'SKIP'; $p.Reason = "already processed ($($have -join ', ') present)"
        }
    }
}

$toRun  = @($plan | Where-Object { $_.Action -eq 'RUN' })
$toSkip = @($plan | Where-Object { $_.Action -eq 'SKIP' })

Write-Host "`n--- PLAN ---" -ForegroundColor Yellow
$plan | Format-Table -AutoSize Session, Anatomy, N1H, N31P, HasMprage, HasB1, ExtRef, Action, Reason | Out-String -Width 250 | Write-Host

Write-Host ("To run : {0}  ({1} brain, {2} calf)" -f $toRun.Count,
    @($toRun | Where-Object Anatomy -eq 'brain').Count,
    @($toRun | Where-Object Anatomy -eq 'calf').Count) -ForegroundColor Green
Write-Host ("Skipped: {0}" -f $toSkip.Count) -ForegroundColor Yellow

# ---- VM preflight (brain only) ----------------------------------------
$brainCount = @($toRun | Where-Object Anatomy -eq 'brain').Count
if ($brainCount -gt 0 -and -not $Skip1H -and -not $DryRun) {
    Write-Host "`nPreflight: checking VM 'ubuntu-vm' for brainseg.sh ..." -ForegroundColor Cyan
    # Must use a LOGIN shell (bash -lc), exactly as processU01_Brain.ps1's
    # Invoke-Vm does: brainseg.sh is only on the PATH set up by the profile.
    & ssh -o ConnectTimeout=15 ubuntu-vm "bash -lc 'command -v brainseg.sh'" 2>&1 | Out-Null
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  VM unreachable or brainseg.sh missing -- brain sessions will be SKIPPED." -ForegroundColor Red
        Write-Host "  Start the VirtualBox VM and re-run to process brain sessions." -ForegroundColor Red
        foreach ($p in $toRun) {
            if ($p.Anatomy -eq 'brain') { $p.Action = 'SKIP'; $p.Reason = 'VM unreachable (brainseg.sh)' }
        }
        $toRun  = @($plan | Where-Object { $_.Action -eq 'RUN' })
        $toSkip = @($plan | Where-Object { $_.Action -eq 'SKIP' })
    } else {
        Write-Host "  VM OK." -ForegroundColor Green
    }
}

if ($DryRun) {
    Write-Host "`n-DryRun set: nothing processed." -ForegroundColor Yellow
    return
}

# ---- Execute -----------------------------------------------------------
$manifest = @()
foreach ($p in $toSkip) {
    $manifest += [pscustomobject]@{
        Session = $p.Session; Anatomy = $p.Anatomy; Status = 'SKIPPED'
        ExtRef  = $p.ExtRef;  Minutes = 0;          Detail  = $p.Reason
    }
}

$i = 0
foreach ($p in $toRun) {
    $i++
    $hdr = "[{0}/{1}] {2}  ({3}{4})" -f $i, $toRun.Count, $p.Session, $p.Anatomy,
           $(if ($p.ExtRef) { ', extRef' } else { '' })
    Write-Host "`n$("=" * 78)" -ForegroundColor Cyan
    Write-Host $hdr -ForegroundColor Cyan
    Write-Host ("=" * 78) -ForegroundColor Cyan

    $driver = if ($p.Anatomy -eq 'calf') {
        Join-Path $scriptDir 'processU01_Calf.ps1'
    } else {
        Join-Path $scriptDir 'processU01_Brain.ps1'
    }

    $args = @{ SessionDir = $p.Path }
    if ($Plot)     { $args['Plot']    = $true }
    if ($Skip1H)   { $args['Skip1H']  = $true }
    if ($Skip31P)  { $args['Skip31P'] = $true }
    if ($ReuseAnat -and $p.Anatomy -eq 'brain') { $args['ReuseAnat'] = $true }
    if ($p.ExtRef) { $args['ExtRef']  = $true; $args['RefConc'] = $RefConc }

    $log  = Join-Path $LogRoot ("{0}_{1}.log" -f $p.Session, $stamp)
    $t0   = Get-Date
    $status = 'OK'; $detail = ''

    try {
        # Each driver is self-contained; a failure in one session must not
        # stop the batch, so ErrorActionPreference is relaxed for the call
        # and the transcript is kept per session.
        Start-Transcript -Path $log -Force | Out-Null
        try {
            & $driver @args
        } finally {
            Stop-Transcript | Out-Null
        }
    } catch {
        $status = 'FAILED'
        $detail = $_.Exception.Message -replace '\s+', ' '
        Write-Host "  FAILED: $detail" -ForegroundColor Red
        Write-Host "  log: $log" -ForegroundColor Red
    }

    $mins = [math]::Round(((Get-Date) - $t0).TotalMinutes, 1)
    if ($status -eq 'OK') {
        Write-Host ("  done in {0} min" -f $mins) -ForegroundColor Green
    }

    $manifest += [pscustomobject]@{
        Session = $p.Session; Anatomy = $p.Anatomy; Status = $status
        ExtRef  = $p.ExtRef;  Minutes = $mins;      Detail  = $detail
    }
}

# ---- Summary -----------------------------------------------------------
$elapsed = (Get-Date) - $batchStart
Write-Host "`n$("=" * 78)" -ForegroundColor Cyan
Write-Host "BATCH SUMMARY" -ForegroundColor Cyan
Write-Host ("=" * 78) -ForegroundColor Cyan
$manifest | Sort-Object Status, Anatomy, Session |
    Format-Table -AutoSize Session, Anatomy, Status, ExtRef, Minutes, Detail |
    Out-String -Width 250 | Write-Host

$manifestCsv = Join-Path $LogRoot "batch_manifest_$stamp.csv"
$manifest | Export-Csv -Path $manifestCsv -NoTypeInformation -Encoding utf8
Write-Host ("OK: {0}   FAILED: {1}   SKIPPED: {2}" -f
    @($manifest | Where-Object Status -eq 'OK').Count,
    @($manifest | Where-Object Status -eq 'FAILED').Count,
    @($manifest | Where-Object Status -eq 'SKIPPED').Count) -ForegroundColor Green
Write-Host ("Total run time: {0:hh\:mm\:ss}" -f $elapsed) -ForegroundColor Green
Write-Host "Manifest: $manifestCsv" -ForegroundColor Green
Write-Host "Logs:     $LogRoot" -ForegroundColor Green
