<#
processU01.ps1 - U01 SVS pipeline driver.

Architecture:
  - DICOMs and .dat read directly by MATLAB on Windows.
  - Anatomic prep NIfTI outputs (T1, T1_SPM8BIAS, svsmask) staged into the
    VirtualBox shared folder so brainseg.sh in the VM can see them.
  - HD-BET runs on the Windows host (PyTorch), writing _hdbet outputs
    into the same staging dir; brainseg.sh's own HD-BET block then skips.
  - FSL stages (FAST, FIRST, fslmaths, ROI stats) still run in the VM
    via brainseg.sh; outputs land alongside the staged inputs.
  - MATLAB quant outputs (*.mat) written directly to <SessionDir>\processed\.
  - Final step: copy non-.mat outputs from staging back to processed,
    then remove staging dir.

Usage:
    .\processU01.ps1 -SessionDir 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR\<scanID>'
	
Note: Make sure VM is running when this program is run.
#>
param(
    [Parameter(Mandatory)][string]$SessionDir,
    [string]$StagingRoot = 'C:\Users\CAMIPM-NW\Documents\VirtualBoxSharedFolders\CAMIPM_Ubuntu_VM2\TEMP',
    [string]$FinalOutDir = (Join-Path $SessionDir 'processed'),
    [switch]$KeepStaging,
    [switch]$Plot,         # -Plot turns on pars.plt (and pars.pltSpec for 1H)
    [switch]$ExtRef,       # -ExtRef adds the ~8 ppm phantom singlet to the 31P basis
                           #   and quantifies against it (in addition to aATP).
                           #   B1+ correction is OFF inside extRefQuant31P.
                           #   A missing B1 series is non-fatal.
    [double]$RefPpm  = 8.0,# external-reference chemical shift (PCr = 0)
    [double]$RefConc = 500 # external phantom concentration [mM]
)

if ($Plot) { $pltArg = 'true' } else { $pltArg = 'false' }

# 5th arg to processU01_Brain_31P: an extRef struct literal, or an empty
# struct() which leaves the feature off (enable defaults to false).
if ($ExtRef) {
    $extRefArg = "struct('enable',true,'refPpm',$RefPpm,'refConc_mM',$RefConc)"
} else {
    $extRefArg = "struct()"
}

$ErrorActionPreference = 'Stop'

$pipelineStart = Get-Date
Write-Host ("Pipeline start: {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineStart) -ForegroundColor Green

# ---- Configuration -----------------------------------------------------
$matlab     = 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe'
$hdbet      = 'hd-bet'        # hd-bet 2.x executable (must be on PATH)
$hdbetDev   = 'cuda'           # 'cuda' / 'cpu' / 'mps' (hd-bet 2.x device flag)
$vmHost     = 'ubuntu-vm'
$scanID     = Split-Path $SessionDir -Leaf
$svsCodeDir = $PSScriptRoot   # this script lives in svs_code/
$camipmDir  = Join-Path (Split-Path $svsCodeDir -Parent) 'CAMIPM'   # sibling of svs_code/

# When plotting, append (per MATLAB call) a snippet that first SAVES every open
# figure to $FinalOutDir as .fig + .png, then waits for them to be closed before
# the function returns (so -batch MATLAB doesn't exit and destroy them).  The
# save is programmatic because -batch has no desktop / File>Save As dialog.
if ($Plot) {
    $pltWait = " drawnow;" +
        " figDir='$FinalOutDir'; sid='$scanID';" +
        " figs=findall(0,'Type','figure');" +
        " for fi=1:numel(figs), fg=figs(fi);" +
        "  nm=get(fg,'Name'); if isempty(nm), nm=sprintf('fig%d',fi); end;" +
        "  nm=regexprep(nm,'[^\w-]','_');" +
        "  base=fullfile(figDir,sprintf('%s_%s_%02d',sid,nm,fi));" +
        "  try, savefig(fg,[base '.fig']); catch, end;" +
        "  try, exportgraphics(fg,[base '.png'],'Resolution',150); catch, end;" +
        " end;" +
        " fprintf('\n*** %d figure(s) saved to %s.\n*** Close all figures in MATLAB to continue pipeline ***\n',numel(figs),figDir);" +
        " while ~isempty(findall(0,'Type','figure')), pause(0.5); end;"
} else {
    $pltWait = ''
}

# Windows -> VM path mapping (only the staging area needs to round-trip).
$PathMap = [ordered]@{
    'C:\Users\CAMIPM-NW\Documents\VirtualBoxSharedFolders\CAMIPM_Ubuntu_VM2\' = '/media/sf_CAMIPM_Ubuntu_VM2/'
}

# ---- Helpers -----------------------------------------------------------
function Convert-ToVmPath([string]$winPath) {
    foreach ($k in $PathMap.Keys) {
        if ($winPath.StartsWith($k, [StringComparison]::OrdinalIgnoreCase)) {
            $tail = $winPath.Substring($k.Length) -replace '\\','/'
            return $PathMap[$k] + $tail
        }
    }
    throw "No VM path mapping for: $winPath"
}

function Invoke-Matlab([string]$Stmt) {
    Write-Host ">> MATLAB: $Stmt" -ForegroundColor Cyan
    & $matlab -batch $Stmt -wait
    if ($LASTEXITCODE -ne 0) { throw "MATLAB failed (exit $LASTEXITCODE)" }
}

function Invoke-Vm([string]$Cmd) {
    Write-Host ">> VM: $Cmd" -ForegroundColor Cyan
    $escaped = $Cmd -replace "'", "'\''"
    ssh $vmHost "bash -lc '$escaped'"
    if ($LASTEXITCODE -ne 0) { throw "VM command failed (exit $LASTEXITCODE)" }
}

# ---- Setup -------------------------------------------------------------
$stageDir = Join-Path $StagingRoot $scanID
New-Item -ItemType Directory -Force -Path $stageDir   | Out-Null
New-Item -ItemType Directory -Force -Path $FinalOutDir | Out-Null

# ---- Stage 1: anatomic prep (MATLAB on Windows) ----------------------------
# Reads DICOMs from $SessionDir on Z:\, writes NIfTIs to $stageDir.
Invoke-Matlab "addpath(genpath('$camipmDir')); addpath(genpath('$svsCodeDir')); prep_anatomic_U01_Brain_1H('$SessionDir','$stageDir'); rmpath(genpath('$camipmDir')); rmpath(genpath('$svsCodeDir'));"

$t1Bias  = (Get-ChildItem $stageDir -Filter '*_SPM8BIAS.nii' | Select-Object -First 1).FullName
$svsMask = Join-Path $stageDir 'svsmask.nii'
if (-not $t1Bias)              { throw "Stage 1: no *_SPM8BIAS.nii in $stageDir" }
if (-not (Test-Path $svsMask)) { throw "Stage 1: no svsmask.nii in $stageDir"   }

# ---- Stage 1b: HD-BET on the host --------------------------------------
# Writes <stem>_hdbet.nii.gz into $stageDir so the HD-BET block inside
# brainseg.sh sees it on disk and skips its own (slow) call.
# hd-bet 2.x requires .nii.gz input; our staged T1 is .nii, so gzip a
# sibling copy, run hd-bet, then delete it.
# We do not write _hdbet_mask: brainseg.sh only references the mask in
# its -U/-T branches, which we never trigger.
$t1BiasStem = [IO.Path]::GetFileNameWithoutExtension($t1Bias)   # e.g. T1_SPM8BIAS
$hdbetOut   = Join-Path $stageDir "${t1BiasStem}_hdbet.nii.gz"

if (Test-Path $hdbetOut) {
    Write-Host "HD-BET output already exists, skipping: $hdbetOut" -ForegroundColor Yellow
} else {
    $t1BiasGz = "$t1Bias.gz"
    Write-Host ">> gzip $t1Bias -> $t1BiasGz" -ForegroundColor Cyan
    $inStrm  = [System.IO.File]::OpenRead($t1Bias)
    $outStrm = [System.IO.File]::Create($t1BiasGz)
    $gzStrm  = New-Object System.IO.Compression.GZipStream($outStrm, [System.IO.Compression.CompressionMode]::Compress)
    try     { $inStrm.CopyTo($gzStrm) }
    finally { $gzStrm.Dispose(); $outStrm.Dispose(); $inStrm.Dispose() }

    Write-Host ">> HD-BET: $t1BiasGz -> $hdbetOut" -ForegroundColor Cyan
    # CPU variant (slower, slightly cleaner with TTA off): add --disable_tta
    & $hdbet -i $t1BiasGz -o $hdbetOut -device $hdbetDev
    $hdbetExit = $LASTEXITCODE
    Remove-Item -Force $t1BiasGz -ErrorAction SilentlyContinue
    if ($hdbetExit -ne 0) { throw "HD-BET failed (exit $hdbetExit)" }
}
if (-not (Test-Path $hdbetOut)) { throw "Stage 1b: HD-BET did not produce $hdbetOut" }

# ---- Stage 2: brainseg.sh in the VM ------------------------------------
# With _hdbet outputs pre-staged, brainseg.sh logs "Skipping hd-bet." and
# proceeds straight to FAST -> FIRST -> merge -> SVS mask -> ROI stats.
$t1BiasVm  = Convert-ToVmPath $t1Bias
$svsMaskVm = Convert-ToVmPath $svsMask
Invoke-Vm "brainseg.sh -f -b N -S '$svsMaskVm' '$t1BiasVm' '$scanID'"

$segStem = [IO.Path]::GetFileNameWithoutExtension($t1Bias)
$segTsv  = Get-ChildItem $stageDir -Filter "${segStem}_hdbet_fast3_first_seg_svsmask.tsv" |
              Select-Object -First 1 -ExpandProperty FullName
if (-not $segTsv) { throw "Stage 2: brainseg.sh did not produce expected .tsv" }
Write-Host "Seg TSV: $segTsv" -ForegroundColor Green

# ---- Stage 3: 1H NAD + absolute quant ----------------------------------
# Reads *svssel*.dat from $SessionDir\datfiles, tissue fractions from $segTsv,
# writes <scanID>_brain_1H_abs.mat directly to $FinalOutDir.
Invoke-Matlab "addpath('$svsCodeDir'); processU01_Brain_1H('$SessionDir','$FinalOutDir','$segTsv',$pltArg); rmpath('$svsCodeDir');$pltWait"

# ---- Stage 4: 31P quant ------------------------------------------------
# Reads *31p*.dat from $SessionDir\datfiles, uses f_CSF from $segTsv, writes
# <scanID>_brain_31P_abs.mat directly to $FinalOutDir.
Invoke-Matlab "addpath('$svsCodeDir'); processU01_Brain_31P('$SessionDir','$FinalOutDir','$segTsv',$pltArg,$extRefArg); rmpath('$svsCodeDir');$pltWait"

# ---- Stage 5: NAD summary to console -----------------------------------
# Reload the two abs.mat files and print NADH2/H4/H6 (1H, uM) and
# NADplus/NADH (31P, mM) to the window.
$mat1H  = Join-Path $FinalOutDir "${scanID}_brain_1H_abs.mat"
$mat31P = Join-Path $FinalOutDir "${scanID}_brain_31P_abs.mat"
$summaryStmt = "m1=load('$mat1H','absQ'); m31=load('$mat31P','absQ');" +
               " fprintf('\n=== NAD summary: ${scanID} ===\n');" +
               " for k=1:numel(m1.absQ.metabolites), n=m1.absQ.metabolites(k).name;" +
               "  if any(strcmp(n,{'NADH2','NADH4','NADH6','Trp'}))," +
               "   fprintf('  1H  %-8s : %.4g uM\n',n,m1.absQ.metabolites(k).conc_mM*1e3); end; end;" +
               " for k=1:numel(m31.absQ.metabolites), n=m31.absQ.metabolites(k).name;" +
               "  if any(strcmp(n,{'NADplus','NADH'}))," +
               "   fprintf('  31P %-8s : %.4g uM\n',n,m31.absQ.metabolites(k).conc_mM*1e3); end; end"
Invoke-Matlab $summaryStmt

# ---- Stage final: ship staging outputs to Z: and clean up --------------
Write-Host "Copying staging outputs to $FinalOutDir" -ForegroundColor Cyan
Copy-Item -Path (Join-Path $stageDir '*') -Destination $FinalOutDir -Recurse -Force

if (-not $KeepStaging) {
    Remove-Item -Recurse -Force $stageDir
} else {
    Write-Host "Staging preserved at $stageDir" -ForegroundColor Yellow
}

$pipelineEnd = Get-Date
$elapsed     = $pipelineEnd - $pipelineStart
Write-Host ("Pipeline end:   {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineEnd) -ForegroundColor Green
Write-Host ("Total run time: {0:hh\:mm\:ss} ({1:N1} min)" -f $elapsed, $elapsed.TotalMinutes) -ForegroundColor Green
Write-Host "Done: $FinalOutDir" -ForegroundColor Green
