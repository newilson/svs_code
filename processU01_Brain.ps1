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
  - -Skip1H drops the 1H quant and the FIRST / brainseg.sh segmentation
    that feeds it; anatomic prep + HD-BET still run with -ExtRef (V_M).
    The 31P f_CSF (QC only) comes from a seg TSV already in processed\,
    if there is one.  -Skip31P drops the 31P quant.
  - -ReuseAnat seeds staging with the anatomic outputs (NIfTI/TSV/TXT) of
    an earlier run in processed\, so prep, HD-BET, FIRST and every
    brainseg.sh module find their results on disk and skip.  Results made
    by older anatomic code are reused as-is.

Usage:
    .\processU01.ps1 -SessionDir 'Z:\Data\Spectroscopy\Brain 7T - 31P\U01_RR\<scanID>'
    .\processU01.ps1 -SessionDir '...\<scanID>' -Ref31P PCr
    .\processU01.ps1 -SessionDir '...\<scanID>' -Ref31P External -RefConc 500
    .\processU01.ps1 -SessionDir '...\<scanID>' -Skip31P    # 1H only
    .\processU01.ps1 -SessionDir '...\<scanID>' -Skip1H     # 31P only
    .\processU01.ps1 -SessionDir '...\<scanID>' -Skip31P -ReuseAnat

Note: Make sure VM is running when this program is run (not needed with -Skip1H).
#>
param(
    [Parameter(Mandatory)][string]$SessionDir,
    [string]$StagingRoot = 'C:\Users\CAMIPM-NW\Documents\VirtualBoxSharedFolders\CAMIPM_Ubuntu_VM2\TEMP',
    [string]$FinalOutDir = (Join-Path $SessionDir 'processed'),
    [switch]$KeepStaging,
    [switch]$Plot,         # -Plot turns on pars.plt (and pars.pltSpec for 1H)
    [switch]$Skip1H,       # -Skip1H skips the 1H quant and its segmentation stages
    [switch]$Skip31P,      # -Skip31P skips the 31P quant
    [switch]$ReuseAnat,    # -ReuseAnat reuses the anatomic outputs already in processed\
    [switch]$ExtRef,       # -ExtRef adds the ~8 ppm phantom singlet to the 31P basis
                           #   and quantifies against it (in addition to aATP).
                           #   A missing B1 series is non-fatal.
    [switch]$ApplyB1,      # -ApplyB1 applies the transmit efficiency sin(relB1*flip)
                           #   and the reciprocity receive factor to the external
                           #   quant. OFF by default: the map is 1H, not 31P, so it
                           #   is the right shape but not the right magnitude.
    [double]$RefPpm  = 8.0,# external-reference chemical shift (PCr = 0)
    [double]$RefConc = 500,# external phantom concentration [mM]
    [double]$RefArea = 200,# phantom cross-section [mm^2]; V_ref = area x slab (hdr)
    [ValidateSet('aATP','PCr','External')]
    [string]$Ref31P  = 'aATP'  # which reference drives absQ.metabolites
)

if ($Plot) { $pltArg = 'true' } else { $pltArg = 'false' }

if ($Skip1H -and $Skip31P) { throw "-Skip1H and -Skip31P together leave nothing to process." }

# -Ref31P External is incoherent without the external quant, so turn it on.
if ($Ref31P -eq 'External' -and -not $ExtRef) {
    Write-Host "-Ref31P External implies -ExtRef; enabling it." -ForegroundColor Yellow
    $ExtRef = [switch]$true
}

# Anatomic prep + HD-BET: needed by the 1H segmentation and by the -ExtRef V_M.
$runAnat = (-not $Skip1H) -or $ExtRef

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
if ($runAnat) { New-Item -ItemType Directory -Force -Path $stageDir | Out-Null }
New-Item -ItemType Directory -Force -Path $FinalOutDir | Out-Null

# ---- Stage 0: -ReuseAnat, seed staging from an earlier run ---------------
# Files already in staging (left by a failed run) are not overwritten.
if ($ReuseAnat -and $runAnat) {
    $prior  = @(Get-ChildItem $FinalOutDir -File | Where-Object { $_.Name -match '\.(nii|nii\.gz|tsv|txt)$' })
    $copied = 0
    foreach ($f in $prior) {
        $dst = Join-Path $stageDir $f.Name
        if (-not (Test-Path $dst)) { Copy-Item $f.FullName $dst; $copied++ }
    }
    Write-Host "-ReuseAnat: staged $copied anatomic file(s) from $FinalOutDir" -ForegroundColor Yellow
}

# ---- Stage 1: anatomic prep (MATLAB on Windows) ----------------------------
# Reads DICOMs from $SessionDir on Z:\, writes NIfTIs to $stageDir.
if ($runAnat) {
    Invoke-Matlab "addpath(genpath('$camipmDir')); addpath(genpath('$svsCodeDir')); prep_anatomic_U01_Brain_1H('$SessionDir','$stageDir'); rmpath(genpath('$camipmDir')); rmpath(genpath('$svsCodeDir'));"

    $t1Bias  = (Get-ChildItem $stageDir -Filter '*_SPM8BIAS.nii' | Select-Object -First 1).FullName
    $svsMask = Join-Path $stageDir 'svsmask.nii'
    if (-not $t1Bias)              { throw "Stage 1: no *_SPM8BIAS.nii in $stageDir" }
    if (-not (Test-Path $svsMask)) { throw "Stage 1: no svsmask.nii in $stageDir"   }
}

# ---- Stage 1b: HD-BET on the host --------------------------------------
# Writes <stem>_hdbet.nii.gz into $stageDir so the HD-BET block inside
# brainseg.sh sees it on disk and skips its own (slow) call.
# hd-bet 2.x requires .nii.gz input; our staged T1 is .nii, so gzip a
# sibling copy, run hd-bet, then delete it.
# We do not write _hdbet_mask: brainseg.sh only references the mask in
# its -U/-T branches, which we never trigger.
if ($runAnat) {
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
}

# ---- Stage 1c: FIRST without BrStem ------------------------------------
# Pre-staged like HD-BET, so brainseg.sh logs "Skipping FIRST." (NOTES 18).
# BrStem dropped; list must end on a boundary-corrected structure.
if (-not $Skip1H) {
    $firstStructs = 'L_Accu,L_Amyg,L_Caud,L_Hipp,L_Pall,L_Puta,L_Thal,' +
                    'R_Accu,R_Amyg,R_Caud,R_Pall,R_Puta,R_Thal,R_Hipp'

    $t1BiasVm   = Convert-ToVmPath $t1Bias
    $svsMaskVm  = Convert-ToVmPath $svsMask
    $firstIn    = ($t1BiasVm -replace '\.nii$','') + '_hdbet'
    $firstOut   = "${firstIn}_first"
    $firstDone  = Join-Path $stageDir "${t1BiasStem}_hdbet_first.nii.gz"

    if (Test-Path $firstDone) {
        Write-Host "FIRST output already exists, skipping: $firstDone" -ForegroundColor Yellow
    } else {
        Invoke-Vm "run_first_all -b -s '$firstStructs' -i '$firstIn' -o '$firstOut'"
        Invoke-Vm "immv '${firstOut}_all_fast_firstseg' '$firstOut'"
        Invoke-Vm "print_FIRST_labels.sh '${firstOut}.txt'"
        Invoke-Vm "rm -fr ${firstOut}*.vtk ${firstOut}*.bvars ${firstOut}.com* ${firstOut}.logs ${firstIn}_to_std_sub.mat ${firstOut}_all_fast_origsegs* ${firstOut}-L_* ${firstOut}-R_*"
    }
    if (-not (Test-Path $firstDone)) { throw "Stage 1c: run_first_all did not produce $firstDone" }
}

# ---- Stage 2: brainseg.sh in the VM ------------------------------------
# With _hdbet and _first outputs pre-staged, brainseg.sh logs "Skipping hd-bet."
# and "Skipping FIRST.", and proceeds to FAST -> SVS mask -> ROI stats.
# With -Skip1H, the seg TSV from an earlier full run in $FinalOutDir (if any).
if (-not $Skip1H) {
    Invoke-Vm "brainseg.sh -f -b N -S '$svsMaskVm' '$t1BiasVm' '$scanID'"

    $segStem = [IO.Path]::GetFileNameWithoutExtension($t1Bias)
    $segTsv  = Get-ChildItem $stageDir -Filter "${segStem}_hdbet_fast3_first_seg_svsmask.tsv" |
                  Select-Object -First 1 -ExpandProperty FullName
    if (-not $segTsv) { throw "Stage 2: brainseg.sh did not produce expected .tsv" }
} else {
    $segTsv  = Get-ChildItem $FinalOutDir -Filter '*_hdbet_fast3_first_seg_svsmask.tsv' |
                  Select-Object -First 1 -ExpandProperty FullName
}
if ($segTsv) {
    Write-Host "Seg TSV: $segTsv" -ForegroundColor Green
} else {
    Write-Host "No seg TSV in $FinalOutDir; 31P f_CSF (QC only) will be NaN." -ForegroundColor Yellow
}

# ---- Stage 3: 1H NAD + absolute quant ----------------------------------
# Reads *svssel*.dat from $SessionDir\datfiles, tissue fractions from $segTsv,
# writes <scanID>_brain_1H_abs.mat directly to $FinalOutDir.
if (-not $Skip1H) {
    Invoke-Matlab "addpath('$svsCodeDir'); processU01_Brain_1H('$SessionDir','$FinalOutDir','$segTsv',$pltArg); rmpath('$svsCodeDir');$pltWait"
}

# ---- Stage 4: 31P quant ------------------------------------------------
# Reads *31p*.dat from $SessionDir\datfiles, uses f_CSF from $segTsv, writes
# <scanID>_brain_31P_abs.mat directly to $FinalOutDir.
# 5th arg to processU01_Brain_31P: an extRef struct literal, or an empty
# struct() which leaves the feature off (enable defaults to false).
# With -ExtRef the external quant needs imscribe_engine (CAMIPM) to intersect
# the HD-BET brain mask with the ISIS slab for V_M.  Without it on the path
# that call throws, gets caught, and V_M silently falls back to a head mask
# that includes skull and scalp -- ~1.4x too large.  Added only for -ExtRef so
# the plain internal-reference run keeps its lean path.
if (-not $Skip31P) {
    if ($ExtRef) {
        if ($ApplyB1) { $b1Arg = 'true' } else { $b1Arg = 'false' }
        $extRefArg = "struct('enable',true,'refPpm',$RefPpm,'refConc_mM',$RefConc,'refArea_mm2',$RefArea,'applyB1',$b1Arg," +
                     "'brainMaskNii','$hdbetOut','t1Nii','$t1Bias')"
        $camipmAdd = "addpath(genpath('$camipmDir')); "
        $camipmRm  = " rmpath(genpath('$camipmDir'));"
    } else {
        $extRefArg = "struct()"
        $camipmAdd = ''
        $camipmRm  = ''
    }
    Invoke-Matlab "${camipmAdd}addpath('$svsCodeDir'); processU01_Brain_31P('$SessionDir','$FinalOutDir','$segTsv',$pltArg,$extRefArg,'$Ref31P'); rmpath('$svsCodeDir');$camipmRm$pltWait"
}

# ---- Stage 5: NAD summary to console -----------------------------------
# Reload the abs.mat file of each quant that ran and print NADH2/H4/H6/Trp
# (1H, uM) and NADplus/NADH (31P, uM) to the window.
$summaryStmt = "fprintf('\n=== NAD summary: ${scanID} ===\n');"
if (-not $Skip1H) {
    $mat1H = Join-Path $FinalOutDir "${scanID}_brain_1H_abs.mat"
    $summaryStmt += " m1=load('$mat1H','absQ');" +
                    " for k=1:numel(m1.absQ.metabolites), n=m1.absQ.metabolites(k).name;" +
                    "  if any(strcmp(n,{'NADH2','NADH4','NADH6','Trp'}))," +
                    "   fprintf('  1H  %-8s : %.4g uM\n',n,m1.absQ.metabolites(k).conc_mM*1e3); end; end;"
}
if (-not $Skip31P) {
    $mat31P = Join-Path $FinalOutDir "${scanID}_brain_31P_abs.mat"
    $summaryStmt += " m31=load('$mat31P','absQ');" +
                    " for k=1:numel(m31.absQ.metabolites), n=m31.absQ.metabolites(k).name;" +
                    "  if any(strcmp(n,{'NADplus','NADH'}))," +
                    "   fprintf('  31P %-8s : %.4g uM\n',n,m31.absQ.metabolites(k).conc_mM*1e3); end; end;"
}
Invoke-Matlab $summaryStmt

# ---- Stage final: ship staging outputs to Z: and clean up --------------
if ($runAnat) {
    Write-Host "Copying staging outputs to $FinalOutDir" -ForegroundColor Cyan
    Copy-Item -Path (Join-Path $stageDir '*') -Destination $FinalOutDir -Recurse -Force

    if (-not $KeepStaging) {
        Remove-Item -Recurse -Force $stageDir
    } else {
        Write-Host "Staging preserved at $stageDir" -ForegroundColor Yellow
    }
}

$pipelineEnd = Get-Date
$elapsed     = $pipelineEnd - $pipelineStart
Write-Host ("Pipeline end:   {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineEnd) -ForegroundColor Green
Write-Host ("Total run time: {0:hh\:mm\:ss} ({1:N1} min)" -f $elapsed, $elapsed.TotalMinutes) -ForegroundColor Green
Write-Host "Done: $FinalOutDir" -ForegroundColor Green
