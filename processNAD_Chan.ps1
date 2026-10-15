<#
processNAD_Chan.ps1 - NAD_Chan (858813) brain SVS pipeline driver.

Same architecture as processU01_Brain.ps1, 1H only, two scans per session:
  slab - *svssel*_sl_NAD.dat  / _sl_WAT.dat   (whole slab)
  vox  - *svssel*_vox_NAD.dat / _vox_WAT.dat  (occipital voxel)

  - DICOMs and .dat read directly by MATLAB on Windows.
  - Anatomic prep NIfTI outputs (T1, T1_SPM8BIAS, one svsmask per scan)
    staged into the VirtualBox shared folder so brainseg.sh in the VM can
    see them.
  - HD-BET runs on the Windows host; FIRST runs in the VM with the reduced
    structure list; brainseg.sh then runs once per mask off the same
    segmentation, giving one tissue-fraction TSV per scan.
  - MATLAB quant outputs (*.mat) written directly to <SessionDir>\processed\.
  - Final step: copy non-.mat outputs from staging back to processed,
    then remove staging dir.

Usage:
    .\processNAD_Chan.ps1 -SessionDir 'C:\Users\CAMIPM-NW\Documents\858813_NAD_smallVox\858813_S02_858813_S02'
    .\processNAD_Chan.ps1 -SessionDir '...\<session>' -Plot
    .\processNAD_Chan.ps1 -SessionDir '...\SXXX_SXXX_RZ' -CombineRuns

Note: Make sure VM is running when this program is run.
#>
param(
    [Parameter(Mandatory)][string]$SessionDir,
    [string]$StagingRoot = 'C:\Users\CAMIPM-NW\Documents\VirtualBoxSharedFolders\CAMIPM_Ubuntu_VM2\TEMP',
    [string]$FinalOutDir = (Join-Path $SessionDir 'processed'),
    [string]$ExamPattern = 'RESEARCH*',   # exam folder under $SessionDir holding the DICOM series
    [switch]$KeepStaging,
    [switch]$Plot,         # -Plot turns on pars.plt and pars.pltSpec
    [switch]$CombineRuns   # concatenate every water-suppressed run of a scan
                           #   instead of using only the highest-MID one
)

if ($Plot)        { $pltArg     = 'true' } else { $pltArg     = 'false' }
if ($CombineRuns) { $combineArg = 'true' } else { $combineArg = 'false' }

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

# scan -> DICOM series used for its imscribe mask
$Scans = [ordered]@{
    'slab' = 'SVSSEL*SL_NAD*'
    'vox'  = 'SVSSEL*VOX_NAD*'
}

# When plotting, append (per MATLAB call) a snippet that first SAVES every open
# figure to $FinalOutDir as .fig + .png, then waits for them to be closed before
# the function returns (so -batch MATLAB doesn't exit and destroy them).
function Get-PltWait([string]$tag) {
    if (-not $Plot) { return '' }
    return " drawnow;" +
        " figDir='$FinalOutDir'; sid='${scanID}_$tag';" +
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

# ---- Stage 1: anatomic prep (MATLAB on Windows) ------------------------
# T1 -> NIfTI -> SPM8 bias once; one imscribe mask per scan (svsmask_<scan>.nii).
$prepCalls = ($Scans.Keys | ForEach-Object {
    "prep_anatomic_U01_Brain_1H('$SessionDir','$stageDir','ExamPattern','$ExamPattern','SvsPattern','$($Scans[$_])','MaskName','svsmask_$_.nii');"
}) -join ' '
Invoke-Matlab "addpath(genpath('$camipmDir')); addpath(genpath('$svsCodeDir')); $prepCalls rmpath(genpath('$camipmDir')); rmpath(genpath('$svsCodeDir'));"

$t1Bias = (Get-ChildItem $stageDir -Filter '*_SPM8BIAS.nii' | Select-Object -First 1).FullName
if (-not $t1Bias) { throw "Stage 1: no *_SPM8BIAS.nii in $stageDir" }
foreach ($s in $Scans.Keys) {
    $m = Join-Path $stageDir "svsmask_$s.nii"
    if (-not (Test-Path $m)) { throw "Stage 1: no $m" }
}

# ---- Stage 1b: HD-BET on the host --------------------------------------
# Writes <stem>_hdbet.nii.gz into $stageDir so the HD-BET block inside
# brainseg.sh sees it on disk and skips its own call.
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
    & $hdbet -i $t1BiasGz -o $hdbetOut -device $hdbetDev
    $hdbetExit = $LASTEXITCODE
    Remove-Item -Force $t1BiasGz -ErrorAction SilentlyContinue
    if ($hdbetExit -ne 0) { throw "HD-BET failed (exit $hdbetExit)" }
}
if (-not (Test-Path $hdbetOut)) { throw "Stage 1b: HD-BET did not produce $hdbetOut" }

# ---- Stage 1c: FIRST without BrStem ------------------------------------
# Same pre-staging as processU01_Brain.ps1 (processU01_Brain_1H_NOTES.md 18).
$firstStructs = 'L_Accu,L_Amyg,L_Caud,L_Hipp,L_Pall,L_Puta,L_Thal,' +
                'R_Accu,R_Amyg,R_Caud,R_Pall,R_Puta,R_Thal,R_Hipp'

$t1BiasVm   = Convert-ToVmPath $t1Bias
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

# ---- Stage 2: brainseg.sh in the VM, once per mask ---------------------
# brainseg.sh names its mask outputs <seg>_svsmask.* whatever the mask file is
# called, so each scan's set is renamed to <seg>_svsmask_<scan>.* before the next.
$segStem = "${t1BiasStem}_hdbet_fast3_first_seg"
$segTsv  = @{}
foreach ($s in $Scans.Keys) {
    Get-ChildItem $stageDir -Filter "${segStem}_svsmask.*" | Remove-Item -Force
    $maskVm = Convert-ToVmPath (Join-Path $stageDir "svsmask_$s.nii")
    Invoke-Vm "brainseg.sh -f -b N -S '$maskVm' '$t1BiasVm' '$scanID'"

    foreach ($ext in 'nii.gz','tsv','txt') {
        $src = Join-Path $stageDir "${segStem}_svsmask.$ext"
        if (-not (Test-Path $src)) { throw "Stage 2: brainseg.sh did not produce $src" }
        Move-Item -Force $src (Join-Path $stageDir "${segStem}_svsmask_$s.$ext")
    }
    $segTsv[$s] = Join-Path $stageDir "${segStem}_svsmask_$s.tsv"
    Write-Host "Seg TSV ($s): $($segTsv[$s])" -ForegroundColor Green
}

# ---- Stage 3: 1H NAD + absolute quant, per scan ------------------------
# Writes <scanID>_<scan>_1H_abs.mat directly to $FinalOutDir.
foreach ($s in $Scans.Keys) {
    Invoke-Matlab "addpath('$svsCodeDir'); processNAD_Chan_1H('$SessionDir','$FinalOutDir','$($segTsv[$s])','$s',$pltArg,'',$combineArg); rmpath('$svsCodeDir');$(Get-PltWait $s)"
}

# ---- Stage 3b: fit figures ---------------------------------------------
# Writes <scanID>_<scan>_1H_fit.png (data, fit, residual per region) to $FinalOutDir.
$plotCalls = ($Scans.Keys | ForEach-Object { "plotU01_2026_fits('1H','$_','$SessionDir',false);" }) -join ' '
Invoke-Matlab "addpath('$svsCodeDir'); $plotCalls rmpath('$svsCodeDir');"

# ---- Stage 4: NAD summary to console -----------------------------------
$summaryStmt = " fprintf('\n=== NAD summary: ${scanID} ===\n');" +
               " scans={" + (($Scans.Keys | ForEach-Object { "'$_'" }) -join ',') + "};" +
               " for s=1:numel(scans), m=load(fullfile('$FinalOutDir',sprintf('%s_%s_1H_abs.mat','$scanID',scans{s})),'absQ');" +
               "  for k=1:numel(m.absQ.metabolites), n=m.absQ.metabolites(k).name;" +
               "   fprintf('  %-5s %-8s : %.4g uM\n',scans{s},n,m.absQ.metabolites(k).conc_mM*1e3); end; end"
Invoke-Matlab $summaryStmt

# ---- Stage 5: ship staging outputs and clean up ------------------------
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
