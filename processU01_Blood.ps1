<#
processU01_Blood.ps1 - U01 blood (whole-blood SVS on Bruker high-field NMR)
                       pipeline driver.

Architecture:
  - Blood is a single compartment with no anatomic segmentation, so none of
    the brain pipeline's staging / HD-BET / VirtualBox / brainseg.sh machinery
    is needed.  The driver just runs the MATLAB quant function and prints
    a summary.
  - Bruker data are read directly by MATLAB on Windows.  Each acquisition is
    its own labeled Bruker experiment subfolder under <SessionDir>,
    containing fid/ser + acqus:
        selse_4p7ppm_T2relax  water T2 spin-echo series (run first)
        selse_9p5ppm          downfield NAD/Trp
        selse_4p7ppm          water reference
        selse_0ppm            TSP reference
    (31p is present but not processed yet.)
  - MATLAB quant outputs (*.mat, *.tsv) are written directly to <FinalOutDir>
    (default <SessionDir>\processed\).

Usage:
    .\processU01_Blood.ps1 -SessionDir 'C:\...\<scanID>'
    .\processU01_Blood.ps1 -SessionDir 'C:\...\<scanID>' -Plot
#>
param(
    [Parameter(Mandatory)][string]$SessionDir,
    [string]$FinalOutDir = (Join-Path $SessionDir 'processed'),
    [switch]$Plot          # -Plot turns on pars.plt (and pars.pltSpec for 1H)
)

if ($Plot) { $pltArg = 'true' } else { $pltArg = 'false' }

$ErrorActionPreference = 'Stop'

$pipelineStart = Get-Date
Write-Host ("Pipeline start: {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineStart) -ForegroundColor Green

# ---- Configuration -----------------------------------------------------
$matlab     = 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe'
$scanID     = Split-Path $SessionDir -Leaf
$stem       = if ($scanID -match '_blood$') { $scanID } else { "${scanID}_blood" }   # matches utils/bloodOutStem.m
$svsCodeDir = $PSScriptRoot   # this script lives in svs_code/

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

# ---- Helpers -----------------------------------------------------------
function Invoke-Matlab([string]$Stmt) {
    Write-Host ">> MATLAB: $Stmt" -ForegroundColor Cyan
    & $matlab -batch $Stmt -wait
    if ($LASTEXITCODE -ne 0) { throw "MATLAB failed (exit $LASTEXITCODE)" }
}

# ---- Setup -------------------------------------------------------------
New-Item -ItemType Directory -Force -Path $FinalOutDir | Out-Null

foreach ($e in @('selse_4p7ppm_T2relax', 'selse_9p5ppm', 'selse_4p7ppm', 'selse_0ppm')) {
    $d = Join-Path $SessionDir $e
    if (-not (Test-Path $d)) { throw "Scan folder not found: $d" }
}

# ---- Stage 1: water T2 from the spin-echo series -----------------------
# Writes <scanID>_blood_waterT2.mat/.tsv; Stage 2 reads the .mat.
Invoke-Matlab "addpath('$svsCodeDir'); processU01_Blood_waterT2('$SessionDir','$FinalOutDir',$pltArg); rmpath('$svsCodeDir');$pltWait"

# ---- Stage 2: 1H NAD/Trp + absolute quant (water and TSP references) ---
Invoke-Matlab "addpath('$svsCodeDir'); processU01_Blood_1H('$SessionDir','$FinalOutDir',$pltArg); rmpath('$svsCodeDir');$pltWait"

# ---- Stage 3: summary to console ---------------------------------------
$mat1H = Join-Path $FinalOutDir "${stem}_1H_abs.mat"
$summaryStmt = "m1=load('$mat1H','absQ','absQtsp','waterCheck');" +
               " fprintf('\n=== Blood 1H summary: ${scanID} ===\n');" +
               " fprintf('  water T2 %.2f +/- %.2f ms (R^2 %.4f); ref/T2-series area ratio %.3f\n'," +
               "  m1.waterCheck.T2_ms,m1.waterCheck.T2_ci_ms,m1.waterCheck.R2,m1.waterCheck.refOverPred);" +
               " fprintf('  %-8s %12s %12s\n','peak','water [uM]','TSP [uM]');" +
               " for k=1:numel(m1.absQ.metabolites), n=m1.absQ.metabolites(k).name;" +
               "  if any(strcmp(n,{'NADH2','NADH6','Trp'}))," +
               "   fprintf('  %-8s %12.4g %12.4g\n',n,m1.absQ.metabolites(k).conc_mM*1e3,m1.absQtsp.metabolites(k).conc_mM*1e3); end; end"
Invoke-Matlab $summaryStmt

$pipelineEnd = Get-Date
$elapsed     = $pipelineEnd - $pipelineStart
Write-Host ("Pipeline end:   {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineEnd) -ForegroundColor Green
Write-Host ("Total run time: {0:hh\:mm\:ss} ({1:N1} min)" -f $elapsed, $elapsed.TotalMinutes) -ForegroundColor Green
Write-Host "Done: $FinalOutDir" -ForegroundColor Green
