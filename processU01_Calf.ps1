<#
processU01_Calf.ps1 - U01 calf (skeletal muscle) SVS pipeline driver.

Architecture:
  - Calf is a single muscle compartment: there is NO anatomic segmentation,
    so none of the brain pipeline's staging / HD-BET / VirtualBox /
    brainseg.sh machinery is needed.  The driver just runs the two MATLAB
    quant functions and prints a summary.
  - .dat files are read directly by MATLAB on Windows from
    <SessionDir>\datfiles\.
  - MATLAB quant outputs (*.mat) are written directly to <FinalOutDir>
    (default <SessionDir>\processed\).

Usage:
    .\processU01_Calf.ps1 -SessionDir 'Z:\Data\Spectroscopy\...\<scanID>'
    .\processU01_Calf.ps1 -SessionDir '...\<scanID>' -Plot
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

# ---- Stage 1: 1H NAD + absolute quant ----------------------------------
# Reads *svssel*.dat from $SessionDir\datfiles, single muscle compartment,
# writes <scanID>_calf_1H_abs.mat directly to $FinalOutDir.
Invoke-Matlab "addpath('$svsCodeDir'); processU01_Calf_1H('$SessionDir','$FinalOutDir',$pltArg); rmpath('$svsCodeDir');$pltWait"

# ---- Stage 2: 31P quant ------------------------------------------------
# Reads *31p*.dat from $SessionDir\datfiles, PCr internal reference, writes
# <scanID>_calf_31P_abs.mat directly to $FinalOutDir.
Invoke-Matlab "addpath('$svsCodeDir'); processU01_Calf_31P('$SessionDir','$FinalOutDir',$pltArg); rmpath('$svsCodeDir');$pltWait"

# ---- Stage 3: NAD summary to console -----------------------------------
# Reload the two abs.mat files and print NADH2/H4/H6 (1H, uM) and
# NADplus/NADH (31P, mM) to the window.
$mat1H  = Join-Path $FinalOutDir "${scanID}_calf_1H_abs.mat"
$mat31P = Join-Path $FinalOutDir "${scanID}_calf_31P_abs.mat"
$summaryStmt = "m1=load('$mat1H','absQ'); m31=load('$mat31P','absQ');" +
               " fprintf('\n=== NAD summary: ${scanID} ===\n');" +
               " for k=1:numel(m1.absQ.metabolites), n=m1.absQ.metabolites(k).name;" +
               "  if any(strcmp(n,{'NADH2','NADH4','NADH6'}))," +
               "   fprintf('  1H  %-8s : %.4g uM\n',n,m1.absQ.metabolites(k).conc_mM*1e3); end; end;" +
               " for k=1:numel(m31.absQ.metabolites), n=m31.absQ.metabolites(k).name;" +
               "  if any(strcmp(n,{'NADplus','NADH'}))," +
               "   fprintf('  31P %-8s : %.4g uM\n',n,m31.absQ.metabolites(k).conc_mM*1e3); end; end"
Invoke-Matlab $summaryStmt

$pipelineEnd = Get-Date
$elapsed     = $pipelineEnd - $pipelineStart
Write-Host ("Pipeline end:   {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineEnd) -ForegroundColor Green
Write-Host ("Total run time: {0:hh\:mm\:ss} ({1:N1} min)" -f $elapsed, $elapsed.TotalMinutes) -ForegroundColor Green
Write-Host "Done: $FinalOutDir" -ForegroundColor Green
