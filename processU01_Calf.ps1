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
    .\processU01_Calf.ps1 -SessionDir '...\<scanID>' -Ref31P PCr
    .\processU01_Calf.ps1 -SessionDir '...\<scanID>' -Ref31P External -RefConc 500
    .\processU01_Calf.ps1 -SessionDir '...\<scanID>' -Skip31P    # 1H only
    .\processU01_Calf.ps1 -SessionDir '...\<scanID>' -Skip1H     # 31P only
#>
param(
    [Parameter(Mandatory)][string]$SessionDir,
    [string]$FinalOutDir = (Join-Path $SessionDir 'processed'),
    [switch]$Plot,         # -Plot turns on pars.plt (and pars.pltSpec for 1H)
    [switch]$Skip1H,       # -Skip1H skips the 1H quant
    [switch]$Skip31P,      # -Skip31P skips the 31P quant
    [switch]$ExtRef,      # -ExtRef adds the ~8 ppm phantom singlet to the 31P basis
                           #   and quantifies against it (in addition to PCr).
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

# 4th arg to processU01_Calf_31P: an extRef struct literal, or an empty
# struct() which leaves the feature off (enable defaults to false).
if ($ExtRef) {
    if ($ApplyB1) { $b1Arg = 'true' } else { $b1Arg = 'false' }
    $extRefArg = "struct('enable',true,'refPpm',$RefPpm,'refConc_mM',$RefConc,'refArea_mm2',$RefArea,'applyB1',$b1Arg)"
} else {
    $extRefArg = "struct()"
}

$ErrorActionPreference = 'Stop'

$pipelineStart = Get-Date
Write-Host ("Pipeline start: {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineStart) -ForegroundColor Green

# ---- Configuration -----------------------------------------------------
$matlab     = 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe'
$scanID     = Split-Path $SessionDir -Leaf
$svsCodeDir = $PSScriptRoot   # this script lives in svs_code/
$camipmDir  = Join-Path (Split-Path $svsCodeDir -Parent) 'CAMIPM'          # sibling of svs_code/
$dixonDir   = Join-Path (Split-Path $svsCodeDir -Parent) 'dixonProcessing' # sibling of svs_code/

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
if (-not $Skip1H) {
    Invoke-Matlab "addpath('$svsCodeDir'); processU01_Calf_1H('$SessionDir','$FinalOutDir',$pltArg); rmpath('$svsCodeDir');$pltWait"
}

# ---- Stage 2: 31P quant ------------------------------------------------
# Reads *31p*.dat from $SessionDir\datfiles, PCr internal reference, writes
# <scanID>_calf_31P_abs.mat directly to $FinalOutDir.
# With -ExtRef the external quant needs imscribe_engine + dicom2nifti (CAMIPM)
# and readB0FWgre (dixonProcessing) to turn the Dixon into a muscle mask and
# intersect it with the ISIS slab for V_M.  Without them on the path that call
# throws, gets caught, and V_M silently falls back to a B1 threshold leg mask
# that includes fat and marrow and covers only part of the slab.  Added only
# for -ExtRef so the plain internal-reference run keeps its lean path.
if (-not $Skip31P) {
    if ($ExtRef) {
        # dixonProcessing goes on with '-end': its bundled fwtoolbox ships
        # matlab_bgl/test/assert.m, which would otherwise shadow the built-in
        # assert for every function in this MATLAB session.
        $extraAdd = "addpath(genpath('$camipmDir')); addpath(genpath('$dixonDir'),'-end'); "
        $extraRm  = " rmpath(genpath('$camipmDir')); rmpath(genpath('$dixonDir'));"
    } else {
        $extraAdd = ''
        $extraRm  = ''
    }
    Invoke-Matlab "${extraAdd}addpath('$svsCodeDir'); processU01_Calf_31P('$SessionDir','$FinalOutDir',$pltArg,$extRefArg,'$Ref31P'); rmpath('$svsCodeDir');$extraRm$pltWait"
}

# ---- Stage 3: NAD summary to console -----------------------------------
# Reload the abs.mat file of each quant that ran and print NADH2/H4/H6
# (1H, uM) and NADplus/NADH (31P, uM) to the window.
$summaryStmt = "fprintf('\n=== NAD summary: ${scanID} ===\n');"
if (-not $Skip1H) {
    $mat1H = Join-Path $FinalOutDir "${scanID}_calf_1H_abs.mat"
    $summaryStmt += " m1=load('$mat1H','absQ');" +
                    " for k=1:numel(m1.absQ.metabolites), n=m1.absQ.metabolites(k).name;" +
                    "  if any(strcmp(n,{'NADH2','NADH4','NADH6'}))," +
                    "   fprintf('  1H  %-8s : %.4g uM\n',n,m1.absQ.metabolites(k).conc_mM*1e3); end; end;"
}
if (-not $Skip31P) {
    $mat31P = Join-Path $FinalOutDir "${scanID}_calf_31P_abs.mat"
    $summaryStmt += " m31=load('$mat31P','absQ');" +
                    " for k=1:numel(m31.absQ.metabolites), n=m31.absQ.metabolites(k).name;" +
                    "  if any(strcmp(n,{'NADplus','NADH'}))," +
                    "   fprintf('  31P %-8s : %.4g uM\n',n,m31.absQ.metabolites(k).conc_mM*1e3); end; end;"
}
Invoke-Matlab $summaryStmt

$pipelineEnd = Get-Date
$elapsed     = $pipelineEnd - $pipelineStart
Write-Host ("Pipeline end:   {0:yyyy-MM-dd HH:mm:ss}" -f $pipelineEnd) -ForegroundColor Green
Write-Host ("Total run time: {0:hh\:mm\:ss} ({1:N1} min)" -f $elapsed, $elapsed.TotalMinutes) -ForegroundColor Green
Write-Host "Done: $FinalOutDir" -ForegroundColor Green
