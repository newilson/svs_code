<#
summarizeNAD_Chan.ps1 - NAD_Chan cohort table.

Collects <session>\processed\<session>_<scan>_1H_abs.mat (scan = slab, vox)
from every session folder under -DataRoot into one TSV via
collectNAD_Chan_TSV.m.  Run after processNAD_Chan.ps1.

Usage:
    .\summarizeNAD_Chan.ps1
    .\summarizeNAD_Chan.ps1 -DataRoot 'C:\...\858813_NAD_smallVox' -OutPath 'C:\...\NAD_Chan_1H.tsv'
#>
param(
    [string]$DataRoot = 'C:\Users\CAMIPM-NW\Documents\858813_NAD_smallVox',
    [string]$OutPath  = (Join-Path $DataRoot 'NAD_Chan_1H.tsv')
)

$ErrorActionPreference = 'Stop'

$matlab     = 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe'
$svsCodeDir = $PSScriptRoot

$stmt = "addpath('$svsCodeDir'); collectNAD_Chan_TSV('$DataRoot','$OutPath'); rmpath('$svsCodeDir');"
Write-Host ">> MATLAB: $stmt" -ForegroundColor Cyan
& $matlab -batch $stmt -wait
if ($LASTEXITCODE -ne 0) { throw "MATLAB failed (exit $LASTEXITCODE)" }
Write-Host "Done: $OutPath" -ForegroundColor Green
