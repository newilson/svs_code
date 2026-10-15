<#
testBrainSeg.ps1 - Compare brainseg results across multiple T1w acquisitions.

For each GRE3D* and MPRAGE* DICOM series under SessionDir:
  1) DICOM -> NIfTI and SPM bias correction (MATLAB on Windows).
  2) brainseg.sh in the VM, using svsmask_1mm.nii or svsmask_2mm.nii
     based on the bias-corrected scan's voxel resolution.
  3) Compare GM/WM/CSF fractions inside the SVS voxel; rank scans by
     similarity to MPRAGE_P2* (ground truth).

Usage:
    .\testBrainSeg.ps1
    .\testBrainSeg.ps1 -ForceRerun

Note: VM must be running. Reuses the SSH alias 'ubuntu-vm' and the
shared-folder mapping under C:\...\CAMIPM_Ubuntu_VM2\.
#>
param(
    [string]$SessionDir = 'C:\Users\CAMIPM-NW\Documents\VirtualBoxSharedFolders\CAMIPM_Ubuntu_VM2\Data\Anatomic\S000_S000\RESEARCH_RAVIN_20260423_151710_226000',
    [string]$OutDir     = (Join-Path 'C:\Users\CAMIPM-NW\Documents\VirtualBoxSharedFolders\CAMIPM_Ubuntu_VM2\Data\Anatomic\S000_S000\RESEARCH_RAVIN_20260423_151710_226000' 'niftis'),
    [switch]$ForceRerun
)

$ErrorActionPreference = 'Stop'

# ---- Configuration -----------------------------------------------------
$matlab = 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe'
$vmHost = 'ubuntu-vm'

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
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

# ---- Stage 1: dicom2nifti + bias correction (MATLAB) -------------------
$forceArg = if ($ForceRerun) { 'true' } else { 'false' }
Invoke-Matlab "addpath(genpath('Z:\code\CAMIPM')); addpath('Z:\code\svs_code'); prep_anat_test('$SessionDir','$OutDir',$forceArg);"

$manifestPath = Join-Path $OutDir 'manifest.csv'
if (-not (Test-Path $manifestPath)) { throw "Stage 1: manifest.csv not produced in $OutDir" }
$manifest = Import-Csv $manifestPath

# ---- Stage 2: brainseg.sh per scan in the VM ---------------------------
foreach ($row in $manifest) {
    $segDone = Join-Path $OutDir "$($row.scanName)_SPM8bias_hdbet_fast3_first_seg_svsmask.tsv"
    if ((Test-Path $segDone) -and (-not $ForceRerun)) {
        Write-Host "[skip] brainseg already done for $($row.scanName)" -ForegroundColor DarkGray
        continue
    }
    $t1Vm   = Convert-ToVmPath $row.biasNii
    $maskVm = Convert-ToVmPath $row.svsMask
    $force  = if ($ForceRerun) { '-x ' } else { '' }
    Invoke-Vm "brainseg.sh $force-f -b N -S '$maskVm' '$t1Vm' '$($row.scanName)'"
}

# ---- Stage 3: compare SVS voxel composition (MATLAB) -------------------
Invoke-Matlab "addpath('Z:\code\svs_code'); compare_svs_segs('$OutDir');"

Write-Host "Done. Comparison written to $(Join-Path $OutDir 'svs_seg_comparison.csv')" -ForegroundColor Green
