# processNAD_Chan — design notes

NAD_Chan (858813) brain 1H NAD+ project. Driver `processNAD_Chan.ps1`,
per-scan quant `processNAD_Chan_1H.m`. Same pipeline architecture as
`processU01_Brain.ps1` (see `processU01_Brain_1H_NOTES.md`), 1H only.

## 1. Data layout and protocol

`C:\Users\CAMIPM-NW\Documents\858813_NAD_smallVox\<session>\RESEARCH_NWILSON_*\`
holds both the DICOM series and the .dat files (no `datfiles\`).

| scan | .dat | geometry | TE / TR | avgs |
|---|---|---|---|---|
| slab | `*svssel_latest_sl_NAD/_WAT.dat` | 200 x 200 x 40 mm, slice-selective | 13 / 500 ms | 512 |
| vox | `*svssel_latest_vox_NAD/_WAT.dat` | ~24 x 24 x 36 mm occipital, 3D voxel | 26 / 500 ms | 2048 |

Both: 32-channel receive, dwell 250 us (4 kHz with oversampling, 1024 pts —
twice U01's spectral resolution), water reference 8 avgs, TR 10 s, same TE
as its metabolite scan. `wiplong(1)` is 1 for the slab and 0 for the voxel:
the voxel is a genuine 3D localization, so non-brain volume inside its mask
is real partial volume (unlike the U01 slab; see U01 notes on slice
localization).

Sessions (Sep 2026): 858813_S02, 858813_S30, SXX_SXX_RR, SXXX_SXXX_RZ.
MATHUR_25_09_12 in the same folder is slab-only with an MP2RAGE — not used.
The `svsglu` upfield scans are not processed yet.

SXXX_SXXX_RZ has two 1024-avg voxel runs (MID01184/01185, identical
headers). Default file selection is the highest-MID run; `-CombineRuns`
concatenates the averages of every run (`dat2mat_svsNAD` accepts a cell of
runs and refuses mismatched headers). Checked: model-free NADH2/water is
6.11e-6 (run 1), 6.55e-6 (run 2), 6.26e-6 combined.

## 2. Segmentation: one FAST/FIRST, one mask per scan

`prep_anatomic_U01_Brain_1H` is called once per scan with
`ExamPattern 'RESEARCH*'`, `SvsPattern 'SVSSEL*SL_NAD*' / 'SVSSEL*VOX_NAD*'`
and `MaskName svsmask_slab.nii / svsmask_vox.nii`. HD-BET and FIRST
(reduced list, U01 notes §18) run once; `brainseg.sh -S` then runs once per
mask and skips everything already on disk except the mask step.

`brainseg.sh` names its mask outputs `<seg>_svsmask.*` regardless of the mask
file name, so the driver removes any existing `<seg>_svsmask.*` before each
call and renames the result to `<seg>_svsmask_<scan>.*` afterwards.

## 3. Launching the driver from a script

Under Windows PowerShell 5.1 with `$ErrorActionPreference = 'Stop'`,
redirecting the driver's streams (`& .\processNAD_Chan.ps1 ... *> log`)
turns the first stderr line of a native command into a terminating error.
HD-BET was killed mid-prediction that way, and its orphaned multiprocessing
worker kept the pipe open so the script hung. Run it from a console, or as a
child process with file-level redirection:

```powershell
Start-Process powershell.exe -ArgumentList '-NoProfile','-File','processNAD_Chan.ps1','-SessionDir',$s `
    -RedirectStandardOutput out.log -RedirectStandardError err.log -NoNewWindow -Wait
```

## 4. Fit settings

Slab: the U01 brain settings unchanged (fits as clean as U01; slab H6/H2
1.24 / 1.13 / 1.16 / 1.11 for S02 / S30 / SXX / SXXX).

Voxel (settled below: Lorentzian +-20, single region). NADH2 SNR ~15 (SXXX ~11), water FWHM ~16.5 Hz vs
23-26 Hz in the slab, NADH2 ~24 Hz vs 32-34 Hz. H6/H2 across the four
subjects (self-consistency only, not accuracy):

| voxel variant | S02 | S30 | SXX | SXXX | SD |
|---|---|---|---|---|---|
| U01 settings | 1.43 | 0.69 | 1.14 | 1.22 | 0.31 |
| complex Lorentzian (mode 3) | 1.16 | 0.69 | 1.15 | 0.97 | 0.22 |
| lambdaWidth 0.3 | 1.42 | 0.97 | 1.03 | 1.26 | 0.21 |
| ph_range +-20 | 1.90 | 0.86 | 1.07 | 1.10 | 0.46 |
| Lorentzian + ph_range +-20 | 1.16 | 0.80 | 1.01 | 0.95 | 0.15 |

Tried on S02/S30 only, no effect: block_align.size 256, region-1
knotSpacing 2, Trp in the width prior, lambdaWidth 0 (also destabilizes the
slab: H6/H2 up to 3.1).

When mode 3 is used for the metabolites, `pars.watfit.mode` must be set to 5
(and `watfit.ph_range` to [-60 60]) or the water fit silently follows
`pars.fit.mode`; the voxel water is Voigt to within 1% of peak.

**Adopted (2026-09-23): Lorentzian + ph_range +-20 for the voxel**, water
kept Voigt +-60.

### Line broadening and single region (Lorentzian +-20 base, 2026-09-23)

`lb_met` 3 / 6 / 10 Hz x layout {two regions (current); one region
[8.7 10.6] knotSpacing 0.43 (7 coef); one region knotSpacing 0.63 (6 coef)}.
"Trp railed" = Trp center on its 10.05/10.15 bound or width on 10/120 Hz.
SXXX r1 vs r2 = NADH2/water difference between the two 1024-avg runs fitted
separately (one same-session pair — noise-level evidence only).

| layout | lb | H6/H2 S02 S30 SXX SXXX | SD | Trp railed | NADH2 uM | SXXX r1 vs r2 |
|---|---|---|---|---|---|---|
| 2 regions | 3 | 1.16 0.80 1.01 0.95 | 0.15 | 3/4 | 234 184 286 234 | +21% |
| 2 regions | 6 | 1.22 0.76 1.16 0.92 | 0.22 | 2/4 | 207 257 228 237 | +12% |
| 2 regions | 10 | 1.38 0.83 1.24 1.04 | 0.24 | 0/4 | 265 271 267 201 | +81% |
| 1 region, 7 coef | 3 | 1.08 0.81 1.00 0.93 | 0.11 | 0/4 | 239 183 290 238 | +17% |
| 1 region, 7 coef | 6 | 1.26 0.75 1.29 0.96 | 0.26 | 0/4 | 208 255 222 229 | +9% |
| 1 region, 7 coef | 10 | 1.39 0.92 1.26 1.10 | 0.20 | 0/4 | 266 254 260 205 | +91% |
| 1 region, 6 coef | 3 | 1.07 0.76 0.98 0.92 | 0.13 | 0/4 | 235 193 289 238 | +21% |
| 1 region, 6 coef | 6 | 1.17 0.77 1.17 0.94 | 0.19 | 0/4 | 209 253 216 225 | +3% |
| 1 region, 6 coef | 10 | 1.43 0.91 1.28 1.09 | 0.22 | 1/4 | 262 256 249 195 | +81% |

Extra line broadening does not help: H6/H2 spread grows, NADH2 moves with
the filter (S30 184 -> 257 -> 271 uM, i.e. the result depends on lb), and at
10 Hz the SXXX runs disagree by 80-90% (the HSVD coarse phase shifts with lb,
e.g. 40 deg on one SXXX run). Residual falls with lb only because the data
are smoother. Keep lb_met = 3.

**Adopted (2026-09-24): single region [8.7 10.6], knotSpacing 0.43,
lb_met 3** for the voxel; the driver reproduces the trial fits exactly.

A single region stops Trp railing (0/4 vs 3/4), tightens H6/H2 slightly and
leaves NADH2 essentially unchanged; the baseline is smooth and monotone over
the whole window. Voxel Trp/H2 is 0.68-0.86 in S02/SXX; an earlier reading
of that as "3-4x high" is withdrawn in §6.

Open: in all four voxels the NADH2 peak is sharper than the fitted line
(residual spike at 9.38 ppm, data 10-20% above the fit) under every variant,
while the water line fits cleanly. Voxel Trp is at the noise level; with two
regions it railed on a bound in 3/4 subjects (none with the single region).

## 5. Outputs

Per session (`<session>\processed\`): `<session>_<scan>_1H_abs.mat`,
`<session>_<scan>_1H_fit.png` (Stage 3b, `plotU01_2026_fits`, data + fit +
residual per fit region), the per-scan seg TSV/mask, and the anatomic
intermediates.

Cohort (`<parent>\NAD_Chan_1H.tsv`, written by `summarizeNAD_Chan.ps1` ->
`collectNAD_Chan_TSV`, not by the per-session driver): one row
per session x scan with the four concentrations, tissue fractions, CSF / GM /
WM / GM+WM volumes inside the mask, mask volume, SNR_NADH2 and water FWHM.

SNR is not computed by the processing pipeline itself. `SNR_NADH2` in the
TSV = fitted NADH2 peak height / std of the real spectrum over 11.5-13.5 ppm
after removing a linear trend (after the pipeline's 3 Hz lb). The
`characterize_data_differences.m` noise window (1-3 ppm) is outside this
protocol's 4 kHz window (~3-16.4 ppm).

First pass, Sep 2026 (voxel with U01 settings; superseded by the Lorentzian
+-20 rerun, see the TSV):

| | NADH2 | NADH6 | NADH4 | Trp (uM) | GM+WM mL | SNR |
|---|---|---|---|---|---|---|
| S02 slab | 298 | 394 | 390 | 79 | 351 | 83 |
| S02 vox | 184 | 299 | 288 | 264 | 17.5 | 14 |
| S30 slab | 326 | 393 | 287 | 80 | 338 | 97 |
| S30 vox | 130 | 102 | 234 | 43 | 20.6 | 13 |
| SXX slab | 302 | 374 | 208 | 88 | 334 | 118 |
| SXX vox | 260 | 337 | 305 | 247 | 18.0 | 19 |
| SXXX slab | 357 | 425 | 245 | 72 | 279 | 93 |
| SXXX vox | 167 | 232 | 190 | 55 | 17.5 | 9 |

CSF + GM + WM equals the mask volume in all four voxels (20.7 mL; 24.3 mL
for S30), i.e. no non-brain volume inside the voxel for these subjects.
Voxel Trp at TE 26 carries a 3.9x T2 correction (T2 19 ms) on a peak at the
noise level; the 247-264 uM values are not credible.

Current results (2026-09-24; voxel = Lorentzian +-20, single region):

| | NADH2 | NADH6 | NADH4 | Trp (uM) | GM+WM mL | SNR |
|---|---|---|---|---|---|---|
| S02 slab | 298 | 394 | 390 | 79 | 351.1 | 83 |
| S02 vox | 239 | 293 | 299 | 341 | 17.5 | 14 |
| S30 slab | 326 | 393 | 287 | 80 | 338.3 | 97 |
| S30 vox | 183 | 169 | 254 | 102 | 20.6 | 13 |
| SXX slab | 302 | 374 | 208 | 88 | 333.7 | 118 |
| SXX vox | 290 | 328 | 311 | 326 | 18.0 | 19 |
| SXXX slab | 357 | 425 | 245 | 72 | 278.9 | 93 |
| SXXX vox | 238 | 251 | 200 | 100 | 17.5 | 9 |

Voxel Trp is ~100 uM in S30/SXXX (slab 72-88) but 326-341 uM in S02/SXX.
The NADH2 sharp-top residual at 9.38 ppm persists, largest in S30.

## 6. Trp linewidth and coil combination (2026-09-24)

Fitted Trp vs NADH2 FWHM (Hz). Slab (Voigt) Trp is Gaussian (eta ~1) and
broader than NAD in all four subjects; the voxel forces Trp Lorentzian.

| | slab H2 | slab Trp | vox H2 | vox Trp | vox Trp/H2 |
|---|---|---|---|---|---|
| S02 | 34.2 | 42.6 | 24.9 | 73.2 | 0.86 |
| S30 | 32.1 | 39.9 | 24.7 | 36.7 | 0.34 |
| SXX | 30.3 | 41.4 | 25.9 | 46.3 | 0.68 |
| SXXX | 36.7 | 41.5 | 22.0 | 14.2 | 0.25 |

SXXX voxel Trp is narrower than its NAD and does not reproduce between the
two runs (run 1: 11.6 Hz, Trp/H2 0.21; run 2: 99.6 Hz, Trp/H2 1.88) — noise.
The earlier reading that S02/SXX voxel Trp was "3-4x high" rested only on a
concentration prior (slab Trp/H2 scaled with a literature T2 of 19 ms) and is
withdrawn. Also, a Lorentzian of the same FWHM carries ~1.48x the area of a
Gaussian, with tails spread across the 1.9 ppm window, so voxel and slab Trp
are not fitted with comparable lineshapes.

Coil combination (`coilCombinationNoPC`): weights |s_i|^2 from the water
reference first point, 5% floor, no noise weighting or prewhitening.
Cross-validated on the four voxels (noise covariance from the 11.5-13.5 ppm
band of the odd averages, SNR measured on the even-average spectrum):
current / optimal (C^-1 s) = 0.88-0.90, noise-weighted MRC / optimal =
0.90-0.94. A noise covariance taken from the FID tail instead predicted
0.38-0.60 and was wrong — the empirical spectrum SNR reversed it; the tail
carries a coil-correlated component that does not appear in the metabolite
band.

## 7. Voigt with the Gaussian width fixed from water (2026-09-24) — tried and rejected

Rationale: water sees the same field distribution as the metabolites and is
~half Gaussian in both slab (eta 0.44-0.55) and voxel (0.49-0.57); slab NAD
eta ~0 is the ill-determined direction (U01 notes §14), not physics. Tested
with a temporary option that fitted the 1-component Voigt water reference
before the metabolites and fixed each NAD peak's Gaussian FWHM to water's
(7.9-9.6 Hz in the voxels), with eta = wG/wV following from the total width.
**The option was deleted after this test** (per user); S02 slab + voxel
outputs are bit-identical with it removed.

Voxel test (single region, ph +-20, NAD from water, Trp free Voigt). The
Voigt fit is start-dependent: from the default phase start of 0 it stayed at
phase ~0 in 3/4 subjects with an SSR up to 22.5% above the minimum reached
from +-10 deg (not caused by the inert eta parameters: starting them at 0.5
changes nothing). Best of starts -10/0/+10 for both models:

| | Lorentzian | Voigt, wG from water |
|---|---|---|
| H6/H2 SD (4 subjects) | 0.144 | 0.151 |
| NADH2 residual spike, noise SD | 4.6 / 6.2 / 4.6 / 3.3 | 4.5 / 6.2 / 4.6 / 3.3 |
| NADH2 uM S02/S30/SXX/SXXX | 234 / 183 / 290 / 238 | 220 / 177 / 281 / 226 |
| SXXX run1 vs run2 NADH2/water | +17% | +13% |
| SSR spread across phase starts | < 1% | up to 22.5% |

The lineshape family is not what clips NADH2. In all four voxels the data
NADH2 peak is 20-35% taller than any fit and the data dip below the fitted
baseline just downfield (9.45-9.55 ppm); the fit cannot raise H2 without
overshooting that dip. The per-peak phases act as a shape fudge (Lorentzian
S30: H2 +20, H6 -20, both on the bound). The under-converged Voigt fits with
all phases held at 0 gave the tightest H6/H2 (SD 0.072), which points at the
per-peak phase freedom rather than the lineshape. Slab fits (ph +-60) do not
show the zero-start trap.

Restricting the per-peak phases (one shared phase, or a tight range around
the coarse correction) was proposed as the next test and declined by the
user (2026-09-24). The voxel stays on Lorentzian +-20, single region.
