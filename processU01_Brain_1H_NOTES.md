# processU01_Brain_1H — design notes, sweeps and rejected alternatives

Working notes for `processU01_Brain_1H.m` (brain 1H NAD+/Trp downfield fit).
These were previously inline comments in the .m file; they were moved here to
keep the processing script distributable. **The .m file is the source of truth
for what the code does; this file records why.**

Cohort referred to throughout as "the 12 controls" = the 12 `SXXX` subjects in
`DFS_CONTROL_STUDY`. "SpecTickle" is the time-domain HSVD program used for
cross-comparison; note there is **no gold standard** for NAD quantification, so
agreement with SpecTickle is corroboration at best, never proof.

Contents:
1. [Fit window and the two-region split](#1-fit-window-and-the-two-region-split)
2. [Baseline spline degrees of freedom](#2-baseline-spline-degrees-of-freedom)
3. [HSVD soft priors](#3-hsvd-soft-priors-width--baseline)
4. [Phase handling](#4-phase-handling)
5. [Lineshape kernel (disabled)](#5-lineshape-kernel-disabled)
6. [HSVD coarse seeding](#6-hsvd-coarse-seeding)
7. [Relaxation constants](#7-relaxation-constants)
8. [areasConv vs areas](#8-areasconv-vs-areas)
9. [File naming / pick_latest_pair](#9-file-naming--pick_latest_pair)
10. [Open issues and do-nots](#10-open-issues-and-do-nots)
11. [Aug 2026: the coarse-phase regression, bisected](#11-aug-2026-the-coarse-phase-regression-bisected)
12. [Aug 2026: coarse phase — the first-order term is noise](#12-aug-2026-coarse-phase--the-first-order-term-is-noise-coarsefirstorder)
13. [Aug 2026: the width prior across three cohorts](#13-aug-2026-the-width-prior-across-three-cohorts-lambdawidth)
14. [Aug 2026: Voigt lineshape — reparameterization, and what eta does](#14-aug-2026-voigt-lineshape--reparameterization-and-what-eta-does)
15. [Aug 2026: lambdaWidth revisited after the reparameterization](#15-aug-2026-lambdawidth-revisited-after-the-reparameterization---0005)
16. [Aug 2026: why U01 reads ~1.5x higher than DFS and Pilot-Sinc](#16-aug-2026-why-u01-reads-15x-higher-than-dfs-and-pilot-sinc)
17. [Aug 2026: the cross-over volunteer — it is the coil](#17-aug-2026-the-cross-over-volunteer----it-is-the-coil-not-the-protocol)
18. [Sep 2026: FIRST drops BrStem](#18-sep-2026-first-drops-brstem-new-stage-1c-in-processu01_brainps1)
19. [Sep 2026: `-Ref31P` — selectable 31P primary reference](#19-sep-2026--ref31p--selectable-31p-primary-reference)
20. [Sep 2026: external-reference B1+ from the 31P B1CALIB series](#20-sep-2026-external-reference-b1-from-the-31p-b1calib-series)
21. [Sep 2026: coarse-HSVD ridge and lambdaWidth re-selection](#21-sep-2026-coarse-hsvd-ridge-and-lambdawidth-re-selection)
22. [Sep 2026: HSVD results saved for comparison with VARPRO](#22-sep-2026-hsvd-results-saved-for-comparison-with-varpro)
23. [Sep 2026: HSVD vs VARPRO across every data set (`compare_hsvd_varpro`)](#23-sep-2026-hsvd-vs-varpro-across-every-data-set-compare_hsvd_varpro)

---

## 1. Fit window and the two-region split

### Why two regions

Two independent fit regions — the NAD+ cluster and the narrow Trp singlet at
~10.1 ppm — each running its own coarse → phase → fine fit (`coarseMode
'perRegion'`).

Trp's narrow window lets the local spline baseline flex enough to absorb the
broad macromolecular hump under Trp without inflating Trp, while the NAD window
keeps a stiff baseline that does not nibble the broad NADH4 peak. A single wide
region cannot do both at once:

- a stiff baseline inflates Trp ~2x
- a flexible one eats NADH4
- an explicit broad component (`BroadMM`, name-detected by `dat2mat_svsNAD`) is
  multimodal/fragile across sessions (verified on S007/S076)

**Retired single-region alternative** (Jul 2026; H4/H2 = 2.44 vs 1.083 expected,
Trp/H2 CV 45%):

```matlab
pars.fit.regions(1).fitRange = [8.7 10.3];
pars.fit.regions(1).peaks    = {'NADH2','NADH6','NADH4','Trp'};
pars.fit.regions(1).name     = 'NAD/Trp';
```

### Region 1 lower edge: the 8.7 → 8.6 → 8.7 round trip

**Widened 8.7 → 8.6 (Jul 2026, first pass).** Under `[8.7 10.0]` / `[9.95 10.3]`
the fraction of each peak's Lorentzian mass falling OUTSIDE its own fit window
was NADH2 3.2%, NADH6 5.0%, but **NADH4 15.3%** — NADH4 sits only ~1 linewidth
from the 8.7 edge, so its left wing was never measured and `areasConv`
extrapolated into it (`areasConv/areas` pinned at ~1.49 in 11/12 controls while
H2/H6 scattered 0.79–1.43). That is the bulk of H4/H2 = 2.44 against 1.083
expected. At 8.6 the H4 tail outside drops to ~8.5%.

**8.6 and not lower.** All 12 controls carry a narrow peak at 8.50–8.56 ppm
(lw 27–79 Hz, position consistent with ATP2 / adenine H8 — SpecTickle's table
puts ATP2 at 8.44) whose HEIGHT averages 3.3x NADH4's, max 5.0x. No peak is
assigned there, and the spline is far too stiff to absorb something that sharp,
so an 8.4 edge would drag the baseline up and push the error into H4 — worse
than the clipping it fixes. An 8.6 edge excludes that peak in all 12 (worst
intruder drops to 0.1x H4 height). **Going below 8.6 requires adding an explicit
nuisance peak near 8.52 first.**

**REVERTED to 8.7 (Jul 2026, second pass) — current setting.** Measuring
SpecTickle's own baseline showed the cost of 8.6: the number of cubic-spline
coefficients needed to follow the baseline to 2% of the NAD peak height is 8–12
(median 8.5) over `[8.6 10.0]` but only 5–11 (median 6) over `[8.7 10.0]`. The
extra 0.1 ppm is where the baseline turns sharply upward into the 8.5 ppm
ATP2/adenine peak, so including it buys a little less H4 clipping at the price
of a baseline the spline cannot follow — and the misfit lands on NADH4, the very
peak the widening was meant to fix.

*Known cost of reverting:* ~15% of H4's Lorentzian mass again falls outside the
window (vs ~8.5% at 8.6) and `areasConv` extrapolates into it. Watch H4/H2.
This is partly moot — see [scope note](#scope-nadh4-may-be-dropped) below.

### Region 2 width

Trp is intrinsically broad (widthL railed at its 0.2019 bound in 4/12), so a
0.35 ppm window left it degenerate with its own local baseline. 0.65 ppm gives
window:width ~3.2:1 instead of ~1.7:1.

### `pars.fit.ppm_range`

With `pars.fit.regions` defined this value only seeds `pars.fit.focusRange`, and
the focus weighting is inert unless `focusWeight > 1`
(`dat2mat_svsNAD.m:1576`) — so it is **currently a no-op**. Kept consistent with
the region edges so it does not silently down-weight them if `focusWeight` is
ever raised. Previous value `[8.7 10.0]`.

### SCOPE: NADH4 may be dropped

Jul 2026, per user: "not as worried about H4, may not use it." NADH2 is the
primary reporting peak (matches the workbook `[NAD-H2](mM)` column); the
self-consistency check then reduces to H6/H2 alone. **Score future tuning on
H6/H2, NADH2 CV, and NADH2 agreement — not on H4/H2.** Keep NADH4 as a MODELED
peak in region 1 regardless (removing it leaks its signal into H6 and the
baseline); just don't report it. This also settles the 8.6-vs-8.7 debate in
favor of 8.7.

---

## 2. Baseline spline degrees of freedom

### The knotSpacing sentinel

`knotSpacing = 2` is **not** a 2 ppm spacing. `makeSplineBasis` uses
`nIntervals = max(2, round(width/knotSpacing))`
(`curvefitAuto_varproBaseline.m:852`), and `round(width/2)` is 0 or 1 for any
window under 5 ppm, so the `max(2,...)` floor always wins. Every region
therefore gets exactly 2 intervals — a **5-coefficient** cubic spline —
whatever its width.

Coefficient count = `nIntervals + 3` for the cubic B-spline basis.

Consequence: widening a region does **not** add baseline parameters; it spreads
the same 5 over more ppm, i.e. it makes the baseline locally **stiffer**.

### Region 1: knotSpacing 0.43 → 6 coefficients

Measured against SpecTickle's unassigned HSVD components
(`analyze_spectickle_baseline.py`): over `[8.7 10.0]` the median subject needs 6
spline coefficients, i.e. 3 intervals over 1.3 ppm. `round(1.3/0.43) = 3` → 6
coefficients, up from the 5 the sentinel gives.

`0.325` (4 intervals, 7 coefficients) covers 9/12 subjects instead of 8/12 but
risks the spline nibbling NADH4, which is the failure mode this fit is most
sensitive to — so start at the median and only loosen if H4/H2 demands it.

### Region 2: knotSpacing 2 → 5 coefficients (current)

History, in order:

1. **0.175 → 7 coefficients.** Set to reproduce the knot INTERVAL of the old
   `[9.95 10.3]` window (`round(0.65/0.175) = 4` intervals ⇒ 7 coefficients).
   **That reasoning was wrong** — it matched the previous setting instead of
   asking what the baseline actually needs. Fitting SpecTickle's own unassigned
   components over `[9.8 10.45]` needs only 4 coefficients (4 in 10/12 subjects,
   5 in 2/12) — every subject needed FEWER than the 7 that override produced.
   The surplus flexibility is one half of the Trp inflation: an over-free spline
   can dive under a too-broad Trp Voigt, and the two split the pedestal
   arbitrarily.

2. **Sentinel 2 → 5 coefficients.** The `max(2,...)` floor, the closest the
   spline basis can get to the 4 actually required.

3. **0.22 → 6 coefficients — TRIED AND REJECTED.** It won on residual (0.0094 vs
   0.0127 in the Trp window) but the extra degree of freedom is spent
   unphysically. `check_baseline_monotonicity.py`:

   | config | monotone | mean % of window rising | dip depth (% of Trp height) | dip distance from Trp |
   |---|---|---|---|---|
   | 6 coef | **0/12** | 20.8% | 18.1% | 0.10 ppm |
   | 5 coef | 8/12 | 1.1% | **0.0%** | — (no interior minimum in any subject) |
   | SpecTickle | 6/12 | 5.8% | 6.5% | 0.12 ppm |

4. **Back to sentinel 2 → 5 coefficients (current).**

**Why the dip is inadmissible.** Inside `[9.8 10.45]` every baseline contributor
is centered BELOW the window — NAD H2/H6/H4 are not modeled in region 2, and
SpecTickle's component census puts ~105% of the Trp-window baseline in the tails
of components centered below 9.8. A sum of absorption-mode Lorentzian tails
approaching from below is **monotone decreasing** across the window. An interior
minimum requires a negative-amplitude component or residual dispersive content,
neither of which would track the fitted Trp position in every subject. The notch
is peak/baseline degeneracy: the spline digs under Trp so Trp can be larger.

**Lesson worth generalizing: a residual gain bought by an inadmissible baseline
shape is degeneracy absorbing mismatch, not model improvement.** A physical
constraint violated 12/12 outranks a small residual gain in a nonlinear fit.

At 5 coef there is no interior minimum in any of the 12 (the few non-monotone
cases are a slight upturn at a window EDGE, depth 0.0%). Trp = 0.221 mM
(CV 19.4%) vs 0.293 mM at 6 coef.

**Do NOT re-select region-2 DOF on Trp-window residual**, and **do NOT re-select
on agreement with SpecTickle**: per-subject Trp:NADH2 correlates NEGATIVELY
between the pipelines (r = −0.79 over all 12, −0.37 excluding S095), so
mean-level agreement is coincidental. The partition stays genuinely uncertain
(~±27%, Trp ~0.20–0.35 mM); 5 coef is the physically admissible end of that band.

**Better fix than a DOF cap**, if more flexibility is ever needed: constrain the
region-2 baseline *shape* (monotonicity / sign of curvature) rather than starving
it of coefficients.

Region 2's baseline is almost entirely the tail of components centered below
9.8 ppm, so inside the window it is smooth and nearly monotonic; there is no
local structure for extra knots to capture.

---

## 3. HSVD soft priors (width + baseline)

The coarse HSVD step already measures both of the things the fine fit otherwise
has to invent: the baseline under the peaks (its unassigned components — the
same partition SpecTickle draws) and each peak's linewidth. Both are fed in as
**penalties rather than constraints**, so the data can still overrule them.

**Why soft rather than a hard cap:** a bound forces the width to sit exactly on
the limit and dumps the remaining discrepancy into the amplitude, which is the
quantity we actually care about. The 45 Hz Trp cap tried first did exactly that
— it moved Trp from ~2.9x SpecTickle to ~3x BELOW it, trading one bias for
another. A penalty lets a genuinely broad peak stay broad if the data insist,
while still stopping a peak from running away to swallow baseline.

Parameters:

- `lambdaPrior` — penalty is `||w.*(B*beta - hsvdBaseline)||^2`, weighted exactly
  like the data term, so 1.0 = "trust the prior as much as the data".
- `lambdaWidth` — normalized so 1.0 = "a 100% relative width error costs as much
  as a 100% relative data misfit". Voigt total FWHM (Olivero–Longbothum) is
  compared against the HSVD linewidth; peaks HSVD did not find are skipped.

**Sweep (Jul 2026, 4×4 over {0, 0.05, 0.1, 0.3}, all 12 controls, scored on
mean |NAD ratio − 1|):**

```
  lamP  lamW   H6/H2         H4/H2         mean|dev|  CV H2   CV H4
  0.00  0.00   1.456+-0.392  0.517+-0.115    0.470    15.7%   35.4%
  0.00  0.30   1.210+-0.105  0.895+-0.198    0.201    11.4%   24.0%   <-- used
  0.30  0.00   2.257+-1.008  2.457+-1.801    1.458    14.3%   63.9%
  0.30  0.30   1.166+-0.165  1.276+-0.614    0.311    16.2%   58.4%
```

**Finding 1 — the WIDTH prior does essentially all the work, and it saturates.**
Going from `lambdaWidth` 0 → 0.05 halves mean|dev| (0.470 → 0.209); everything
from 0.05 to 2.0 is flat within noise (0.209 → 0.198). So the value barely
matters, only that the prior exists. 0.3 sits mid-plateau.

**Finding 2 — the BASELINE prior is actively HARMFUL and is therefore OFF.** At
`lambdaWidth = 0` it makes things much worse than no prior at all (0.470 →
0.935/1.126/1.458 as `lambdaPrior` rises), and at every fixed `lambdaWidth` it
degrades the result monotonically (0.201 → 0.251 → 0.264 → 0.311). The reason is
the one `nadCoarsePhase`'s own header warns about: HSVD's `lw_max` criterion
truncates broad wings, so signal that genuinely belongs to the NAD peaks is left
in the unassigned "baseline". Pulling the spline toward that baseline hands peak
area to the baseline — NADH4's between-subject CV blows up to 63–90%. **The HSVD
baseline is a good description of what SpecTickle CALLS baseline, not of what is
actually under the peaks.**

The infrastructure stays (`lambdaPrior = 0` disables it cleanly) because the same
hook would work with a better baseline estimate — e.g. one from an HSVD run whose
`lw_max` does not truncate wings, or a measured MM spectrum.

---

## 4. Phase handling

### `ph_range = [-60 60]` — wide, and why that is still unresolved

**Original note (May 2026).** Wide phase range is intentional. Tightening
`ph_range` to `[-15 15]` (and/or making `lineShapeOpt` symmetric / lowering its
`maxSide`) makes the per-component figures LOOK cleaner — pure-absorptive peaks,
baseline no longer floats above the data — but tested, it makes
NADH2:NADH6:NADH4 LESS self-consistent across averages/sessions. The
interpretation at the time: the wide phase + asymmetric lineShape kernel are
absorbing real systematic distortion (eddy currents, B0 asymmetry, residual
phase drift); take that freedom away and the distortion is forced into the peak
amplitudes, where it pulls each NAD proton unequally.

**AMENDED Jul 2026 — the above diagnosis looks wrong.** The lineShape kernel it
refers to is now DISABLED, so the "wide phase + kernel together" reasoning no
longer describes the config. More importantly: measured on the 12 controls, the
coarse-fit phase correction left a residual ramp of ~74 deg/ppm (~118 deg across
this window), and the applied `pc1` correlated only r = 0.34 with what needed
correcting. The per-peak phases the fine fit then landed on rose monotonically
with ppm (NADH4 +2, NADH6 +34, NADH2 +37, Trp +47 deg, Trp railing at +60 in
4/12) — the signature of an **uncorrected linear ramp**, not of eddy currents.

So the wide `ph_range` was most likely absorbing a bad first-order phase estimate
(unweighted polyfit through 4 peaks at `dat2mat_svsNAD.m:1366`, dominated by
low-SNR Trp at the longest lever arm) rather than real physics. **Re-test
narrowing `ph_range` only AFTER the coarse phase step is fixed; the May 2026
test was confounded by this.**

### `coarseMode` — was `'perRegion'`, now `'union'` (Aug 2026)

**SUPERSEDED, see §11.** The driver now sets `coarseMode = 'union'`: one HSVD
and one phase shared by both regions, with the baseline partition still done per
region. The reason is the flip side of the argument below — a single-peak Trp
region cannot fit a slope at all, so `pc1 = 0` was not a considered choice but a
fallback, and on U01 the Trp region got **no** phase correction whatsoever
(`pc0 = pc1 = 0`) whenever HSVD failed to seed Trp (1/11 sessions). Region 1
meanwhile carried `pc1` of −2363 to +295 deg. The same spectrum was being phased
two different ways.

The original argument for `perRegion`, kept because it still describes a real
cost of `union`:

Chosen for the baseline reasons in §1, but it also improves the phase step for an
independent reason: with `perRegion` the Trp region contains a single peak, and
`nadCoarsePhase` then takes `pc0 = ph0, pc1 = 0` — a pure zeroth-order phase with
no ramp. That

- removes Trp from the NAD ramp fit, where as the lowest-SNR peak at the longest
  lever arm it dominated the slope, and
- stops Trp railing at the ±60 deg `ph_range` bound, which it did in 4/12
  controls under one region.

*Cost:* the NAD ramp is now fit from 3 peaks spanning only ~0.5 ppm, so the slope
is less well conditioned — watch `coarse.phResidDeg` for region 1.

### Phase reconstruction gotcha

The pivot in `shiftSpectrumPhase` is `center_freq` = **4.72** (water reference),
NOT `exc_freq` = 9.7. `pc0` is therefore the phase extrapolated to 4.7 ppm, which
is why `pc0`/`pc1` look implausibly large (hundreds of deg). Applied phase =
`pc0 + 2*pc1/span_full*(ppm - 4.72)`; with the right pivot this reproduces the
stored spectrum to 1e-16. Real slopes are −18 to −123 deg/ppm.

Also: `output.met.fit.spec` already holds the phase-corrected data
region-by-region (`specFull(inds_r) = y_r`), so no reconstruction is needed for
plotting — but **region 2 overwrites region 1 on their `[9.8 10.0]` overlap**, and
the two regions carry independent phase corrections, so there is a genuine step
at 9.8 ppm in the stitched arrays.

---

## 5. Lineshape kernel (disabled)

Shared lineshape kernel DISABLED (Jul 2026). The kernel and the Voigt Gaussian
width are collinear — both describe non-Lorentzian line distortion — and the
kernel was winning: across the 12 controls `widthG` sat railed at its LOWER bound
(0.0336 ppm) in 11/12 fits for NADH2 and 8/12 for NADH6. The fit was therefore
effectively Lorentzian-convolved-with-a-free-kernel, with the Voigt's Gaussian
term pinned at zero and contributing nothing.

Large voxels at 7T should produce genuine Gaussian character from B0
inhomogeneity, so the Voigt is the physically motivated model; turning the kernel
off lets `widthG` float and actually carry it.

**Pre-registered test, now settled:** with the kernel off, `widthG` *still* rails
at its lower bound in 18/36 NAD peak fits (all three peaks railed in 4/12
subjects). Per the criterion above, the Gaussian term is not earning its
parameter — **mode 5 (Voigt) should drop to a pure Lorentzian.** This also
matches SpecTickle's success with 1–2 Lorentzians per peak. *Not yet actioned.*

Settings kept (disabled) rather than deleted so they survive for A/B tests:
`nSide = 4`, `asymmetric = true`, `maxSide = 0.50`.

---

## 6. HSVD coarse seeding

Replaces the coarse VARPRO. **The VARPRO path is no longer commented out** — as
of Aug 2026 it is a live local function `coarseVarpro` in `dat2mat_svsNAD`,
selected by `hsvdCoarse.enable = false`. Between 041082f and then, that flag was
parsed but never branched on and the code it named was commented out, so setting
it `false` silently did nothing. See §11 for the A/B: HSVD wins on U01, VARPRO
wins by 3x on Pilot-Sinc.

`K` and `npts` match SpecTickle's defaults (`hsvd_order` 60, `hsvd_np` 1024).
The pole solve does **not** match: `utils/hsvd_fit.m` uses a ridge where
SpecTickle uses plain least squares, which widens every HSVD linewidth by
ln(1+λ)·bw/π Hz. Removing it was tested and rejected — see §21.

`lw` / phase acceptance windows follow SpecTickle's peak table for 7T 1H brain:
NAD+ 10–90 Hz, TRYP 10–120 Hz. **Note** the archived SpecTickle runs actually
used 80 Hz for TRYP, which was demonstrably truncating Trp (a co-located
86–120 Hz component was being thrown to baseline in S092/S095) — 120 Hz is the
current table value and the safer one for seeding.

`phaseTolDeg = 100` is the current SpecTickle table value and **matters**: the
archived runs used ±180 (criterion effectively off), and because our
`pars.peaks.range` windows are wider than SpecTickle's `ppm_tol` they then admit
inverted junk. Verified on S095, whose H6 window `[9.05 9.25]` pulls in a large
−164 deg component at 9.064 ppm: at ±180 the H6 seed phase comes out −178 deg
(peak effectively inverted) and wrecks the ramp fit; at ±100 that component is
rejected and the seed lands at +37 deg. Worst-case residual phase across the 12
controls: 90 deg → 34 deg.

Vectors are in `pars.peaks.name` order and get subset per fit region.

---

## 7. Relaxation constants

Literature relaxation / water content: Swago 2025, NMR Biomed 5324.

**Trp: indole NH proton (1H) at 10.1 ppm — T1/T2 NOT FINALIZED.** Jul 2026,
changed from 100/10 ms to **75/19 ms** to match the values SpecTickle uses
(`DFS_QUANT.xlsx` constants row: `T1 TRP (ms)` = 75, `T2 TRP(ms)` = 19), so that
a Trp concentration from this pipeline is directly comparable with SpecTickle's
rather than differing by a constant relaxation factor.

At TE = 13 ms / TR = 600 ms this raises Trp by ~1.6x (`exp(-13/10)/exp(-13/19)` =
1.35 from T2, and the T1 term moves 0.9975 → 0.9997).

**These are placeholders pending a measured Trp T1/T2 — revisit before reporting
absolute Trp.** NAD+ values are unchanged and NADH2's 205.6/33.6 already match
SpecTickle's constants row exactly.

---

## 8. areasConv vs areas

Use **kernel-aware areas**. `areas` is base-Voigt (`ampl .* width`) only and
under-counts whenever `pars.fit.lineShapeOpt` is enabled; `areasConv` integrates
the convolved component (see `dat2mat_svsNAD` / `curvefitAuto_varproBaseline`).
When the water fit has no kernel and the metabolite fit does, swapping
`areas` → `areasConv` raises the reported metabolite concentrations by the kernel
factor (~10–40% with the settings that were current when this was measured).

**`areasConv` is phase-invariant (Jul 2026, settled — do not reconsider).**
Quantitation must not depend on the fitted phase: the reported quantity is the
magnitude of the complex amplitude, evaluated as the frequency-domain area of the
**zero-phase** model peak (so it still captures kernel / Voigt shape that a bare
amplitude would miss). `curvefitAuto_varproBaseline.m` zeroes the per-peak phase
block before the extended-axis integration.

Effect when introduced: Trp 0.196 → 0.221 mM, NAD H6/H2 SD 0.105 → 0.081, overall
mean|dev| 0.201 → 0.185.

---

## 9. File naming / pick_latest_pair

Two naming conventions are recognized, because the U01 sessions and the
DFS_CONTROL_STUDY export label the same two scans differently:

| | U01 | DFS export |
|---|---|---|
| water-suppressed | `_NAD_` | `OSrem` |
| water reference | `_Water_` | `waterref` |

The DFS names carry no `_NAD_`/`_Water_` token at all — the metabolite file is
e.g. `meas_MID00129_..._svssel_sinc2x2_2x9_7ppm_OSrem.dat` and the reference is
`meas_MID00128_..._svssel_sinc2x2_waterref.dat` — so the old matcher found
neither and asserted.

Note `waterref` deliberately does not match the `_Water_` test (no trailing
underscore), so the two conventions cannot cross-match; the water-suppressed test
is applied first regardless.

**Known failure: Pilot-Sinc-DFMRS.** The water-suppressed test now also accepts
`ppm`, and the Pilot-Sinc water references are named
`..._svssel_WATERREF_2x4_7ppm.dat` — they contain `ppm`, so the WS test fires
first (it is an `elseif` chain) and *every* file lands in the WS bucket, leaving
`nWAT = 0` and an assert. Even with that fixed, highest-MID would pair the 9.7 ppm
metabolite scan with the water reference belonging to the 10.1 ppm scan. That
study needs explicit pairing — `*2X9_7ppm*` and `*WATERREF_2x4_7ppm*`, which is
what the original `process_repeatability.m` used and what reproduces the archived
results. Do not point `pick_latest_pair` at Pilot-Sinc without fixing this.

**Previous version (U01 naming only), kept for reference:**

```matlab
for k = 1:numel(allF)
    fn = allF(k).name;
    tok = regexp(fn, 'MID(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end          % require MID# to disambiguate
    midNum = str2double(tok{1});

    if contains(fn, '_NAD_', 'IgnoreCase', true)
        nNAD = nNAD + 1;
        nadFiles(nNAD) = string(fn);
        nadMid(nNAD)   = midNum;
    elseif contains(fn, '_Water_', 'IgnoreCase', true)
        nWAT = nWAT + 1;
        watFiles(nWAT) = string(fn);
        watMid(nWAT)   = midNum;
    end
end

assert(nNAD > 0, 'No water-suppressed (*_NAD_*) .dat in %s',   datDir);
assert(nWAT > 0, 'No water-reference (*_Water_*) .dat in %s',  datDir);
```

**Session layout.** Sessions are normally `<SessionDir>\datfiles\*.dat`, but some
exported/archived studies (e.g. DFS_CONTROL_STUDY) keep the .dat files directly
in the session folder. Previous behavior (required `datfiles\`):

```matlab
datDir = fullfile(SessionDir, 'datfiles');
assert(isfolder(datDir), 'datfiles\\ subfolder not found in %s', SessionDir);
```


**Sep 2026: `MID` is not always there — the gate now accepts `MID` or `FID`.**
Every U01_RR session exports as `meas_MID#####_FID#####_...`, but
`859875_U01\Brain\859875_S01_09012026_brain` exports with the `meas_MID#####_`
prefix stripped:

```
859875_S01_09012026_brain_FID29330_svssel_20260101_Water_8avg_tra.dat
859875_S01_09012026_brain_FID29331_svssel_20260101_NAD_9_7ppm_1024avg_tra.dat
859875_S01_09012026_brain_FID29348_isis_svs_20260901_31P_1D.dat
```

The disambiguation gate was `regexp(fn, 'MID(\d+)', ...)` with
`if isempty(tok), continue; end`, so **both files were skipped before the
`_NAD_` test ever ran** and the assert reported "No water-suppressed (_NAD_ /
ppm) .dat" for a directory that plainly contains one. Misleading message: the
failure was the MID gate, not the token match.

Now `'[MF]ID(\d+)'` in all four drivers (`processU01_Brain_1H:237`,
`processU01_Brain_31P:397`, `processU01_Calf_1H:211`, `processU01_Calf_31P:363`).
Leftmost-match keeps the old behavior on standard names — `meas_MID01248_FID12878`
still resolves to 1248, not 12878 — so existing sessions are unaffected.

Caveat: MID and FID are separate counters, so a directory mixing both naming
styles would compare across scales and the FID-named file would always win the
`max()`. Real directories are uniform, so this is not currently reachable; if a
mixed directory ever appears, prefer MID when any file has one.

---

## 10. Open issues and do-nots

**Do not:**

- Re-select region-2 spline DOF on Trp-window residual, or on agreement with
  SpecTickle (§2).
- Narrow `ph_range` before the coarse phase step is fixed (§4). Part of that step
  *has* now been fixed (§11), so this is worth revisiting — but only on the HSVD
  path, which is what U01 uses.
- Assume a setting that helps one study helps another. `knotSpacing`, the
  lineshape kernel, `lambdaPrior` and the coarse estimator all went opposite ways
  on U01 vs Pilot-Sinc (§11). Tune per acquisition, and say which data a number
  came from.
- Tune anything on the U01 test–retest alone. It is **2 subject pairs**; the
  λ<sub>prior</sub> = 0.5 recommendation derived from it in Aug 2026 was reversed
  by 4 more pairs (§11).
- Reconsider the phase-invariant `areasConv` (§8) — settled.
- Report H4/H2 as a tuning target (§1, scope note).
- Treat agreement with SpecTickle as accuracy. There is no gold standard; the
  H2/H6/H4 ratio test is *self-consistency*, not truth.
- Remove the ridge from the coarse `hsvd_fit` call, on its own or together with
  a `lambdaWidth` re-sweep, or subtract its width offset from the width prior
  target. All three were run on U01 + Pilot-Sinc + DFS and lose to the current
  config (§21).

**Open:**

- `widthG` still rails with the kernel off ⇒ mode 5 should become pure
  Lorentzian (§5). Not yet actioned.
- Trp T1/T2 are placeholders (§7).
- Per-subject concentrations correlate poorly between this pipeline and
  SpecTickle (NADH2 r = 0.016; Trp:NADH2 r = −0.79). At most one is tracking real
  biology. No age effect survives in either arm (n = 12, only 2 subjects > 40 y).
- S167 is a clear failure: H6/H2 = 5.33 and region-1 `pc1` = +967 deg against a
  cohort median of −679 deg (unwrapping-branch flip), despite small phase
  residuals — **small `phResidDeg` does NOT prove the ramp is right.**
- Going below 8.6 ppm on the region-1 edge needs an explicit nuisance peak near
  8.52 first (§1).

**Current results at the config in the .m file** (5 coef region 2, all 12
controls): NADH2 0.402 mM (CV 11.4%), NADH6 0.467 (13.8%), NADH4 0.384 (28.3%),
Trp 0.221 (19.4%), H6/H2 1.160 ± 0.081. Every between-subject CV is lower than
SpecTickle's (23.0 / 21.0 / 28.6 / 29.7%).

**Helper scripts** (in `DFS_CONTROL_STUDY`, not in this repo):
`reprocess_second_fit.m` (batch runner), `compare_second_vs_spectickle.py`,
`analyze_spectickle_baseline.py`, `check_baseline_monotonicity.py`,
`plot_fit_comparison.py`.

---

## 11. Aug 2026: the coarse-phase regression, bisected

Prompted by reprocessing the Pilot-Sinc-DFMRS 9.7 ppm repeatability data (4
subjects x 2 visits) with the current pipeline and finding NAD+ test–retest CV
of 5.3% against 1.1% in the archived `PROCESSED/AbsoluteQuant_9_7.tsv`.

### The archived result is real and exactly reproducible

`git worktree` at **994b652** (the commit that first contains
`demo_absoluteQuant_svs.m`; the TSV predates it by 6 days) run today reproduces
that TSV **byte for byte** — same md5, under R2026a. So the 1.1% is not an
artifact, not lost to an uncommitted state, and not MATLAB-version dependent.

### MATLAB R2024a vs R2026a: only the old water fit was sensitive

Comparing the archived `ResultsMatlab2024a` / `ResultsMatlab2026a` folders
(same code, two MATLAB releases):

| | mean \|change\| | max |
|---|---|---|
| metabolite areas | 0.04–0.08% | 0.35% |
| **water reference areas** | **5.2–7.0%** | **14.4%** |

Because concentration is metabolite/water, that propagated into CVs erratically —
S139 NADH6 went 4.3% → 20.3% purely from the MATLAB version. That was
`process_repeatability.m`'s `watfit.method='fid'`, `nComp=7`. **The current
pipeline is immune**: `watfit.method='fit'`, `nComp=1` is well-conditioned enough
that R2024a and R2026a agree to five significant figures on every scan and every
test–retest CV. Any archived result predating that change carries an unknown
version-dependent offset.

### Bisect (harness fixed at the 994b652 script; only the fitting code varies)

Pilot-Sinc 9.7 ppm test–retest CV%, mean over the 4 pairs:

| commit | what changed | H2 | H6 | H4 | NAD+ mean |
|---|---|---|---|---|---|
| 994b652 | reproduces published | 1.89 | 3.23 | 1.84 | **1.13** |
| 7ac5156 | 31P + bruker helpers | 1.90 | 3.32 | 1.90 | 1.83 |
| 273f7f4 | 31P/1H fitting | 3.44 | 1.53 | 2.25 | **1.17** |
| **1a4640b** | **external referencing** | 4.30 | 4.82 | **12.16** | **3.30** |
| 041082f | coarse fit → time domain | 6.01 | 5.50 | 12.43 | **5.09** |
| de76847 | bug fixes and cleanup | 5.53 | 4.78 | 11.73 | 4.62 |

Two defects, roughly equal weight.

### Defect 1 — `conj(Vc)` vs the phase correction sign

`1a4640b` is a 1948-line commit about 31P external referencing that never touches
`dat2mat_svsNAD.m`. It added 54 lines to `utils/curvefitAuto_varproBaseline.m`,
of which two (`baselineOpt.phaseInit`, `lineShapeOpt.init`) are gated on fields
this pipeline never sets. The third is ungated and sits in the Voigt evaluator:

```matlab
- yfit = real((Vc / N0) * exp(-1i*ph));
+ yfit = real((conj(Vc) / N0) * exp(-1i*ph));
```

The fix is **correct**: `fadf` returns w(z) with Im > 0 for u > 0, while the
sigma=0 Lorentzian branch and `compositel_complex_unit` use the NMR convention
with Im < 0 for x > p. Without it the phase parameter flipped sign
discontinuously at sigma=0 and was inverted relative to mode 3.

But nothing downstream was updated. At that commit the coarse phase came from
`curvefitAuto_varproBaseline` running mode 5 *itself*, so the measured `ph0`
inverted while `shiftSpectrumPhase(ws, [pc0 pc1], ...)` kept adding it — the
correction roughly **doubled** the phase error instead of canceling it. Fits
still converged with clean residuals, but repeat scans landed in different basins.

Proven three ways at 1a4640b: reverting `conj` alone → 1.17 (identical to
273f7f4); keeping `conj` and negating the applied phase → **1.15**; doing
neither → 3.30.

**Resolution — the sign is per estimator, not global:**

| coarse estimator | phase source | correct sign | evidence |
|---|---|---|---|
| HSVD (041082f+) | `angle(A)` of components | `+pc0` | negating at HEAD: 4.62 → 5.75 (worse) |
| VARPRO (`coarseVarpro`) | the conj'd Voigt's own `ph` | **`-pc0`** | not negating: 1.49 → 6.61 (worse) |

`041082f` incidentally repaired defect 1 by replacing the phase source. `conj(Vc)`
stays — it is correct, and costs nothing once the sign matches it.

### Defect 2 — the coarse estimator swap (still open on U01)

`041082f` replaced the VARPRO coarse fit with time-domain HSVD, +1.79 points.
That is not a bug, it is a substitution that regressed on this acquisition.
A/B in current code, VARPRO coarse with the corrected sign:

| dataset | HSVD coarse | VARPRO coarse (`-pc0`) |
|---|---|---|
| **Pilot-Sinc** NAD+ mean TRT | 4.62% | **1.49%** |
| **U01** NAD+ mean TRT (n=2 pairs) | 8.7% | 8.3% |
| **U01** NAD+ mean group CV (n=6) | **11.9%** | 18.3% |
| **U01** NADH6 group CV (n=6) | **10.7%** | **32.8%** |

**They disagree, so the driver keeps HSVD.** Pilot-Sinc has a far steeper
baseline ramp (sinc excitation centered at 9.7 ppm; the spline baseline carries
~70–100% of the peak height), which the VARPRO coarse tracks and HSVD does not.
Use `hsvdCoarse.enable = false` + `coarseMode 'perRegion'` if reprocessing
Pilot-Sinc; leave U01 alone.

Residual: even fixed, Pilot-Sinc NADH4 only recovers 11.73 → 6.00 (2.25 before),
so something in 041082f/de76847 still costs H4 specifically. Not chased.

### Trp (added after the fact — it was NOT in the bisect)

The archived Pilot-Sinc analysis fit **three peaks and one region**
(`peaks.name = {'NADH2','NADH6','NADH4'}`, `ppm_range = [8.7 10.0]`, no
`pars.fit.regions`). Every bisect row above is therefore NAD-only. The 4-peak /
two-region `demo_absoluteQuant_svs.m` in the repo today is the 273f7f4 version,
written *after* those results — do not read its `regions(2) = [9.95 10.30]` as
describing the archived config.

Trp *was* fit in the reprocessing runs, which used `processU01_Brain_1H`
(4 peaks, NAD `[8.7 10.0]` / Trp `[9.8 10.6]`). Trp on the 8 Pilot-Sinc scans:

| config | Trp TRT CV% (4 pairs) | group mean ± SD (µM) | group CV% | Trp width/ub |
|---|---|---|---|---|
| union, HSVD, λ=0 | 14.3 | 98.5 ± 37.2 | 37.8 | 28% |
| union, HSVD, λ=0.3 | 20.2 | 123.8 ± 65.8 | 53.2 | 38% |
| union, HSVD, λ=0.5 | 20.0 | 128.7 ± 70.9 | 55.1 | 44% |
| perRegion, HSVD, λ=0 | 16.9 | 93.8 ± 33.5 | 35.8 | 28% |
| **perRegion, VARPRO, λ=0** | **4.1** | 71.9 ± 19.8 | **27.5** | **12%** |

Same direction as NAD+ on this study: the VARPRO coarse wins, and by a lot
(4.1% vs 14.3%, per subject 1.6 / 3.9 / 8.4 / 2.3%). λ<sub>prior</sub> hurts Trp
here exactly as it does on U01.

**But treat the VARPRO Trp as precision without demonstrated accuracy.** The
group mean falls 27% and the fitted Trp width drops to 12% of its bound (~14 Hz);
on inspection the model puts a narrow line on the 10.15 ppm step and lets the
spline carry most of it. It reproduces beautifully because it is consistently
fitting the same small feature — which is not the same as fitting Trp correctly,
and there is no gold standard to appeal to.

union beats perRegion for Trp under HSVD (14.3 vs 16.9), which supports the
phase-sharing argument in §4, but the effect is small next to the coarse
estimator.

### `R1 res%` is not comparable across coarseMode when regions overlap

Region 1 `[8.7 10.0]` and region 2 `[9.8 10.6]` overlap by 0.2 ppm. Under
`perRegion` the two carry independent phases — measured on S004_01092025,
region 1 `pc0/pc1 = -264.7/386.8` vs region 2 `-17.7/0.0` — and
`specFull(inds_r) = y_r` means region 2 **overwrites** region 1 on the overlap.
The stored `output.met.fit.spec` therefore has a genuine step at 9.80 ppm (§4
already noted this), and any residual computed as `spec(region-1 window) − model`
is inflated over ~15% of its points. Under `union` the phases are identical and
no step exists.

So a `perRegion` vs `union` residual comparison measures the stitch, not the fit.
Concentrations are unaffected — each region is fit against its own correctly
phased data. `perRegion`-vs-`perRegion` comparisons (e.g. the U01 HSVD/VARPRO
A/B) remain valid because the artifact is common to both.

### Code changes made

In `dat2mat_svsNAD.m`:

- `coarseVarpro` restored as a live local function behind `hsvdCoarse.enable`,
  applying `[-pc0 -pc1]` for the reason above.
- Hard error on `coarseMode 'union'` + `hsvdCoarse.enable = false` — union
  re-partitions a component set the VARPRO path never produces.
- `union` mode now re-partitions the single HSVD component set **per region**
  (`unionCoarseForRegion`), so the Trp region's baseline still carries the NAD
  tails. Previously union stored one coarse struct for every region, whose grid
  never matched a narrower region, so `lambdaPrior` was silently dropped.
- Assignment/reconstruction factored into `hsvdAssign` / `hsvdSpectra` /
  `hsvdPartition` so the two callers cannot drift.
- A `warning` instead of a silent no-op when the baseline prior's grid mismatches.

Regression check: current code with HSVD coarse reproduces de76847 to the digit
(5.53/4.78/11.73/4.62), so none of the above changed HSVD-path behavior.

### Reproducing any of this

Scratch harnesses (session scratchpad, not the repo): `pilot_sinc_cfg.m` runs
either study with arbitrary overrides of `lambdaPrior` / `knotSpacing` / kernel /
`coarseMode` / `lambdaWidth` / `hsvdCoarse.enable`, writing renamed shadow copies
of the driver so the repo file is never edited; `report_bisect.py` tabulates
`AbsoluteQuant_9_7.tsv` across code states. The bisect used `git worktree add`
per commit with the 994b652 script pinned as the harness and its `addpath` and
`outDir` redirected — **the published TSV must not be overwritten**; verify its
md5 (`b0127d9fed50995d624ee34420473379`) before and after.

## 12. Aug 2026: coarse phase — the first-order term is noise (`coarseFirstOrder`)

**Outcome: `pars.fit.coarseFirstOrder` added, DEFAULT `false`.** The coarse
phase is now zeroth order only, for **both** coarse estimators. Setting it
`true` restores the previous behavior exactly (verified bit-identical).

### `pc0`/`pc1` were never comparable across studies

`shiftSpectrumPhase` applies

    phi(x) = pc0 + 2*pc1/|ppm(end)-ppm(1)| * (x - pivot)

and the pipeline sets `pc1 = slope*|dppm|/2`, `pc0 = line(pivot)`, so this
collapses to `phi(x) = slope*x + intercept` — **independent of the pivot**.
`pc0` is therefore just the fitted line read off at `center_freq`, and

| study | `center_freq` | peak span | extrapolation |
|---|---|---|---|
| U01_2026 | **9.70 ppm** | 8.88–10.15 | +0.5 ppm (inside) |
| Pilot-Sinc | **4.72 ppm** | 8.88–10.15 | **−4.4 ppm (outside)** |

(`center_freq` comes from `seqname` / `length(wipdbl)` near the top of
`dat2mat_svsNAD`.) On Pilot-Sinc `pc0` is an extrapolation 4.4 ppm beyond the
data, so slope error is multiplied by 4.4. That fully explains an apparent
"111° within-pair phase instability" that does not exist: measured as
`slope × extrap`, e.g. −60.0 °/ppm × −4.38 ppm = +263 against an observed
`pc0` of 262.4. Within pairs, `Δslope × 4.4` predicts `Δpc0` to a few degrees
on all four Pilot-Sinc pairs.

**Do not judge coarse-phase quality by `pc0` or `pc1`.** Compare `phi(x)` over
the fit window.

### Why the first-order term was dropped

The slope is fit from 3–4 peaks spanning ~1.25 ppm whose individual phases
scatter by ~20°, giving a slope uncertainty near 16 °/ppm. Measured slopes ran
**−69 to +20 °/ppm** across scans, and `pc1` reached −934° (U01 S138_04102026).
It is fitting noise, and it is then applied across the whole window.

Median within-pair |Δphi| over 8.9–10.15 ppm, 6 test–retest pairs:

| study | HSVD | HSVD `pc1=0` | VARPRO | VARPRO `pc1=0` |
|---|---|---|---|---|
| U01 | 16.3° | 10.0° | 8.0° | **4.2°** |
| Pilot-Sinc | 10.7° | **3.7°** | 2.8° | **1.1°** |

2.5–2.9× better for both estimators, in both studies. It never hurt.

### Effect on quantitation — large for VARPRO, neutral for HSVD

Test–retest CV% (U01 = 2 pairs, Pilot-Sinc = 4):

| config | U01 H2/H6/H4 | PS H2/H6/H4 |
|---|---|---|
| archive (VARPRO) | 14.6 / 9.7 / 11.6 | 4.0 / 2.7 / 6.6 |
| archive + `pc1=0` | **4.0 / 6.2 / 6.1** | 5.6 / **1.1** / 5.0 |
| ref (HSVD) | 6.2 / 13.0 / 11.4 | 5.7 / 6.3 / 7.9 |
| ref + `pc1=0` | 7.8 / 11.2 / 12.0 | 4.7 / 6.5 / 8.8 |

`pc1=0` rescues the VARPRO/archive config on U01 (H2 14.6% → 4.0%) while also
improving Pilot-Sinc H6 to 1.1%. For the HSVD path it is **outcome-neutral**
despite phase reproducibility improving 2.9× — so phase is not what limits the
current pipeline's test–retest; do not chase it further there.

Caveat: archive+`pc1=0` still has poor U01 **group** CV (28.4% H6 vs ~12–13%
for the ED configs), and its U01 test–retest rests on 2 pairs — see §10 and
`feedback-dont-tune-on-two-pairs`.

### Ruled out along the way

- **HSVD component multiplicity is NOT the cause of anything on Pilot-Sinc.**
  `hsvdAssign` complex-sums every qualifying component (`A = sum(amp(sel))`,
  `phDeg = angle(A)`), so multiplicity *could* destabilize phase. Measured: on
  every Pilot-Sinc scan each peak collects exactly **one** component, coherence
  `|sum(amp)|/sum(|amp|)` = 1.00. U01 shows some multiplicity (NADH4 in 50% of
  scans) with coherence still 0.97.
- **The sinc excitation profile is not an explanation for anything** here. It
  has been simulated and tested experimentally and is fine for Trp, H2 and H6;
  possibly slight for H4. Do not invoke it.

### Reproducing

Scratch harnesses (session scratchpad, not the repo): `ed_cfg.m` runs either
study from an **isolated copy of the tree** (`scratchpad/edtree`) — path
shadowing does **not** work for anything in `utils/`, because
`dat2mat_svsNAD.m:22` does `addpath([thisPath filesep 'utils'])` on every call
and re-prepends the repo copy ahead of any shadow dir. Two smoke tests silently
ran the unmodified fitter before an execution marker caught it. `ed_score.m`
scores test–retest / group CV / H2:H6; `applied_phase.m` compares `phi(x)` over
the peaks; `hsvd_components.m` counts components per peak; `hsvd_phasefit.m`
dumps slope / extrapolation / `pc0` / `pc1` per scan.

## 13. Aug 2026: the width prior across three cohorts (`lambdaWidth`)

> **SUPERSEDED ON THE VALUE -- see section 15.** The driver now runs
> `lambdaWidth = 0.005`. Everything below about *why* the prior exists, why Trp
> is exempt, and how the three cohorts disagree still stands; only the chosen
> lambda changed, after the Voigt reparameterization (section 14) removed the
> 10 Hz component floor that had been acting as an accidental regularizer.

**Outcome (at the time): `lambdaWidth = 0.3` with `widthPriorPeaks = [true true true false]`.**
NAD+ widths pinned to the HSVD coarse estimate; Trp fitted freely.

### The three cohorts disagree, and the disagreement is real

| | U01_2026 | Pilot-Sinc | DFS_CONTROL_STUDY |
|---|---|---|---|
| n | 10 (1 excluded) | 8 | 12 |
| repeat pairs | 2 | 4 | **0** |
| sequence | `svssel_20260101` / `_latest` | `svssel_latest` sinc 9.7 | `svssel_sinc2x2` 9.7 |
| what it can measure | TRT + group CV | TRT + group CV | group CV + self-consistency only |
| prefers width prior | mixed | **OFF** | **ON, strongly** |

Pilot-Sinc with the prior off: H2/H6 group CV 11.3/11.2 -> **6.4/6.4**, test-retest
H6 6.5 -> 5.1%. DFS with the prior off: NAD+ self-consistency MAD from 1.0
0.226 -> **0.418**, H4/H2 1.025 -> **0.572**, NADH4 concentration **-54%**.

DFS is the largest single cohort and the only one with no repeats, so group CV
and the H6/H2 + H4/H2 ratios are all it can speak to -- but on those it is
emphatic. Pilot-Sinc is the only cohort with 4 repeat pairs. **They are
measuring different things and both are right.**

### DFS isolation -- `lambdaWidth` is the whole effect

One change at a time off the current driver, 12 controls, segmentation reused
from each subject's original `.mat` so only the fit differs:

| arm | H6/H2 | H4/H2 | MAD | NADH4 vs July |
|---|---|---|---|---|
| `second` (Jul 2026 ref) | 1.160 | 0.953 | 0.185 | — |
| current driver, prior off | 1.408 | 0.572 | **0.418** | **-54%** |
| + `lambdaWidth = 0.3` | 1.191 | **1.025** | **0.226** | +10% |
| + `coarseMode = perRegion` | 1.380 | 0.583 | 0.399 | -53% |
| + region 2 `[9.8 10.45]` | 1.408 | 0.572 | 0.418 | -54% |
| + `coarseFirstOrder = true` | 1.427 | 0.538 | 0.444 | -59% |
| SpecTickle (quoted) | 1.110 | 1.011 | 0.149 | — |

The region-2 edge has **exactly zero** effect on NAD+ (it only moves Trp).
`perRegion` is within noise. Restoring the fitted `pc1` is slightly *worse*,
which independently supports section 12's `coarseFirstOrder = false`.

### Why Trp is exempt

The HSVD width prior for Trp is present in only 91% (U01) / 75% (Pilot-Sinc) of
scans, so the clamped width is **bimodal** -- ~60 Hz when HSVD finds Trp, ~40 Hz
when it does not -- while the free width is one tight distribution. Pinning
therefore *degrades* Trp width stability (U01 25.4 -> 10.9, Pilot-Sinc 24.0 ->
7.2 group CV when freed) and halves reported Trp. NAD+ priors are available in
100% of scans and are stable, so they act as a consistent anchor.

### The cost that is being accepted

`lambdaWidth >= 0.1` is a **hard constraint, not a soft prior** -- the fitted
width equals the prior to within 0% at every lambda from 0.1 to 3. So NAD+
linewidths and, through them, areas are HSVD's, not the data's. HSVD has no
baseline model, so it absorbs the pedestal under each peak into that peak's
width and runs biased high: NADH2 +21%/+9%, NADH6 +2%/-1%, NADH4 +22%/+31%
(U01/Pilot-Sinc). Concentrations sit 5-23% higher than a free fit.

**This is an accuracy cost accepted for a precision gain, and no metric here can
say whether it is correct** -- there is no gold standard. What is measured:
freeing the widths degrades DFS self-consistency 2.3x, and self-consistency is
the closest available proxy for accuracy.

Also ruled out: reducing baseline flexibility does **not** recover free-width
reproducibility. A 2-D sweep (ED/ppm 2/3/4/6 x lambdaWidth 0/0.3, 19 sessions)
found free-width group CV flat across ED and *worse* at every ED than the plain
knotSpacing baseline (U01 NADH6 ~30-44% vs 18.6%). ED is a global DOF count, not
local stiffness at peak scale.

### DFS specifics

No segmentation TSVs ship with that dataset; `reprocess_second_fit.m` reuses
`segFracs` from each subject's original `*_brain_1H_abs.mat` through a synthetic
TSV. That harness asserts **exactly one** non-generated `*_brain_1H_abs.mat` per
subject, so writing a new output stem into the subject folders breaks future
runs unless the stem is added to its `GENERATED` list. Scratch runs should
redirect `outDir` instead.

### 13a. Why the three cohorts report different concentrations

Median under the current driver (`ks_lw0.3` / `lw03trp`), all four peaks:

| cohort | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01_2026 | 583.4 | 683.4 | 545.7 | 163.1 |
| Pilot-Sinc | 325.2 | 406.4 | 352.8 | 60.8 |
| DFS | 409.4 | 498.0 | 397.5 | 114.6 |

`conc ~ ratio x wTerm / R_M`. Decomposed:

| U01 / Pilot-Sinc | conc | = ratio | x wTerm | / R_M |
|---|---|---|---|---|
| NADH2 | 1.79 | **1.44** | 1.05 | 1.09 |
| NADH6 | 1.68 | **1.36** | 1.05 | 1.09 |
| NADH4 | 1.55 | **1.24** | 1.05 | 1.12 |
| **Trp** | **2.68** | **2.34** | 1.05 | **1.00** |

**It is the raw met/water amplitude ratio — the measurement, upstream of every
correction.** Water protocol is identical in all three (`TR_wat = 10000`,
`TE_wat = 13`, matched Tx voltage, same sequence family), and relaxation is
correctly applied and far too small to matter.

Confirmed within-subject: several subjects appear in more than one cohort.

| subj | U01 | Pilot-Sinc | DFS |
|---|---|---|---|
| S092 NADH2 | 627, 742, 526, 496 | 340, 353 | 406 |
| S092 Trp | 161, 162, 365, 215 | 68, 66 | 91 |
| S095 NADH2 / Trp | 577/181, 539/155, 590/270 | — | 328 / 90 |
| S138 NADH2 / Trp | 718/154, 612/150 | — | 415 / 93 |
| S145 NADH2 / Trp | — | 328/49, 280/49 | 453 / 117 |

Two separate effects, do not conflate them:

**1. Pilot-Sinc's segmentation is wrong.** No seg TSVs exist for it, so the
harness assumes a fixed 50 WM / 40 GM / 10 CSF voxel and every subject gets
`fCSF = 0.100`. Measured values for the same subjects are 0.149 (S092) and
**0.328** (S145). S145 is the clean demonstration: raw ratio identical to DFS
(2.72-3.18 vs 2.77) but concentration 280-328 vs 453 uM. **Pilot-Sinc absolute
concentrations are biased low and must not be pooled with the other two.** They
remain valid for test-retest, since both scans of a pair carry the same error.
Note `wTerm` is per-scan, so this bias is metabolite-INDEPENDENT (1.05 for all
four) and cannot explain Trp's extra discrepancy.

**2. U01 measures ~1.6-1.7x more signal in the same subjects**, in the raw
ratio. Cause not yet identified. U01's own session-to-session spread is also
large: S092 gives 496-742 uM NADH2 and **161-365 uM Trp** across four sessions.

### Trp is the cleanest cross-study probe

**`R_M(Trp) = 0.504 in all three cohorts** despite TR 500 / 600 / 1000 ms,
because T1(Trp) = 75 ms is fully relaxed at every TR used. Trp therefore carries
**no relaxation confound at all**, unlike NAD+ whose R_M moves 0.58-0.72 with TR.
Any Trp difference between cohorts is measurement or segmentation, nothing else.

And Trp is the *most* discrepant peak: 2.68x U01/Pilot-Sinc and 1.89x
DFS/Pilot-Sinc, against 1.13-1.79x for the NAD+ peaks. Reporting NAD+ only hides
this. **Always tabulate all four.**

## 14. Aug 2026: Voigt lineshape — reparameterization, and what eta does

> Written while the driver still ran `lambdaWidth = 0.3`; it now runs 0.005
> (section 15). The lineshape findings here are independent of that.

**Outcome:** Voigt modes 5/7 reparameterized to (total width, Gaussian
fraction); a latent pure-Lorentzian bug fixed; `baselineOpt.etaInit` and
`watfit.etaInit` added. Driver runs `etaInit = [0 0 0 1]` and
`watfit.etaInit = 0`.

### A latent bug in the pure-Lorentzian branch

`compositeV_complex_unit` had a `sigma == 0` shortcut that produced a line of
**half** the intended width:

```matlab
gamma = wL/2;  xpw = (x-p)/gamma;      % gamma is a HALF width
Vc = 1 ./ (1 + 1i*2*xpw);              % -> FWHM = wL/2
```

Measured: FWHM 12.00 Hz against Olivero-Longbothum's 24.00, ratio exactly
0.5000. The `1/(1+2i*xpw)` form is correct in modes 2/3/6, where `xpw` is
normalized by the FULL width; the Voigt branch copied it but normalizes by
`gamma`. Fixed by dropping the factor 2. **The branch had never executed** --
the 10 Hz width floor guaranteed `sigma > 0` -- so nothing historical is
affected. After the fix the kernel is continuous across the branch to
**2.2e-16** at `wG = 1e-9`, scaling as `wG^2`, in both Voigt modes.

### The old width floor blocked both pure limits

`widthBoundsHz` was applied to `wL` AND `wG` separately, so `[10 60]` really
constrained the TOTAL to [16.4, 98.3] Hz. Railing at that floor, 30 fits/peak:

| peak | wL railed | wG railed |
|---|---|---|
| NADH2 | 0% | **50%** |
| NADH6 | 0% | **67%** |
| NADH4 | **40%** | 20% |
| Trp | **43%** | 10% |

NADH2/H6 wanted to be *more Lorentzian* than allowed; NADH4/Trp wanted to be
*more Gaussian*. One number blocked both.

### Reparameterization

Nonlinear vector is now `[pos, wV, ph, eta]` (mode 5) and `[pos, wV, eta]`
(mode 7), with `eta = wG/wV` in [0,1]:

```
wG = eta*wV
wL = wV*(1-eta^2)/(A+s),   s = sqrt(A^2 - (A^2-B)*(1-eta^2))
```

The `(1-eta^2)/(A+s)` form avoids the cancellation in the algebraically
equivalent `(A-s)/(A^2-B)`. Verified: `max|OL(wL,wG) - wV|/wV = 1.8e-16`,
eta round-trip `1.1e-16`, kernel FWHM matches `wV` to <=2e-4 across eta 0..1
including both pure limits. Output `pars` layout, `parnames`, `coarseVarpro`
and all analysis scripts are unchanged -- the conversion happens at assembly.
`widthBounds` now bounds the total width, as the docstring always claimed.

Synthetic test (30 realizations, 3 peaks): total-width accuracy unchanged,
eta recovery better, railed parameters **54% -> 2%**, and `cond(J)` unchanged
(3.62e3 vs 3.59e3). **The collinearity is not removed, only relocated into
eta.**

### eta carries the entire ill-determined direction

Same data, old vs new parameterization:

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 `wV` OLD-NEW (Hz) | **0.01** | **0.00** | **-0.02** | 0.50 |
| PS `wV` OLD-NEW (Hz) | **0.00** | **0.01** | **0.01** | 0.61 |
| DFS `wV` OLD-NEW (Hz) | **-0.00** | **0.01** | **-0.01** | 0.27 |

eta moved from ~0.30 to ~0.00 for every NAD peak -- a complete change of
lineshape composition -- while the total width moved by 0.01 Hz and its CV%
did not move at all. **The data pins `wV` and says almost nothing about eta.**

This also explains why `lambdaWidth` behaves as a switch: it regularizes `wV`,
the parameter the data already determines, so it overrides a measurement
rather than resolving an ambiguity (section 13).

### Converged lineshape: NAD is Lorentzian, Trp is Gaussian

Median eta, and the fraction of scans landing at eta < 0.01:

| cohort | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 median | 0.000 | 0.001 | 0.001 | **0.874** |
| PS median | 0.006 | 0.020 | 0.003 | **0.882** |
| DFS median | 0.000 | 0.000 | 0.000 | **0.820** |
| U01 % eta<0.01 | 60 | 60 | 70 | 10 |
| PS % eta<0.01 | 50 | 50 | 62 | 25 |
| DFS % eta<0.01 | 75 | 75 | 75 | 25 |

The distribution is **bimodal, not tight**: a spike at exactly 0 plus a tail
(U01 NADH2 reaches 0.618). The departures track SUBJECT across studies --
S092 shows elevated NADH2 eta in both U01 repeat scans (0.464, 0.618) *and*
both Pilot-Sinc scans (0.245, 0.261).

### METABOLITE eta is seed-sensitive; WATER eta is not

`etaInit` seeds tested on all 30 scans:

- **Trp is bistable.** **6 of 30 scans flip solution branch** on the seed alone
  (e.g. S092_01072025 1.000 -> 0.000; DFS S167 1.000 -> 0.000). Seeding Trp at
  0 costs Pilot-Sinc Trp group CV **22.0% -> 37.2%** and raises Trp
  concentration 10-17%. Seeding at 1 picks the branch consistent with the other
  80% of scans. **This is why `etaInit = [0 0 0 1]`.**
- **NAD is indifferent.** Seeds 0.764 (the historical implicit default,
  from `widthL = 0.4*width`), 0, and 0.25 give identical eta medians and
  concentrations within 1-3 uM. The scans that depart from eta = 0 are the same
  ones each time, so no seed moves them.
- **Water is completely insensitive.** Seeds 0 / 0.25 / 0.5 / 0.75 / 1.0:
  **0 of 30** scans change branch, per-scan eta range median **0.000** and max
  **0.008**, water area moves <=0.10%, concentrations identical to the uM.

**But water eta is highly VARIABLE between scans** -- 0.000 to 0.974 -- with
cohort medians **0.868 (U01) / 0.518 (Pilot-Sinc) / 0.434 (DFS)** at nearly
identical total widths (19.7 / 20.3 / 21.9 Hz). U01's water is markedly more
Gaussian, and that is a robust property of those acquisitions, not a fitting
artifact. Within repeat pairs it reproduces to a median 0.034, but 2 of 6 pairs
differ by 0.281 and 0.322.

So water eta is *well-determined but genuinely variable*; metabolite eta is
*poorly determined*. Same model, opposite behavior -- water has enormous SNR
and one isolated peak, the metabolites are weak, overlapping and competing with
a spline.

### Water eta tracks NADH2 only

Correlation of metabolite eta with water eta (n=30 pooled):

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| Pearson / Spearman | **0.50 / 0.47** | 0.34 / 0.19 | -0.21 / -0.29 | -0.42 / -0.31 |

Consistent sign in all three cohorts for NADH2 only. Split test: median water
eta is **0.487** where NADH2 fits pure Lorentzian vs **0.706** where it picks
up a Gaussian component. NADH6 is weak, NADH4 and Trp go the WRONG way -- so a
regularizer pulling all NAD peaks toward water eta would be wrong for H4.
**Not implemented.**

### The cost being paid

Removing the floor made group CV slightly WORSE:

| cohort | gCV H2/H6/H4 OLD -> NEW | MAD OLD -> NEW |
|---|---|---|
| U01 | 13.4/13.2/22.7 -> 18.2/19.1/29.3 | 0.163 -> 0.178 |
| PS | 11.3/11.2/17.2 -> 13.3/12.7/20.3 | 0.211 -> **0.201** |
| DFS | 12.6/15.7/39.5 -> 13.1/16.1/**37.2** | 0.226 -> **0.206** |

**The 10 Hz floor was accidentally acting as a regularizer** -- badly targeted
but effective, clamping the ill-conditioned direction. Kept anyway: the model
is now honest, the bound means what it says, and both pure limits are
reachable. Any future work on reproducibility should regularize **eta**, not
`wV`.

### Trp bistability is not explained by the multiplet hypothesis

Prediction was that if Trp's Gaussian appearance is an unresolved multiplet,
the Lorentzian branch would be a narrow core with the outer lines pushed into
the spline. **Neither holds**: the 3 seed-independent Lorentzian scans give
`wV` = 46.9 / 52.9 / 37.8 Hz against a cohort spread of 37-49 Hz (not
narrower), and baseline-to-peak ratios of 9.67 / 4.08 / 5.58, the *lowest* in
U01 (others 9.8-23.0, less baseline not more). All three are U01, and two
(S092_03032026/03052026) are the same scans carrying anomalous NADH2 eta and
the highest water eta in the set -- whole-scan lineshape outliers rather than
anything Trp-specific. Mechanism remains open.

## 15. Aug 2026: lambdaWidth revisited AFTER the reparameterization -> 0.005

**Supersedes the `lambdaWidth = 0.3` recommendation in sections 13 and 14.**
The driver now runs **`lambdaWidth = 0.005`**.

Sections 13 and 14 were written against the OLD (wL, wG) parameterization with
the 10 Hz component floor. Once that floor was removed (section 14), the
balance moved. Re-run on all 30 scans with the current code (reparameterized
Voigt, `etaInit = [0 0 0 1]`, `watfit.etaInit = 0`):

### Group CV%

| cohort | lambda | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|---|
| U01 | 0 | 15.2 | 17.4 | 39.1 | 46.8 |
| U01 | **0.005** | **15.0** | **15.1** | 29.4 | 46.8 |
| U01 | 0.3 | 19.1 | 18.8 | **29.2** | 46.8 |
| Pilot-Sinc | 0 | 12.5 | 8.3 | 22.5 | 22.0 |
| Pilot-Sinc | **0.005** | **8.6** | **8.0** | **17.9** | 22.0 |
| Pilot-Sinc | 0.3 | 13.2 | 12.8 | 20.1 | 22.0 |
| DFS | 0 | 14.0 | 24.9 | **28.9** | 19.9 |
| DFS | **0.005** | **12.4** | **15.2** | 29.8 | 19.9 |
| DFS | 0.3 | 13.1 | 16.6 | 35.1 | 19.9 |

### Test-retest CV%

| lambda | U01 (2 pairs) H2/H6/H4/Trp | Pilot-Sinc (4 pairs) H2/H6/H4/Trp |
|---|---|---|
| 0 | 13.3 / 16.1 / 32.8 / 25.5 | 12.5 / 10.2 / 23.4 / 6.2 |
| **0.005** | 13.6 / **11.4** / 19.8 / 25.5 | **2.9 / 2.8 / 9.1** / 6.2 |
| 0.3 | **12.5** / 17.7 / **17.9** / 25.5 | 3.9 / 5.9 / 9.1 / 6.2 |

### Self-consistency (MAD) -- the cost

| lambda | U01 | Pilot-Sinc | DFS |
|---|---|---|---|
| 0 | 0.271 | 0.224 | 0.419 |
| **0.005** | 0.246 | 0.218 | 0.232 |
| 0.3 | **0.164** | **0.202** | **0.209** |

### Why 0.005 was chosen

Pilot-Sinc test-retest is the strongest single line of evidence available --
**four pairs**, the most direct repeatability data in the project -- and it
nearly halves at 0.005 (H6 5.9 -> 2.8%, H2 3.9 -> 2.9%). That is the best
repeatability produced anywhere in this work and approaches the archived
1.9 / 3.2% (section 11). Group CV is best at 0.005 for five of nine NAD
peak-cohort combinations; 0.3 wins only U01 NADH4, and that by 0.2 pp.
Concentrations run 5-10% below 0.3, carrying less of HSVD's width bias.

**The cost is self-consistency**, worse in every cohort and worst in U01
(MAD 0.164 -> 0.246). That is one cohort's ratio metric against four pairs of
direct repeatability, and there is no gold standard for the ratio
([[feedback-no-gold-standard-nad]]).

`lambdaWidth = 0` remains the worst option nearly everywhere -- DFS MAD 0.419,
H4/H2 0.612, NADH4 down to 189 uM. The reparameterization did NOT rescue the
no-prior case; the width prior is still doing real work at 0.005.

**Trp is unaffected by lambdaWidth at any value** (46.8 / 22.0 / 19.9 group CV
in all three arms) because `widthPriorPeaks = [true true true false]` exempts
it -- see section 13.

### Trp retested inside the width prior at lambda = 0.005 -- still exempt

Section 13 exempted Trp because its HSVD prior is available in only 75-91% of
scans, making the clamped width bimodal. Retested under the CURRENT code
(reparameterized Voigt, `etaInit = [0 0 0 1]`, lambda = 0.005), with
`widthPriorPeaks = [true true true true]`:

| cohort | arm | Trp conc | Trp gCV% | Trp wV | Trp eta med | Trp width CV% |
|---|---|---|---|---|---|---|
| U01 | exempt | 165 | 46.8 | 42.9 | **0.87** | **11.4** |
| U01 | in prior | 278 | 36.4 | 50.2 | 0.31 | 21.7 |
| Pilot-Sinc | exempt | 55 | **22.0** | 42.4 | **0.99** | **7.0** |
| Pilot-Sinc | in prior | 125 | 44.4 | 54.5 | 0.02 | 23.2 |
| DFS | exempt | 102 | 19.9 | 43.8 | **0.97** | **4.9** |
| DFS | in prior | 197 | 17.3 | 55.5 | 0.16 | 10.6 |

Trp test-retest: U01 **25.5 -> 45.0**, Pilot-Sinc **6.2 -> 23.0** (3.7x worse).
NAD is bit-identical in both arms -- the regions fit independently.

**The prior does not just bias Trp's width, it destroys its lineshape.** eta
collapses from 0.97-0.99 to 0.02-0.16: the HSVD prior wants ~60 Hz, the free
fit wants ~43 Hz, and the only way to reach a wider nominal FWHM while still
matching the observed profile is to abandon the Gaussian branch for a
Lorentzian with heavy wings. Width CV doubles or triples as a result.

So the HSVD width prior is the WRONG TARGET for Trp, not merely an unreliable
one. `widthPriorPeaks = [true true true false]` stands.

## 16. Aug 2026: why U01 reads ~1.5x higher than DFS and Pilot-Sinc

**The three cohorts are not the same experiment.** U01 runs an `M0_NAD` /
`M0_Water` protocol pair; DFS_CONTROL_STUDY and Pilot-Sinc both run the
`sinc2x2` family. The cross-cohort concentration gap is an acquisition
difference in the metabolite-to-water calibration, not biology, not
segmentation, and not the fit. Do not pool absolute concentrations across these
studies. Within-cohort comparisons and test-retest are unaffected.

Investigated on the three subjects scanned in **both** U01 and DFS (S092, S095,
S138 -- 9 U01 sessions vs 3 DFS), all processed from `grid6` under identical
current-code settings (`lambdaWidth = 0.005`). S092 additionally appears in
Pilot-Sinc.

### The protocols

| cohort | metabolite | water reference |
|---|---|---|
| U01 | `svssel_..._M0_NAD_9.7ppm` | `M0_Water` |
| DFS | `svssel_sinc2x2_2x9.7ppm` | `sinc2x2_waterref` |
| Pilot-Sinc | `svssel_sinc2x2_2X9.7ppm` | `WATERREF_2x4_7ppm` |

`output.met.seqname` is the sequence binary and reads `svssel_latest` for both
U01 and DFS -- it does **not** distinguish them. Read `tProtocolName` and
`sWipMemBlock.adFree`.

### The prediction, and it held

If the gap is protocol-driven, the two sinc cohorts must agree with each other.
Segmentation-free met/water amplitude ratio (the fixed 50/40/10 segmentation
makes Pilot-Sinc absolute concentrations unusable for this comparison):

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 / Pilot-Sinc | 1.51 | 1.40 | 1.31 | 2.67 |
| U01 / DFS | 1.51 | 1.53 | 1.64 | 1.86 |
| **Pilot-Sinc / DFS** | **1.01** | **1.09** | 1.25 | 0.70 |

Pilot-Sinc and DFS agree on H2/H6; U01 sits ~1.5x above both. Trp is the
exception in every row -- Pilot-Sinc Trp is separately known to be the least
reliable of the four.

### Ruled out

| checked | result |
|---|---|
| reference voltage, met vs water | matches in **all 30** sessions; the `dat2mat_svsNAD.m:205` warning never fires |
| VOI | 200 x 200 x 40 mm in every session |
| flip angle / TE | 90 deg / 13 ms everywhere |
| TR and relaxation | 500-600 (U01) vs 600 (DFS); `R_M` differs <= 5% |
| tissue composition | `wTerm` differs 1%; U01 has **more** CSF (0.183 vs 0.154), i.e. the wrong direction |
| water reference quality | fit/integral 1.04-1.24 both, equal FWHM, TR 10 s fully relaxed |
| residual water | extrapolated Lorentzian tail at 9.4 ppm = 0.1-0.5% of peak height; measured 6-8 ppm pedestal 0.17 (U01) vs 0.20 (DFS) of H2 |
| fitted linewidths | U01/DFS 0.93-1.02 on all four peaks |
| coil combination | measured, see below. **The old reasoning here was wrong** -- it said the weights come from water and are applied to both scans "so it cannot bias the ratio". They are, but the water reference is ECC'd against *itself* and so combines with zero loss by identity, while the metabolite scan is not. The conclusion survives; the residual is ~3% |
| ECC method (Klose vs Brown) | full pipeline, all 12 DFS, both settings | no bias: medians within 1.6% on H2/H6, water area 0.995. See below |

### Rejected hypothesis: baseline flexibility

DFS acquires 4096 points where U01 acquires 1024 at the same bandwidth (512 vs
128 ms readout), giving 198 vs 50 points in the NAD window. (True for U01 from
Apr 2026 on; the Mar 2026 S092 sessions ran at 4000 Hz / 256 ms and get 99 --
see "U01 is not one protocol" below.) `baselineOpt.lambda`
enters `solve_linear_block` as a fixed `sqrt(lambda)*D` block appended to `npts`
data rows, so the **effective** smoothing scales as 1/npts: DFS reached baseline
ED 4.68 where U01 got 3.95, and the ED-matching lambda is 3.58 (region 1) / 3.45
(region 2), i.e. the 3.96 points ratio.

Swept lambda = 0.25 / 1 / 3.6 / 8 on both cohorts (66 fits, `lamsw`):

| comparison | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 lam=1 (ED 3.95) / DFS lam=3.6 (ED 3.89) | 1.46 | 1.57 | 2.29 | 2.85 |
| U01 lam=0.25 (ED 4.74) / DFS lam=1 (ED 4.68) | 1.43 | 1.54 | 2.12 | 2.98 |
| unmatched, both lam=1 | 1.43 | 1.58 | 2.25 | 2.19 |

**No effect.** DFS concentrations move <10% across a 32x lambda range. With ~6
basis functions over 1.3 ppm the finest spline feature is ~130 Hz, far too
smooth to absorb a 30-40 Hz peak. Do not redo this.

Two things worth keeping from it anyway:

- `baselineOpt.lambda` is the only active penalty **not** normalized against the
  data block. `lambdaWidth` divides by `numel(y)` and `lambdaPrior` scales with
  npts; `lambda` does neither, so it means different things at different digital
  resolutions. Not what caused this, but real.
- `region.fit.fit` is the reconstructed **model**, not the data (peaks +
  baseline - fit is identically zero). The measured spectrum in the window is
  `output.met.fit.spec` on `output.met.fit.ppm`. An earlier SNR calculation in
  `trp_bistable.m` was invalid for exactly this reason.

### Where the difference actually is

Fit-free, peak minus interpolated local background, water-normalized:

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 / DFS | 1.66 | 1.61 | 1.57 | 1.64 |

A **uniform** factor on all four peaks -- a single multiplicative met-vs-water
scale, nothing peak-specific. It is present before any processing: the raw
`|fid(1)|` MET/WAT ratio straight out of the twix files is 1.31x higher in U01.

Note what is *not* different: the raw water-normalized level in the NAD window
is nearly equal between cohorts (total window area U01/DFS = 0.91). U01 peaks
stand out more because its **troughs are deeper** -- less broad background under
the same peak tops. That is why the fit assigns exactly 2.0x the fraction of the
window to peaks in U01 (36.5% vs 18.4%) while the window total matches.

U01 also carries a DC pedestal beyond 10.5 ppm where there is no signal
(water-normalized level 9-13 vs 3-5 in DFS). This inflates the fitted baseline,
not the peaks, and is a separate observation.

### The header difference

**The two protocols play identical RF.** Six pulses, read from the raw protocol
text (`mapVBVD` parses `tName` to NaN, so the struct is not enough):
`burp90`, `slr90`, `rf_sinc`, `ss_rf_pi_ro`, `ss_rf_pi_sl`, `ss_rf_pi_ph`.
Every amplitude, divided by that session's `flReferenceAmplitude`, matches to
six figures across cohorts:

| | burp90 | slr90 | rf_sinc | ss_rf_pi x3 |
|---|---|---|---|---|
| U01 (ref 300) | 1.0645 | 0.30113 | 0.301849 | 1.594550 |
| DFS (ref 287.2) | 1.0645 | 0.30113 | 0.301849 | 1.594551 |

Same `flNominalB0`, same `flGMax`, same baseline string. `alFree` is identical
in all 64 elements across all eight met/water scans. `adFree` is identical in
all 14 elements except positions 2 and 4.

**But the binary is not the same throughout U01.** `tSequenceFileName` is
`%CustomerSeq%\svssel_latest` for DFS, Pilot-Sinc and U01 Mar-Apr, and
`%CustomerSeq%\svssel_20260101` for the six U01 sessions from May 2026 on. The
RF is unchanged across that switch, but a different build can differ in timing,
gradients or ADC handling. And `Config.SequenceString` -- which is what
`output.met.seqname` reports -- reads `svssel_latest` for **all 30 sessions,
including the six that ran the 20260101 binary**. It is useless as a
discriminator twice over. Read `tSequenceFileName` or `tProtocolName`.

Element meanings (from NW, Aug 2026 -- an earlier version of this section had
4 and 5 swapped and drew a wrong conclusion from it):

- `adFree(4)` = **ADC center** (ADC_CENTER)
- `adFree(5)` = **excitation center**
- `adFree(6)` = **excitation bandwidth in ppm** (= 2 everywhere; pulse duration
  is derived from it, so `sinc2x2` and `M0_NAD` play the same 2 ppm pulse)
- `adFree(2)` = minimum gradient rise time. Not a candidate.

| | adFree(2) | adFree(4) = ADC center | adFree(5) = excitation center |
|---|---|---|---|
| U01 metabolite | 50 | 9.7 | **9.7** |
| U01 water ref | 50 | 4.7 | **4.7** |
| DFS metabolite | 40 | **4.7** | **9.7** |
| DFS water ref | 40 | 4.7 | **4.7** |

So **both** cohorts excite the metabolite scan at 9.7 and the water reference at
4.7, with the same pulse. The met/water calibration is constructed identically
in both. The only sequence-level difference is where the receiver sits.

Pilot-Sinc uses an older, shorter `adFree` and `aRFPULSE` layout and cannot be
compared element-wise with the other two.

### Mechanisms tested and excluded

Everything below was measured on the existing data, not argued from the header.

| hypothesis | test | result |
|---|---|---|
| decimation filter attenuates DFS NAD (off ADC center) but not its water (at center) | \|H(f)\|^2 from the metabolite noise floor (FID tail, FFT) | real filter, wrong place. At 8000 Hz the corner is ~2800 Hz; DFS NAD sits at 1240-1600 Hz with power **0.95-0.98** of center. Biases met/water by 1-2.5%, not 43%. Only the 4000 Hz U01 sessions roll off in band (0.77 at 1400 Hz, 0.49 at 1600) and there NAD is at the center, so it attenuates residual water instead |
| water reference averages do not add coherently (`lPhaseCyclingType` differs) | `\|sum_a x_a\| / sum_a \|x_a\|` at the water peak | **0.989-0.999** in both cohorts, and for both pc=4 and pc=16. `mean()` over averages is correct. No cancellation |
| the water fit is the problem (2.6 pts across the FWHM) | replaced `areaTotal` with a model-free trapz of the water peak | agrees with the fit to 4-9% in both cohorts. Redone on the model-free denominator, U01/DFS = **1.44 / 1.45 / 1.42 / 2.00** -- unchanged |
| `aRFPULSE{7}` looked like a protocol constant that might discriminate | read the raw protocol text | **it is not a pulse.** Index [6] carries `flAmplitudeNL` and nothing else -- no `tName`, no `flAmplitude`, no `bAmplitudeValid`, unlike all six real entries -- in a 256-element `ParamArray`. Its value is absolute, not proportional to the session transmit reference (391.519 across six U01 sessions with refs 282-316; 356.428 across twelve DFS with refs 271-311), so it is a stale value saved into the protocol file. Pilot-Sinc's older protocol has no [6] slot at all. Useful only as a build fingerprint |

### Coil combination: the one asymmetry that is real, and it is ~3%

**Census (all 30 sessions, metabolite and water scans both):** U01 is
**1 channel** in every scan -- single-channel dual-tuned coil, no combination
performed. Pilot-Sinc and DFS_CONTROL_STUDY are **32 channels** in every scan.
No exceptions either way.

So coil combination happens for the sinc cohorts and never for U01, which makes
it the one processing step that is structurally asymmetric between them.

Pipeline order is Klose ECC (`eccopt=0`, `dat2mat_svsNAD.m:321`) and *then*
`coilCombinationNoPC` (`:375`). ECC divides out the reference phase per channel,
so real positive weights `|fid_ref(1,c)|^2` are the correct thing afterwards --
an earlier reading of `coilCombinationNoPC` in isolation called this a missing
phase correction, which was wrong.

The consequence that matters: ECC uses the water reference as its own reference,
so **the combined water is `|nws|` -- real, positive, perfectly coherent by
identity.** It cannot lose anything. The metabolite scan is corrected with the
*water's* phase, so it combines coherently only insofar as the coil phases agree
between the two acquisitions. Any shortfall is met-side only, i.e. a met/water
bias present in the 32-channel cohorts and impossible in U01.

Measured on the metabolite FID first point (high SNR, so effectively unbiased),
coherence `|sum w.v| / sum w|v|`:

| | DFS S092 | S138 | S095 | S179 | S002 | PS S092 |
|---|---|---|---|---|---|---|
| coherence | 0.973 | 0.974 | 0.970 | 0.976 | 0.975 | 0.882 |
| weighted rms phase spread | 13.6 deg | 13.2 | 14.2 | 12.6 | 13.1 | 29.0 |

**~2.5-3% loss for DFS.** Right sign -- it makes DFS read low against U01 -- but
it accounts for about a fourteenth of the 1.43 gap. Not the explanation. The
Pilot-Sinc figure is worse (12%), though that file is the 10.1 ppm protocol
rather than the one `grid_run` processes, so treat it as indicative only.

Two methodological traps here, both worth not repeating:

- **Do not test coherence on the water reference.** It is ECC'd against itself,
  so every channel is real and positive and the coherence is exactly 1.0000 as
  an algebraic identity, not a measurement. Same reason its combined magnitude
  equals `sum w|v|` exactly -- neither the plain sum of channel magnitudes nor
  the RSS, because the weights are neither 1 nor `|S|`.
- **Do not read coherence at the NAD peak.** Per-channel SNR there is 1.2-1.7,
  and `|v_c|` is Rician-inflated while the coherent sum is not, so the raw
  estimator reads ~0.9 under perfect alignment. Noise-debiasing it at that SNR
  overshoots to >1 (1.25-1.33 observed), i.e. unusable in either direction.

Separately, `w ∝ |S|^2` is not the SNR-optimal weighting -- after phase
alignment the matched filter is `w ∝ |S|`, giving `SNR = sqrt(sum SNR_c^2)`.
The pipeline over-weights strong channels. This costs precision in the sinc
cohorts and cancels exactly in met/water, so it does not touch concentrations.

### ECC method: Klose vs Brown, full pipeline, all 12 DFS sessions

Run to test whether the two ECC methods produce a **biased** difference in the
fitted results. `eccopt = 0` (Klose, current) divides out the reference phase at
every time point; `eccopt = -1` (Brown) applies only the constant phase from
`fid_ref(1,c)`. Everything else identical (`eccsw` grid, 24 fits, 0 failures).

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| Klose median (uM) | 399.7 | 484.1 | 377.2 | 101.8 |
| Brown median (uM) | 406.7 | 513.4 | 341.5 | 94.0 |
| **Brown/Klose, median** | **1.012** | **1.006** | **0.845** | **0.938** |
| Brown/Klose, range across sessions | 0.80-1.76 | 0.90-1.17 | 0.64-1.57 | 0.78-1.84 |
| group CV, Klose | 12.4% | 15.2% | 29.8% | 19.9% |
| group CV, Brown | 30.5% | 18.2% | 32.2% | 22.1% |

Water `areaTotal`, Brown/Klose: median 0.995 (range 0.847-1.060).

**No bias between the methods.** Cohort medians differ by <=1.6% on H2/H6 and
the water area -- the denominator of every concentration -- is unchanged. ECC
method cannot account for the cohort gap.

**Klose is still clearly better.** Per-session ratios swing 0.64-1.84 while the
medians hold, i.e. Brown adds scatter rather than a shift: group CV rises on all
four peaks, and self-consistency degrades (MAD of H6/H2 0.085 -> 0.108, H4/H2
0.109 -> 0.144, H4/H2 median 0.890 -> 0.821). Keep `eccopt = 0`.

### U01 is not one protocol

Found while checking the above, and it matters beyond this question. At least
three revisions inside U01_2026:

| | S092 Mar | S138 Apr | May 2026 onward | DFS (all) |
|---|---|---|---|---|
| n sessions | 2 | 2 | 6 | 12 |
| subjects | S092 | S138 | S007, S092, S095 | 12 |
| dwell | 250 us | 125 us | 125 us | 125 us |
| bandwidth / readout | 4000 Hz / 256 ms | 8000 Hz / 128 ms | 8000 Hz / 128 ms | 8000 Hz / 512 ms |
| water ref `lPhaseCyclingType` | 4 | 16 | 16 | 4 |
| `tSequenceFileName` | `svssel_latest` | `svssel_latest` | **`svssel_20260101`** | `svssel_latest` |
| `aRFPULSE[6]` (stale slot, build fingerprint only) | 63.782 | 63.782 | 391.519 | 356.428 |
| **NADH2** (uM) | 535 | 691 | 561 | 400 |
| **NADH6** (uM) | 712 | 866 | 683 | 484 |
| **NADH4** (uM) | 712 | 575 | 493 | 377 |
| **Trp** (uM) | 353 | 145 | 164 | 102 |
| ratio to DFS, H2/H6/H4/Trp | 1.34 / 1.47 / 1.89 / **3.47** | 1.73 / 1.79 / 1.52 / 1.42 | 1.40 / 1.41 / 1.31 / 1.61 | -- |

Medians, `lw0.005` grid6, current settings. **Each U01 column is one or two
subjects, so subject and protocol are fully confounded across columns** -- these
are not three estimates of a protocol effect. The only clean contrast is S092,
who was scanned under both dwells:

| S092 | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| Mar 2026, 250 us (n=2) | 535.0 | 711.6 | 712.4 | 352.6 |
| Jul 2026, 125 us (n=2) | 640.5 | 815.8 | 753.1 | 178.4 |
| **Mar / Jul** | 0.84 | 0.87 | 0.95 | **1.98** |
| DFS (n=1) | 404.0 | 509.7 | 331.0 | 78.7 |

Note the U01 ratio to DFS is between 1.3 and 1.9 on the NAD peaks in every
variant, i.e. the cross-cohort gap is present in all three and is not an
artifact of any one of them. Trp is the unstable one: 3.47 in Mar, 1.42-1.61
after.

Consequences:

- **Test-retest is safe** -- every U01 repeat pair falls inside one variant.
- **Trp is not comparable across the dwell change.** Within S092, the only
  subject scanned under both, Mar/Jul = 0.84 / 0.87 / 0.95 / **1.98** for
  H2/H6/H4/Trp. The NAD peaks move <=16%; Trp doubles. Position-dependent, i.e.
  the sec.12 window-clipping / `areasConv` extrapolation behavior on a coarser
  grid (Trp sits at the window edge), and not a water-scale effect, which would
  be uniform. An earlier group-level version of this comparison put NADH4 at
  1.43; that was subject mix, not protocol -- within S092 it is 0.95.
- The water reference is the most under-sampled object in the pipeline: a ~20 Hz
  line on a 7.8 Hz grid is **2.6 points across the FWHM** (4.7 in the 4000 Hz
  sessions) against 12.1 in DFS. It survives -- the fit tracks the model-free
  integral -- but there is no margin.
- Neither study runs preparing scans on the water reference. At TR 10 s that is
  ~0.7%. Worth knowing, not worth fixing.

### Is `svssel_20260101` the same build as `svssel_latest`?

`svssel_latest` is the working copy and is modified over time; `svssel_20260101`
is meant to be its frozen state at 1 Jan 2026. Two tests, both from the headers.

**Structural** -- a binary declares its protocol parameters, so adding, removing
or renaming one changes the key set regardless of the values a session used.
Populated WIP slots:

| | adFree | alFree |
|---|---|---|
| Pilot-Sinc, Jan 2025, `latest` | 1-9 | 1,3,4,5,7,8,10,11,**16**,64 |
| DFS (S092/S138/S179), `latest` | 1-10,12,13,14 | 1,3,4,5,7,8,10,11,**17**,64 |
| U01 Mar + Apr 2026, `latest` | 1-10,12,13,14 | 1,3,4,5,7,8,10,11,**17**,64 |
| U01 Jun + Jul 2026, `20260101` | 1-10,12,13,14 | 1,3,4,5,7,8,10,11,**17**,64 |

Identical across every 2026 file on both binaries. **Pilot-Sinc differs** -- four
fewer `adFree` elements and `alFree` slot 16 where the others use 17 -- which
matters twice: it confirms `svssel_latest` really did change between Jan 2025 and
2026, and it shows the test has the power to see such a change.

(A 66-vs-68 key count in an earlier version of this comparison cross-cut the
binary boundary -- Mar-`latest` and Jun-`20260101` both 66, Apr-`latest` and
Jul-`20260101` both 68 -- so it was never a binary effect. The WIP slots are
identical in all eight 2026 files; the two extra keys are `sVoI.dInPlaneRot` and
`sVoI.sNormal.dSag`, which the protocol writes only for an angulated slab.)

**Functional** -- S138 Apr (`latest`) vs S092 Jul (`20260101`) have matching
protocol inputs except TR:

| input | latest | 20260101 |
|---|---|---|
| lVectorSize / dwell / TE / flip / averages | 512 / 125 us / 13 ms / 90 / 1024 | same |
| phase cycling / dDeltaFrequency / VOI / adFree 2,4,5,6 | 32 / 9.7 / 200x200x40 / 50,9.7,9.7,2 | same |
| **alTR** | 500 ms | **520 ms** |

| sequence-computed output | latest | 20260101 |
|---|---|---|
| `lPreparingScans` | 16 | 16 |
| `dRefSNR` | 5.656854249e10 | 5.656854249e10 |
| `dRefSNR_VOI` | 1.81019336e10 | 1.81019336e10 |
| ADC points returned | 1056 | 1056 |
| `lScanTimeSec` | 520 | 541 |

Both scan times are exactly `(1024 + 16) * TR` -- same formula, and the only
computed field that moves is the one the changed input forces. A full value diff
over `MeasYaps`, session-specific fields filtered, leaves **exactly one differing
key: `alTR`**.

**Verdict: consistent with identical.** Limit: this only sees what the sequence
writes into the protocol. A change to gradient waveforms, RF phase, ADC timing
within the same point count, or the localization logic would be invisible here.
The definitive check is on the host -- checksum the two files in
`C:\MedCom\MriCustomer\seq\` (`certutil -hashfile <name>.dll MD5`, plus the
`.evp`), or diff the source in the IDEA workspace.

**Separate consequence: Pilot-Sinc ran a genuinely different build**, one more
reason not to trust its absolute scale.

### Not the sinc excitation profile question

The settled result (simulated and experimentally tested) is that the sinc
profile does not distort the peaks **relative to each other** across
8.7-10.6 ppm. That still stands and is not reopened here. A uniform scale factor
on all four peaks is fully consistent with a flat profile sitting at a different
absolute level. These are different claims.

### Where this stands

Every mechanism visible in the header has now been measured and excluded. The
sequence is identical, the excitation is identical, the met/water calibration is
constructed identically, the water reference is faithful, and the factor is
~1.43 uniform on H2/H6/H4 with a fit-free numerator and a model-free
denominator. Whatever it is, it is not in the protocol as written.

The two acquisition variables that remain untested (they cannot be separated
retrospectively because no session varies one alone) are:

1. **ADC center of the metabolite scan** -- 9.7 (U01) vs 4.7 (DFS). Excluded as
   a filter-magnitude effect above, but not excluded as anything else.
2. **vectorSize / dwell / readout length** -- 512 @ 128-256 ms vs 2048 @ 512 ms.

NW is running a volunteer with both protocol pairs back to back, then varying
these one at a time. That is the only thing that decides it; nothing in the
existing data can, and the argument that the U01 pair is self-consistent by
construction is an argument from protocol design, not a measurement.

Until then: do not pool absolute concentrations across U01 and the sinc studies.

Scripts: `u01_vs_dfs.m`, `u01_vs_dfs_bl.m`, `u01_vs_dfs_data.m`,
`ratio_vs_ppm.m`, `lam_run.m` / `lam_score.m`, `ed_match_lambda.m`,
`raw_dims.m`, `exc_pulses.m`, `wip_all.m`, `wsup_probe.m`, `phase_u01_dfs.m`,
`wat_diff.m`, `wat_fit_cmp.m`, `wat_modelfree.m`, `filter_probe.m`,
`dwell_src.m`, `phasecycle.m`, `seq_ident.m`, `pulse7.m`.

**Superseded by §17.** The cross-over volunteer was run on 14 Aug 2026 and
answers this: the gap is the **coil**, not the protocol.

---

## 17. Aug 2026: the cross-over volunteer -- it is the COIL, not the protocol

S138, one session, 2026-08-14, `CoilTest\S138_26_08_14\RESEARCH_NWILSON_20260814`.
Both protocols run back to back on the **dual-tuned** coil and again on the
**32ch 1H** coil, slab matched as closely as the two setups allow.

| tag | coil | protocol | met MID | water MID |
|---|---|---|---|---|
| `DFS_DT`   | dual-tuned | DFSControl | 242 | 241 (4-step) |
| `DFS_DT16` | dual-tuned | DFSControl | 242 | 243 (16-step EXOR) |
| `U01_DT`   | dual-tuned | U01        | 245 | 244 (8 avg) |
| `U01_DT16` | dual-tuned | U01        | 245 | 246 (16-step EXOR) |
| `DFS_32`   | 32ch       | DFSControl | 268 | 267 |
| `U01_32`   | 32ch       | U01        | 270 | 269 |

Everything that could confound is matched: TR 600 ms (met) / 10000 ms (water),
TE 13 ms, 200x200x40 mm slab, 90 deg, 8000 Hz, in all six runs. The protocols
differ only where §16 said they do -- `lVectorSize` 2048 vs 512 (readout 512 vs
128 ms), ADC center `adFree(4)` 4.7 vs 9.7, `adFree(2)` 40 vs 50, 768 vs 1024
averages, `svssel_latest` vs `svssel_20260101` -- and these reproduce the real
cohort protocols exactly (checked against S092/S138/S095 U01 and S002/S138/S179
DFS raw headers). Slab position differed between coils by 7 mm Cor / 9 mm Tra.

**Segmentation cannot enter any number in this section.** One subject, one
session, and TR/TE identical across all six runs, so the Gasparovic
tissue/relaxation correction is a common constant and every ratio below is
exactly the ratio of `met.fit.areasConv` to `wat.fit.areaTotal`. A fixed
0.21/0.39/0.40 TSV was used for all six.

### The water reference does not move

| | dual-tuned | 32ch |
|---|---|---|
| `areaTotal` DFS | 3.813 | 22.81 |
| `areaTotal` U01 | 3.770 | 21.64 |
| **U01/DFS** | **0.989** | **0.949** |

Model-free peak height and in-window integral agree. The denominator is not the
story, in either coil.

### What the pipeline reports (uM)

| run | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| `DFS_DT`   | 483.9 | 566.7 | 348.6 | 135.3 |
| `U01_DT`   | 787.0 | 818.1 | 636.0 | 261.5 |
| `DFS_32`   | 330.2 | 382.9 | 238.1 |  82.5 |
| `U01_32`   | 372.0 | 471.0 | 256.5 |  73.6 |

| ratio | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| protocol U01/DFS, dual-tuned | 1.63 | 1.44 | 1.82 | 1.93 |
| protocol U01/DFS, 32ch       | 1.13 | 1.23 | 1.08 | 0.89 |
| coil DT/32ch, DFS protocol   | 1.47 | 1.48 | 1.46 | 1.64 |
| coil DT/32ch, U01 protocol   | 2.12 | 1.74 | 2.48 | 3.55 |

Taken at face value this says the protocol matters enormously on one coil and
barely on the other, which is already a warning.

### The raw spectra say something different

Water-normalized, phased metabolite spectra, no peak model. On the **dual-tuned
coil the two protocols give the same spectrum**: RMS difference 1.8% of peak
height (1.0% smooth, 1.5% peak-like). On 32ch the difference is 8.0%, and 7.2%
of that is a smooth pedestal offset.

Peak metric that does not depend on the spline baseline: put every run on a
common ppm grid, divide by its own water area, subtract an identical 0.62 ppm
Savitzky-Golay trend, take peak heights. Validated as sampling-insensitive by
the truncation control below (0.99-1.02).

| detrended peak height | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| protocol U01/DFS, dual-tuned | 0.89 | 0.94 | 0.96 | 0.62 |
| protocol U01/DFS, 32ch       | 1.04 | 1.04 | 1.07 | 0.73 |
| coil DT/32ch, DFS protocol   | 1.61 | 1.74 | 1.51 | 1.63 |
| coil DT/32ch, U01 protocol   | 1.39 | 1.57 | 1.36 | 1.39 |

- **Protocol has no effect on NAD** -- 0.89 to 1.07 on both coils.
- **Protocol does affect Trp**: U01 reads 0.62-0.73 of DFS, and Trp is broader
  under U01 (model-free FWHM 37 vs 26 Hz, dual-tuned). This is the one real
  acquisition difference found, and it goes the *opposite* way to the cohort gap.
- **Coil is the effect**: 1.4-1.7x more NAD and Trp per unit water on the
  dual-tuned coil, same subject, same session, same sequence.

### Controls

| control | result |
|---|---|
| readout length -- DFS FID truncated 4096 -> 1024 pts, i.e. U01's 128 ms | fit 483.9 -> 507.7 uM (H2); detrended heights 0.99-1.02. **Not the mechanism** |
| coil combination -- 32ch restricted to its single strongest channel | detrended heights 0.85/0.84/0.90/0.69. Combination is not losing metabolite signal; the receive *weighting* does change NAD/water |
| water-ref phase cycle -- 4/8-step vs 16-step EXOR | raw heights within 5% under both protocols |
| ECC via the other water reference | DFS fit moves 1.6%, **U01 fit moves 31% on NADH2** |
| coarse HSVD order K 60 -> 40 / 80 | `DFS_DT` 3-7%, `U01_32` 3-5%, **`U01_DT` 27/9/37/60%** |

### Conclusion 1 -- the cohort gap is the coil

The dual-tuned volume coil reports 1.5-1.7x more NAD per unit water than the
32ch array on the same brain minutes apart. U01 is dual-tuned in all 10
sessions, DFS and Pilot-Sinc are 32ch in all 20 -- so this alone is the right
size and the right direction for the whole cohort gap.

The mechanism is coil weighting over a 200x200x40 mm **slab**, not a voxel. The
array is surface-weighted and so takes a larger share of its water from the
periphery of the brain, which carries no NAD; the volume coil weights the slab
far more centrally. The single-channel control confirms NAD/water genuinely
depends on the weighting: changing it inside one coil moves the ratio 10-15%.
Which compartment supplies the extra water is settled below -- it is **CSF**,
not extracranial tissue.

The slab was 7/9 mm apart between coils and this cannot be separated from
receive weighting here, but it does not explain the result -- the 32ch slab sat
*further from the skull base*, which would raise its NAD/water, not lower it.

This is the same class of problem as [[project-slice-localization-nonbrain]]:
what is in the slab, weighted by how well it is seen. Tissue segmentation of the
slab does not capture it, and no segmentation correction will.

### Conclusion 2 -- the fit adds a protocol-dependent bias on top

On the dual-tuned coil the two protocols deliver the same spectrum to 1.8% and
the pipeline reports +44 to +93%. That number is manufactured by the fit. The
`U01_DT` fit is the least-constrained case in the whole set -- lowest SNR
(single channel) and coarsest sampling (73 points in the fit window vs 289) --
and it is not converging to a unique solution:

- NADH4's per-peak phase rails at the -60 deg bound in every run, and its width
  rails at the upper bound in `U01_DT16`;
- fitted NAD widths inflate to 44-56 Hz against a model-free FWHM of 21-36 Hz,
  with the extra area taken out of the spline baseline;
- two independent perturbations that do not touch the data -- which water
  reference feeds ECC, and the coarse HSVD order -- move it 27-60%, against
  3-7% for the same perturbations on DFS.

`U01_32` is stable (3-5%), so this is not the U01 protocol per se; it is the U01
protocol at dual-tuned SNR. Every U01 cohort session is exactly that combination.

### What this changes

- The gap is **not** the pulse sequence, the ADC center, the readout length, the
  phase cycling, the water reference, ECC, or coil combination. §16 chased all of
  these and correctly excluded them one by one; the remaining two candidates it
  named (ADC center, vectorSize/readout) are now both excluded directly.
- Absolute concentrations still must not be pooled across U01 and the sinc
  studies, but the reason is now known and is a coil/receive-profile
  normalization problem, not an unexplained protocol artifact.
- The U01 dual-tuned fits carry a real stability problem that is worth fixing on
  its own terms, independent of the cohort question. It inflates U01 NAD by a
  further ~1.4-1.9x on top of the coil effect and it is why the U01 numbers are
  seed- and ECC-sensitive.

### Why water normalization does not cancel the coil

Normalizing to water cancels coil sensitivity **pointwise**, not **on average**:

```
S_met / S_wat  =  int s(r) c_met(r) dr / int s(r) c_wat(r) dr
               =  < c_met(r)/c_wat(r) >   weighted by  s(r) * c_wat(r)
```

It is a sensitivity-weighted *mean of the local met/water ratio*. `s(r)` cancels
only where that local ratio is constant, i.e. only if metabolite and water are
co-distributed. They are not: water is in CSF, brain, and extracranial soft
tissue, NAD only in brain parenchyma, so the local ratio runs from 0 in CSF to
its tissue value in parenchyma. Change the weighting, change the answer.

`s(r)` is the **product** of the transmit weighting and the receive profile, and
the MPRAGE below measures the product. Both enter identically for metabolite and
water, so the argument does not need them separated. For these two setups the
difference is most likely receive-dominated -- both transmit with a volume coil
at similar efficiency (refVolt 300 vs 290) -- but that was not measured here; a
B1+ map would separate them. Transmit is not *exactly* a common weight, because
met TR = 600 ms and water TR = 10000 ms saturate differently with flip angle,
but over 45-180 deg that moves met/water by only 4-6%.

The prescription makes this as bad as it can be: a **200x200x40 mm slab**, not a
voxel, so a large peripheral shell is inside it -- exactly where the two coils
differ most.

### Which compartment: CSF, not extracranial tissue

The discriminating measurement is whether the broad **macromolecule pedestal**
moves with the NAD peaks or against them. Read straight off the water-normalized
spectra in peak-free windows, with no fitted baseline involved:

| dual-tuned / 32ch, DFS protocol | 10.30-10.55 | 9.95-10.02 | 9.60-9.80 | 9.45-9.55 |
|---|---|---|---|---|
| broad signal | 1.74 | 1.50 | 1.47 | 1.45 |

against NAD peaks **1.51-1.74** and Trp **1.63**. Everything downfield falls by
the same factor.

If the array's extra water came from scalp, skull marrow or temporalis muscle,
those tissues carry macromolecules -- pedestal/water would rise on the array
while NAD/water fell. It does not. The lipid probe says the same independently:
`lipid/water` is only 1.12x higher on the array, which under a two-compartment
model (`water = W_brain + W_extra`, lipid tracking `W_extra`) accounts for
**1.2-5.0%** of the effect depending on the assumed baseline extracranial
fraction, not 50%.

So the compartment holds water and **none** of NAD, Trp or macromolecules --
a pure-water compartment, invisible to both probes above, which is why they came
back small.

**But not CSF inside the segmented brain.** The forward check below measures
that directly and it contributes 2.8% (dual-tuned) and 0.0% (32ch), because CSF
and cortical GM are co-distributed in depth and any depth weighting scales them
together. Sulcal CSF is not the answer, and neither is the ventricular/sulcal
balance.

What survives is pure water that is **not** co-distributed with brain tissue and
therefore lies outside the brain mask: extra-axial CSF in the cisterns and the
subarachnoid space beyond the HD-BET boundary, and possibly vitreous if the slab
clips the globes. This is an inference by elimination, not a measurement -- the
tools here cannot weigh a compartment the segmentation does not label.

The direct test is a **small localized voxel vs the slab, on one coil** -- it
measures `B` for that coil with no model at all, and the two candidate
mechanisms predict opposite signs: CSF weighting inside the brain predicts the
voxel reads ~3% LOW (a null), while water outside the brain mask predicts the
voxel reads HIGH, since a central voxel excludes that water entirely and the
correction has no term for it. A null result kills both.

The other direct test is a **multi-TE water reference on each coil**: CSF T2 at
7T is ~1-2 s against ~40-50 ms for tissue water, so a second water reference at
TE ~300 ms isolates it, including the extra-axial CSF the segmentation never
labels. That yields a correction rather than a diagnosis. Cheap controls alongside: a homogeneous
phantom on both coils (a coil-dependent met/water *there* would mean something
instrumental is also involved), and a B1+ map. Strong prediction: a small,
centrally placed, genuinely 3D-localized voxel should collapse the difference
toward zero.

### The two profiles, measured

Both MPRAGEs from this session are `ORIGINAL\PRIMARY\M\ND` -- no normalization
filter -- so they carry the coil weighting directly. Same head, same sequence,
minutes apart. Median signal binned by normalized depth (0 = center of head,
1 = surface):

| depth | 0.13 | 0.23 | 0.33 | 0.43 | 0.53 | 0.63 | 0.73 | 0.83 | 0.93 |
|---|---|---|---|---|---|---|---|---|---|
| dual-tuned | 1.00 | 0.89 | 0.73 | 0.48 | 0.37 | 0.29 | 0.26 | 0.25 | 0.23 |
| 32ch       | 0.86 | 0.98 | 0.90 | 0.85 | 0.88 | 0.82 | 0.76 | 0.83 | 0.86 |
| **32/DT**  | 0.86 | 1.10 | 1.24 | 1.79 | 2.35 | 2.81 | 2.91 | 3.34 | 3.80 |

Surface(>=0.8) / center(<=0.5) signal ratio: dual-tuned **0.33**, 32ch **0.97**
-- the array weights the periphery **3.0x** more relative to the center. It is
obvious by eye: the dual-tuned image is center-bright with a nearly invisible
scalp (7T volume-coil center brightening), the 32ch has a bright scalp rim and a
dimmer brain.

### Forward check -- and what it rules OUT

**Two earlier attempts at this were wrong; both are withdrawn.**

1. A radial `rB = 0.75` threshold gave a "41.6% geometric non-brain share". That
   is just `1 - 0.75^2 = 43.75%` for a disc -- set by an arbitrary threshold,
   measuring nothing anatomical.
2. Replacing it with an Otsu-thresholded head mask was worse in a subtler way:
   **CSF median T1 is 143 against an Otsu threshold of 148**, so the mask
   excluded ~50% of all CSF voxels and zeroed their weight. That depressed the
   apparent CSF weighting for *both* coils by nearly the same amount and made
   the effect look real when it was not.

Redone with a per-slice convex-hull head mask (cannot drop dark tissue; weight
coverage 100.0% on all three classes) the answer inverts:

| dual-tuned | CSF | GM | WM |
|---|---|---|---|
| geometric | 0.2082 | 0.4062 | 0.3856 |
| sensitivity-weighted | 0.1938 | 0.3995 | 0.4067 |

| 32ch | CSF | GM | WM |
|---|---|---|---|
| geometric | 0.2184 | 0.4180 | 0.3636 |
| sensitivity-weighted | 0.2182 | 0.4183 | 0.3635 |

Bias factor `B = ((1-f_eff)/(1-f_geo)) * (sum f_geo c R)/(sum f_eff c R)`, the
error from feeding geometric fractions to Gasparovic when the measurement
carries sensitivity-weighted ones:

**B = 1.028 (dual-tuned), 1.000 (32ch).**

**So CSF weighting INSIDE the brain mask is not the mechanism.** It accounts for
2.8% of a ~46% coil difference. 44.8% of the slab's CSF sits within 5 mm of the
brain surface, so this is not for lack of peripheral CSF to reweight.

The reason is worth stating, because it is general: **within the brain, CSF and
cortical GM are co-distributed in depth.** Any depth-dependent weighting scales
both together, so `f_CSF` barely moves, and ventricular CSF being centrally
over-weighted roughly cancels sulcal CSF being peripherally under-weighted. The
met/water ratio is close to insensitive to the coil profile *inside* the brain.

That is a real constraint on the mechanism: the water responsible must **not be
co-distributed with brain tissue**, i.e. it lies outside the segmented brain
entirely -- extra-axial CSF in the cisterns and the subarachnoid space beyond
the HD-BET boundary, and possibly vitreous if the slab clips the globes. All of
it is peripheral, water-dense, macromolecule-free and lipid-free, which is
exactly the signature the pedestal and lipid probes point to.

**No forward model here is trustworthy** and no number from one is quoted any
more. Every version has been dominated by how the head/brain masks and the
MPRAGE-intensity depth proxy behave, not by the physics. The slab-vs-voxel
experiment measures `B` directly and does not need any of this.

### How similar were the two slabs? (FSL segmentation, both MPRAGEs)

Full standard path, run separately per coil with the MPRAGE and SVS series
selected explicitly: `dicom2nifti` -> SPM8 bias -> `imscribe_engine` svsmask ->
HD-BET on the host -> `brainseg.sh -f -b N -S` (FAST 3-class + FIRST) in the VM.

Inside the SVS slab:

| | f_CSF | f_GM | f_WM | GM/WM | brain volume in slab |
|---|---|---|---|---|---|
| dual-tuned | **0.2082** | 0.4062 | 0.3856 | 1.053 | 377.2 mL |
| 32ch | **0.2184** | 0.4180 | 0.3636 | 1.150 | 323.6 mL |

**CSF differs by 1.0 percentage point, 4.9% relative -- not notably different.**
For scale, the two whole-brain segmentations of the same head differ by 5.2% on
CSF (255.2 vs 242.5 mL) while agreeing to 0.7% on GM, 0.2% on WM and 0.7% on
total brain, so the slab CSF difference is the same size as the method's own
reproducibility on this subject. It is noise, not anatomy.

What the 7/9 mm slab offset *does* show up in is the amount of brain caught:
377.2 vs 323.6 mL, **16.5% less on the 32ch**, whose slab sat 9 mm higher where
the brain cross-section is smaller. The GM/WM ratio shifts the same way (1.053
-> 1.150, more cortex and less deep WM in a more superior slab).

Re-running the pipeline with these real fractions in place of the fixed
0.21/0.39/0.40 TSV moves the dual-tuned results by **-0.1%** and the 32ch by
**+1.9%**, uniformly across all four peaks. The coil ratio goes 1.466 -> 1.436
(NADH2, DFS protocol). So the segmentation difference between the two
acquisitions accounts for **2% of the ~46%**, and every earlier number in this
section stands to within 2%.

Note what the TSV does and does not count: it is the FAST/FIRST segmentation
inside the **HD-BET brain mask**, intersected with the svsmask. Scalp, skull and
muscle are not in it at all -- which is correct for what Gasparovic models, and
is also why the model has no term for extracranial water.

### Supporting evidence and exclusions

| observation | reading |
|---|---|
| downfield **pedestal**/water in peak-free windows, DT/32ch = 1.45-1.74, vs NAD peaks 1.51-1.74 and Trp 1.63 | the whole tissue-borne downfield spectrum moves together. Water is being over-counted; NAD is not being selectively lost. This is also what identifies the compartment as CSF -- see above |
| lipid(1.3 ppm)/water in the water reference, 32ch/DT = **1.12x** | accounts for only 1.2-5.0% of the effect. Rules extracranial tissue **out** as the driver. Weak in any case because 1.3 ppm is 3.4 ppm off the 2 ppm excitation center, so that slab is displaced ~68 mm |
| water-line wing fraction, 32ch/DT = **1.10x**, water FWHM 0.085 vs 0.078 ppm (24.6 vs 22.3 Hz) | array's water comes from a less well-shimmed, more distributed volume. Same size as the lipid probe, and likewise far too small to be the driver |
| ADC clipping on the water reference | none: max\|FID\| 0.249 (DT) / 0.199 (32ch), 1-5 samples near max, no flat-topping |
| coil-combination weight exponent p swept 0 -> 8 (`weights ~ \|fid_ref(1,c)\|^p`) | NAD/water moves only ~7%, Trp ~20%. Reweighting 32 *surface* elements cannot synthesize a volume-coil profile, so the small span is expected and does not argue against the mechanism |
| chemical-shift displacement of the metabolite slab (2 ppm bandwidth / 40 mm = 20 mm per ppm; H2 -6.4, H6 -10.2, H4 -16.8, Trp +8.0 mm, water reference 0) | **not** the driver: the coil ratios are 1.61/1.74/1.51/1.63, with no ordering by displacement, and Trp is displaced the opposite way yet behaves like the NADs |

### What this means for quantitation

The Gasparovic denominator is built from **geometric** tissue fractions, but the
water signal that goes into it carries **sensitivity-weighted** fractions. The
two agree only when `s(r)` is flat over the voxel. Over this slab, with these two
coils, they do not:

- the dual-tuned coil under-counts peripheral CSF relative to the model's
  assumption, so `S_met/S_wat` is too high and **U01 reads high**;
- the 32ch array counts the periphery about right geometrically, so it carries
  more of the CSF water the correction assumes is fully sampled, and reads low.

The segmentation TSV does supply `f_CSF`, so the model does know about CSF -- but
it supplies the **geometric** fraction, which is right for neither coil. That is
the whole failure: the correction and the measurement disagree about how much
CSF was seen, and they disagree by different amounts on each coil.

Neither cohort is "correct". The dual-tuned number is closer to a brain-only
measurement, for the accidental reason that its profile suppresses the
periphery. See [[project-slice-localization-nonbrain]]: this is the same problem
as counting non-brain volume in the svsmask, one step further back -- it is not
only what is in the slab, it is how strongly each part of it is seen.

The fix is not a segmentation correction. It is either a genuinely localized
voxel, or a water reference whose spatial weighting is matched to the metabolite
signal's, or a measured `s(r)` folded into the water model.

Scripts: `ct_hdr.m`, `ct_run.m`, `ct_score.m`, `ct_plot.m`, `ct_cmp_real.m`,
`ct_trunc.m`, `ct_modelfree.m`, `ct_diff.m`, `ct_detrend.m`, `ct_cc.m`,
`ct_rail.m`, `ct_seed.m`, `ct_gain.m`, `ct_prof.m` / `coilCombinationCT.m`,
`ct_why.m`, `ct_profimg.m`, `ct_predict.m`, `ct_pedestal.m`, `ct_prep_anat.m`,
`ct_run_seg.m`, `ct_onefig.m`, `ct_onefig_wat.m`, `ct_nonbrain.m`, `ct_predict2.m`, `ct_effcsf.m` (superseded), `ct_geom2.m`,
`ct_slabview.m`, `ct_voi.m`.

## 18. Sep 2026: FIRST drops BrStem (new Stage 1c in processU01_Brain.ps1)

`859875_S01_09012026_brain` failed in the VM with a message that pointed
nowhere near the actual problem:

```
fsl_sub: error: argument -j/--jobhold: expected one argument
>>> immv .../t1w_SPM8bias_hdbet_first_all_fast_firstseg .../t1w_SPM8bias_hdbet_first
Could not find a supported file with prefix "..._all_fast_firstseg"
Exiting brainseg.sh
```

**The chain, bottom-up.** Task 15 of the `first` array is BrStem, and it
segfaulted:

```
first.e3300.15:  Segmentation fault (core dumped)
first.o3300.15:  mode  2.8823e+17
                 Error: cannot find image ..._first-BrStem_first
```

Every other structure reports `mode` ~180-240; BrStem diverged to 2.9e17
before crashing. `FSLPARALLEL` is empty, so fsl_sub 2.7.5 runs the array
through the shell plugin: one task segfaulting makes fsl_sub exit non-zero
and print no job id, so `HOLDID` at `run_first_all:282` comes back empty.
Line 295 is then `fsl_sub ... -j $HOLDID -t ${COM}2` — unquoted and empty,
so argparse sees `-j` followed by `-t`. That aborted job is `first_post`,
the fslmerge that builds `_all_fast_firstseg`, hence the `immv` failure and
`brainseg.sh` exit 1.

**Why BrStem diverged.** The brain mask is clipped at the inferior edge of
the volume. Geometry is identical across sessions (120x120x160, 2x2x1 mm,
z = I-to-S); the nonzero bounding box of `*_hdbet.nii.gz` is not:

| session | z start | z span |
|---|---|---|
| 859875_S01_09012026 | **0** | 135 |
| CoilTest_32 | 14 | 131 |
| CoilTest_DT | 12 | 132 |

The brainstem runs into slice 0 with no inferior margin, so the shape model
has no boundary to converge against. Dataset-specific, not a VM/FSL config
problem — the two CoilTest sessions segmented fine.

**Why dropping it is free.** The svsmask spans z 51-126. Thalamus sits at
z 53-81, hippocampus at 36-64, so the brainstem is entirely below the slab
and contributes nothing to the tissue fractions the pipeline consumes.

**The trap.** `run_first_all` sets `methname=$btype` *after* the structure
loop (line 286-288), so the merged file is named from whatever the **last**
structure left in `btype`. The default list ends on BrStem (bcorr=1 ->
`btype=fast`) -> `_all_fast_firstseg`. Drop BrStem naively and the list ends
on R_Thal (bcorr=0 -> `btype=none`) -> `_all_none_firstseg`, and
`brainseg.sh:334`'s hardcoded `immv ..._all_fast_firstseg` fails in exactly
the same way. So the list must **end on a boundary-corrected structure**
(Accu / Amyg / Caud / Hipp). Stage 1c ends on `R_Hipp`.

Reordering is otherwise safe: label values are intrinsic to each `_corr`
file (L_Thal=10, L_Caud=11, L_Hipp=17, L_Accu=26, R_Thal=49, R_Caud=50,
R_Hipp=53), not assigned by merge order.

**Where the change lives.** `brainseg.sh` is in
`/home/camipm/Documents/Elliott-Scripts/` and is shared with other studies,
so it is left untouched. Stage 1c instead pre-stages the FIRST result the
same way Stage 1b pre-stages HD-BET: run `run_first_all` with the reduced
list, `immv`, write the label legend, clean up intermediates. `brainseg.sh`
then finds `<stem>_hdbet_first.nii.gz` on disk and logs "Skipping FIRST."
The guard at `brainseg.sh:331` is live because the driver passes `-f`
(doFIRST) and not `-x` (forceRERUN).

Verified on 859875_S01: exit 0 in 71 s, merged file correctly named
`_all_fast_firstseg`, all 14 structures populated, label 16 empty.

| label | voxels | | label | voxels |
|---|---|---|---|---|
| 10 L_Thal | 2524 | | 49 R_Thal | 2524 |
| 11 L_Caud | 1203 | | 50 R_Caud | 1126 |
| 12 L_Puta | 1784 | | 51 R_Puta | 1627 |
| 13 L_Pall | 461  | | 52 R_Pall | 596  |
| 16 BrStem | **0** | | 53 R_Hipp | 974 |
| 17 L_Hipp | 931  | | 54 R_Amyg | 89   |
| 18 L_Amyg | 264  | | 58 R_Accu | 42   |
| 26 L_Accu | 72   | | | |

Caveat, not acted on: at 2 mm in-plane FIRST is outside its comfort zone and
the volumes run large (putamen ~7 mL vs a typical 4.5-5.5). Accumbens is
tiny and L/R asymmetric (72 vs 42 voxels). Fine for GM/WM/CSF fractions in
the slab; do not read the per-structure volumes as morphometry.

## 19. Sep 2026: `-Ref31P` — selectable 31P primary reference

**What the drivers did before.** Neither `processU01_Brain.ps1` nor
`processU01_Calf.ps1` exposed any control over the 31P reference. Their
`-ExtRef / -ApplyB1 / -RefPpm / -RefConc` flags configure only the external
phantom path, and even with `-ExtRef` on, the external numbers landed in
`absQ.external` and never became the reported result.

Both MATLAB stages computed *both* internal references via `utils/refQuant31P.m`
into `absQ.byRef.aATP` and `absQ.byRef.PCr`, then hardcoded **aATP** as primary:

| | primary | [ref] | secondary | [ref] |
|---|---|---|---|---|
| brain | aATP | 2.7 mM | PCr | 3.5 mM |
| calf  | aATP | 8.2 mM | PCr | 33 mM |

So the reported `absQ.metabolites` — what the Stage 5 console summary and
`checkU01_2026_relaxation` read — was referenced to aATP, the species §
[alpha-ATP is under-fit] identifies as the weakest of the three.

**The flag.** `-Ref31P <aATP|PCr|External>` on both drivers, threaded to a new
trailing `ref31P` argument on `processU01_Brain_31P` / `processU01_Calf_31P`
(added to the `arguments` block, default `'aATP'`). It selects only which set
populates `absQ.reference` / `absQ.metabolites`; every available set is still
computed and saved under `absQ.byRef` / `absQ.external`. The choice is recorded
in the new `absQ.primaryRef`.

Notes on the implementation:

- `-Ref31P External` implies `-ExtRef` and turns it on with a console notice —
  selecting the external reference without adding `RefExt` to the basis is
  incoherent.
- External is promoted *after* `extRefQuant31P` runs, so the internal set is
  still the basis for its `ext_over_int` column.
- If the external quant failed (`absQ.external.ok = false`) and External was
  requested, this is a hard error, not a fall back. The sets differ by a scale
  factor, so silently substituting an internal one would mislabel every
  concentration.
- `absQ.metabolites` from the external path carries two extra fields
  (`conc_internal_mM`, `ext_over_int`). `collectU01_2026_TSV` reads
  `absQ.byRef.*` and `absQ.external` directly and is unaffected by the primary.

**Verified** on `S092_03032026_31P-1H-HeadCoil` (output to scratch, brain stage
run twice):

| requested | primaryRef | NADplus | NADH | PCr | aATP |
|---|---|---|---|---|---|
| aATP | aATP | 0.3300 | 0.3275 | 4.788 | 2.700 |
| PCr  | PCr  | 0.2412 | 0.2394 | 3.500 | 1.973 |

`byRef` holds both sets in both runs. The two columns differ by one common
factor of 1.368 (4.788/3.5 = 2.700/1.973 = 0.3300/0.2412), and NADplus/NADH is
1.0076 either way — the documented behavior: changing reference rescales
everything and changes no ratio.

That 1.368 *is* the aATP/PCr disagreement: referencing to aATP = 2.7 makes PCr
read 4.79 mM against an assumed 3.5, or equivalently referencing to PCr = 3.5
makes aATP read 1.97 against an assumed 2.7. Consistent with the aATP under-fit.

**Default is still `aATP`**, i.e. current behavior is unchanged. Flipping it to
PCr is a one-word edit in both `param()` blocks and both `arguments` blocks.

**Not yet exercised: `-Ref31P External`.** The guards are in place but no
end-to-end run has been done, and per § [Ext-ref volumes are the weak link] the
brain V_M and calf phantom volumes are the known-bad part of that path.

## 20. Sep 2026: external-reference B1+ from the 31P B1CALIB series

Replaces the 1H TurboFLASH B1 map as the source of the transmit correction in
`extRefQuant31P`. The map was the right shape but the wrong nucleus, and the
faint D2O/PBS phantom is masked out by `calc_B1map_new` so its relB1 had to be
extrapolated — which is why `applyB1` defaulted to false. The new source is
the on-scanner ISIS_SVS `*B1CALIB*` flip series: 31P, same coil, same session,
measured at the actual peak positions.

**Reused `isis_svs_refVolt` rather than writing a parallel program.** It
already reads spectroscopy DICOM (`readSpecDicom`), fits relB1 **per peak**
(`output.relB1Peak`), takes a configurable `pars.peaks` list, lets
`pars.combine` keep refVolt on PCr, and models the off-resonance rect
excitation profile (`rectMxy`) — which the +8 ppm reference needs far more
than PCr does. A second program would have duplicated the DICOM reader, the
HSVD pole extraction, the `dropPts` pedestal handling and the fit, then drifted.

### What was added

| file | change |
|---|---|
| `utils/rectMxy.m` | **new** — extracted from `isis_svs_refVolt` so both callers share one model. Verified byte-identical results after the move (refVolt 359.7552, relB1 0.8339 / 0.5720). |
| `utils/relB1FromCalib31P.m` | **new** — locates `*B1CALIB*` under `SessionDir\E*`, calls `isis_svs_refVolt` with `{PCr, RefExt}`, returns per-location relB1 + diagnostics. |
| `utils/extRefQuant31P.m` | `extRef.b1calib` supersedes the map for transmit only; excitation factor now per metabolite via `rectMxy` at the **acquisition's** geometry. |
| `processU01_Brain_31P.m` / `_Calf_31P.m` | `b1Source` / `b1CalibFiles` defaults; read acquisition pulse/center/Hz-per-ppm from the header; pass `basisInfo` for per-metabolite shifts. |

Only the tissue peak needs measuring for the tissue side — PCr and every other
metabolite share the ISIS slab, so one relB1 covers them all. Only the phantom,
elsewhere in the coil, needs its own.

### Two things that are easy to get wrong

**The calibration's geometry is not the acquisition's.** On `859875_S01`
the calibration is a 600 µs pulse at 0 ppm; the acquisition is 200 µs at
−4.7 ppm, flip 57°. The long calibration pulse is deliberate — it is how flips
past 90° are reached — and the narrow profile that comes with it is the price.
So the excitation factor is recomputed at the acquisition's own pulse, center
and flip, never carried over from the fit.

**`isis_svs_refVolt` rmpaths `utils` and `mapVBVD` on exit** (line 615-616).
Calling it mid-pipeline therefore stripped `utils` from under the driver and
the next `utils` call failed with *"Undefined function 'extRefQuant31P' for
input arguments of type 'cell'"*. `relB1FromCalib31P` now snapshots and
restores the path with `onCleanup`.

### Size of the effect

859875 acquisition geometry (200 µs, −4.7 ppm, 57°), relB1 0.8339 tissue /
0.5720 phantom. `fac = eta_ref/eta_met`, old = on-resonance sine:

| met | ppm | off Hz | eta_new | fac_new | fac_old | change |
|---|---|---|---|---|---|---|
| PCr | 0.000 | 565 | 0.7250 | 0.6410 | 0.7305 | −12.2% |
| gATP | −2.530 | 261 | 0.7350 | 0.6323 | 0.7305 | −13.4% |
| aATP | −7.560 | −344 | 0.7330 | 0.6340 | 0.7305 | −13.2% |
| NADH | −8.111 | −410 | 0.7310 | 0.6357 | 0.7305 | −13.0% |
| NADplus | −8.299 | −433 | 0.7302 | 0.6364 | 0.7305 | −12.9% |
| bATP | −16.180 | −1381 | 0.6617 | 0.7023 | 0.7305 | −3.9% |

Almost all of the ~13% is the **reference** peak's own profile (eta_ref
0.5388 → 0.4647): at +8 ppm it sits 1528 Hz off a −4.7 ppm center. The tissue
peaks are all within ±565 Hz and vary only ~1% among themselves, so the
per-metabolite refinement matters mainly for bATP. Do not read this as a 13%
"improvement" — it is a term that was previously missing, not a measured gain.

### Status: not yet exercised end to end

**No session has both inputs.** The volumes still come from the 1H B1 map;
across all 37 U01_RR sessions plus 859875, only `S095_07232025` has any B1 map
(and it is `PREP_MOCO_B1MAP_*`, not the `PREP_TFL3D_B1MAP` that `findB1Dir`
prefers), and only `859875_S01_09012026` has B1CALIB. So the concentration
arithmetic with the new eta has not run on real data.

Verified so far on `859875_S01`: the calibration fits (relB1 0.8339 / 0.5720,
refVolt 359.8 V), the acquisition geometry is read correctly (200 µs,
−4.7 ppm, 120.298 Hz/ppm), `basisInfo` arrives with 21 components, the path is
preserved, and the external quant stops cleanly at *"no B1 map series found in
session"* while the internal quant is untouched (NADplus 0.6663 mM).

`applyB1` is left **false**. The objection that made it false is gone, but
flipping the default before a single end-to-end run would change every external
number on trust. `-ApplyB1` on both drivers turns it on explicitly.

Also unresolved, and now the binding constraint: § [Ext-ref volumes are the
weak link]. The transmit term is on firmer ground than the volumes it
multiplies.

### 20a. Revisions

**Phantom volume is geometric, and the fallback is now loud.** V_ref is the
known cross-section times the ISIS slab thickness read from the header:

```
V_ref = refArea_mm2 * slabThickness_mm
```

`refArea_mm2 = 200` (measured off the images; observed 180-220, so ~±10%
straight onto every external concentration). The phantom is a 15 mL tube
standing across the slab. `slabThickness_mm` is per-scan and NOT constant
(40 mm brain, 60 mm calf), so it must never be hard-coded.

The threshold-segmentation fallback recovers only 64-78 mm² against the
measured 200 -- the D2O/PBS fill is faint in the 1H magnitude image and the rim
is lost to partial volume -- so it is wrong by ~3x, not by a little, and biases
every external concentration HIGH. It now emits a warning naming the missing
field instead of quietly returning plausible-looking numbers. Kept as
`V_ref_seg_mL` for QC only.

Exposed as `-RefArea` on both drivers (default 200) so the cross-section is a
real input rather than a buried constant.

**Never assume the calibration and the acquisition share parameters.** They are
not expected to match, and nothing is wrong when they do not -- the two scans
are tuned for different jobs:

| | calibration | acquisition |
|---|---|---|
| goal | reach flips past 90° so the sine curve is sampled through its maximum | excite the whole spectrum uniformly |
| pulse | long (narrow spectral profile is the accepted cost) | short (flat profile) |
| 859875 | 600 µs at 0 ppm, flips 50-130°, TR 30 s, ref ampl 300 V | 200 µs at −4.7 ppm, flip 57°, TR 5 s, ref ampl 331 V |

So take from the calibration only what is a property of the **coil and the
location** -- relB1 -- and recompute every sequence-dependent term with the
acquisition's own header values. `b1.pulseDur` / `b1.excCtr` are returned for
QC and provenance, never for reuse in the quant. Written into the
`relB1FromCalib31P` header where the numbers are produced.

**B1CALIB locating moved out of the drivers into `extRefQuant31P`,** since it
is only ever used for absolute quantitation. `extRefDefaults` now sets
`e.sessionDir`; `extRefQuant31P` does the search and fit itself. A caller can
still pre-supply `extRef.b1calib` (already fitted) or `extRef.b1CalibFiles`
(explicit DICOMs) and both win over the search.

It runs BEFORE the B1-map early returns, so `absExt.b1calib` is recorded even
when the volume path cannot run -- verified on 859875: external `ok=0`
("no B1 map series found in session") while `b1calib.ok=1` with relB1
0.8339 / 0.5720 still saved.

### 20b. The 1H B1 map is now optional

Nothing the external quant needs comes from it any more:

| term | source |
|---|---|
| relB1 (tissue + phantom) | 31P B1CALIB flip series |
| V_ref | `refArea_mm2 * slabThickness_mm`, geometric |
| V_M | explicit `brainVol_mL`, or HD-BET (brain) / Dixon (calf) mask INTERSECT the ISIS slab |

The map is still read when present — its threshold masks remain the last-resort
V_M and give QC voxel counts — but the hard early return on a missing `b1Dir`
is gone, so a session without one runs through to concentrations.

What changed structurally:

- The map read + segmentation is wrapped in `if haveMap`, with every variable
  it would have supplied pre-set (`voxVol_mm3 = NaN`, so any accidental
  voxel-count route yields NaN rather than a quietly wrong number).
- The old "no phantom / no tissue voxels" bail is **replaced by a check on the
  final V_ref and V_M**. The masks are only one of several routes now and no
  longer the expected one, so validating the volumes is what actually matters.
  With no map and no anatomic mask, V_M is NaN and the quant stops naming
  exactly which input to supply.
- `applyB1` with no relB1 from either source now stops rather than multiplying
  every concentration by NaN.
- The map-based relB1 estimator and its "falling back to the 1H map" warning
  are suppressed when there is no map to fall back to.

**First end-to-end external-reference run** (859875_S01, no B1 map,
`applyB1=true`, `-Ref31P External`):

```
V_ref    = 8.000 mL   [cross-section 200 mm2 * 40 mm slab (hdr)]
V_M      = 1100.0 mL  [caller brainVol_mL -- STAND-IN, not measured]
relB1    = 0.8339 tissue / 0.5720 phantom   [31P B1CALIB, 600 us at 0 ppm]
excModel = rect 200 us at -4.7 ppm, flip 57 deg
eta_ref/eta_M = 0.641 [APPLIED]   recv = 0.686 [APPLIED]
```

| met | ext [mM] | int [mM] | ext/int |
|---|---|---|---|
| NADplus | 0.8938 | 0.6663 | 1.34 |
| NADH | 0.2134 | 0.1593 | 1.34 |
| PCr | 7.458 | 5.52 | 1.35 |
| aATP | 3.608 | 2.7 | 1.34 |

**These are not real concentrations.** V_M = 1100 mL was passed in by hand to
exercise the path; 859875 has no HD-BET mask staged for the 31P slab. The
constant 1.34 across every metabolite is the expected signature of a single
scale factor, so what this validates is the plumbing, not the number. A real
V_M moves all four together.

Note the two geometries printed separately and correctly: the correction uses
the acquisition's `rect 200 us at -4.7 ppm, flip 57`, while `b1method` records
that relB1 came from a `600 us at 0 ppm` calibration — the mismatch is expected
(§20a) and is now visible in the output rather than assumed away.

### 20c. Legacy B1 mapping removed outright

The 1H TurboFLASH B1 map path is gone, not merely optional. Every input to the
external quant now has exactly one source and no fallback:

| input | source |
|---|---|
| relB1 (tissue + phantom) | 31P ISIS B1CALIB flip series |
| V_ref | `refArea_mm2 * slabThickness_mm` |
| V_M | explicit `brainVol_mL`, or HD-BET (brain) / Dixon (calf) mask ∩ ISIS slab |

A missing input stops the quant naming what to supply, instead of substituting
something that yields a plausible wrong number. That was the whole problem with
the removed code: every one of its routes failed silently.

Deleted from `utils/extRefQuant31P.m` (781 → 416 lines):

- `readB1` call, `B1map` / `images` / `b1alpha` / magnitude volume
- slab-slice selection (`dropEndSlices`, `sel`, `slthk`, `voxVol_mm3`)
- tissue and phantom threshold segmentation, and the QC overlay figure
- local helpers `b1VoxelDims`, `keepLargestCC`, `imdilate3`,
  `computeTissueMask`, `extrapB1ToMask`, `polyDesign`, `calcRelB1_unmasked`,
  `safeMean`
- the `phantomB1Method` switch (`extrap` / `measure` / `value`) and its relB1
  estimation
- the V_ref threshold fallback (`V_ref_seg_mL`) and the V_M voxel-count route
  with its `slabCov` B1-slab-to-ISIS-slab scaling
- output fields `nVoxRef`, `nVoxTissue`, `voxVol_mm3`, `b1Slab_mm`, `slabCov`,
  `slicesUsed`, `V_ref_seg_mL`, `brainSource`
- the now-dead `outDir` / `scanID` arguments (they only fed the QC figure), so
  the signature is `extRefQuant31P(fitNames, fitAreas, metProps, TR_ms, TE_ms,
  extRef, internalAbsQ)`

Deleted from both drivers: `findB1Dir` / `findB1DirCalf` and the extRefDefaults
fields `b1Dir`, `b1Source`, `dropEndSlices`, `brainThr`, `phantThr`,
`phantMinVox`, `phantomMask`, `phantomB1Method`, `b1ExtrapOrder`,
`phantomB1Value`.

`utils/readB1.m` had no callers left and has been **deleted** from this repo
(`git rm`, so it is recoverable). Copies elsewhere on the machine are
untouched and remain the canonical ones.

Verified behavior-preserving on 859875_S01: identical to the pre-removal run
(V_ref 8.000 mL, NADplus ext 0.8938 mM, relB1 0.8339 / 0.5720, excModel
`rect 200 us at -4.7 ppm, flip 57 deg`). `checkcode` clean on all five files.

Two self-inflicted breakages caught by running it, both from the bulk deletion:
the `relB1_M` / `relB1_ref` / `b1method` initializers were carried off with the
map estimation block (undefined variables whenever the calibration failed), and
`getfielddef` was swallowed by the adjacent `safeMean` deletion. Restored.

### 20d. Oct 2026: relB1 is rescaled to the acquisition's reference voltage

§20a said relB1 is the one thing that crosses from the calibration to the
acquisition. It crosses, but not unchanged: the B1CALIB fit gives relB1 relative
to the calibration's reference voltage (`b1.TxVolt`), and the acquisition runs at
its own. The achieved flip in the acquisition is prescribed flip x relB1_cal x
(V_acq / V_cal), so both tissue and phantom relB1 are multiplied by that ratio
before the excitation factor is computed. The receive factor relB1_ref/relB1_M
is a ratio and does not change.

- drivers: `extRef.acqRefVolt` = `sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude`
  of the acquisition.
- `extRefQuant31P`: `relB1 = relB1_cal * acqRefVolt / b1.TxVolt`; the calibration
  values and the ratio are kept as `relB1_M_calib`, `relB1_ref_calib`,
  `relB1VoltScale`. A missing voltage on either side stops the quant.

859875_S01 brain (calibration 300 V, acquisition 331 V, scale 1.103): relB1
0.834/0.572 -> 0.920/0.631, eta_ref/eta_PCr 0.6410 -> 0.6512, every external
concentration +1.6% (NAD+ 1.404 -> 1.427 mM, ext/PCr-ref 3.32 -> 3.38). Small
here because the eta ratio is nearly flat in relB1 at 57 deg; it would matter
more at larger flips or larger voltage mismatches.

The same session originally had TWO B1CALIB flip series in E100 (0006-0010 at
-4.7 ppm, 0011-0015 at 0 ppm), which `isis_svs_refVolt` correctly refuses to
pool. The -4.7 ppm series was deleted from the data folder (user's call; the
phantom is only visible at 0 ppm). Sessions are expected to carry one series,
so the auto-locator was left as is.

## 21. Sep 2026: coarse-HSVD ridge and lambdaWidth re-selection

**Outcome: no change.** The ridge in the coarse HSVD pole solve stays, and
`lambdaWidth = 0.005` stays. Removing the ridge, alone or with `lambdaWidth`
re-selected across 0–0.3, is worse than the current config on all three cohorts
taken together, and NADH4 degrades at every setting.

### What the ridge is

`utils/hsvd_fit.m` finds the poles from the shifted signal subspace with a
Tikhonov-regularized solve instead of plain least squares:

```
A = (U1'U1 + λI)⁻¹ · U1'U2        λ = (0.05·s_K)² ≈ 0.0024 (heuristic, lines 61-68)
```

`U1` is nearly orthonormal, so `A ≈ A_LS/(1+λ)`: every pole shrinks by the same
factor, frequencies are untouched, and every component linewidth gains
**ln(1+λ)·bw/π Hz** (amplitudes compensate upward). Measured on real data by
rerunning with `λ = 0`: coarse widths move by exactly **−6.0 Hz at 8000 Hz**
(U01 from Apr 2026) and **−3.0 Hz at 4000 Hz** (Pilot-Sinc; U01 S092 March).

This is the only algorithmic difference between our coarse decomposition and
SpecTickle's `funcs/me_hsvd.m` (plain `pinv`). On a synthetic FID,
`hsvd_fit(..., 0)` reproduces `me_hsvd` to 0.02 Hz / 0.1°. Everything after the
decomposition differs regardless — SpecTickle quantifies from the HSVD component
amplitudes; we only use them for seeds, the width prior and the phase.

### What the HSVD width feeds

The HSVD linewidth for each peak becomes `width0_all`, which is used twice
(`dat2mat_svsNAD.m` ~1047):

1. the fine VARPRO fit's **initial width**, and
2. the **width prior target** (`baseR.priorWidth`).

As a seed alone it is inert: Trp, which is exempt from the prior, saw its seed
move 6 Hz and its fitted width move 0.04 Hz. As the prior target it moves
concentrations. Subtracting the ridge offset from the target is therefore
equivalent to the `λ = 0` arm below (same width shift, identical component
assignment, applied phase within 1°).

### A/B: ridge vs `λ = 0`, all other settings current

U01 (10 sessions, 2 pairs) and Pilot-Sinc (8, 4 pairs). The ridge arm was
bit-identical to the repo on all 18 sessions.

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 median conc change | −7.5% | −12.7% | −13.3% | −0.3% |
| U01 group CV% | 15.0 → 15.1 | 15.1 → 16.7 | 29.4 → 44.8 | 46.8 → 46.5 |
| U01 TRT% | 13.6 → 13.7 | 11.4 → 11.3 | 19.8 → 26.4 | 25.5 → 25.6 |
| Pilot-Sinc median conc change | −7.2% | −8.2% | −7.8% | −0.1% |
| Pilot-Sinc group CV% | 8.6 → 12.1 | 8.0 → 8.9 | 17.9 → 17.7 | 22.0 → 22.4 |
| Pilot-Sinc TRT% | 2.9 → 5.0 | 2.8 → 2.4 | 9.1 → 7.7 | 6.2 → 6.5 |

MAD 0.246 → 0.233 (U01), 0.218 → 0.186 (Pilot-Sinc). Peaks found, component
assignment and residual RMS unchanged.

### `lambdaWidth` re-selected with `λ = 0`

U01 (10) + Pilot-Sinc (8) + DFS (12). "current" = ridge on, 0.005.

| lambdaWidth | PS TRT H2/H6/H4 | U01 TRT H2/H6/H4 | group CV H4 U01/PS/DFS | MAD U01/PS/DFS | DFS H4/H2 | DFS H2/H6/H4 µM |
|---|---|---|---|---|---|---|
| **current** | 2.9 / 2.8 / 9.1 | 13.6 / 11.4 / 19.8 | **29.4 / 17.9 / 29.8** | 0.246 / 0.218 / **0.232** | **0.97** | 400 / 484 / 377 |
| 0 | 9.7 / 7.1 / 28.8 | 13.0 / 14.6 / 29.5 | 34.4 / 32.2 / 32.0 | 0.259 / 0.261 / 0.364 | 0.57 | 330 / 417 / 188 |
| 0.001 | 5.7 / 3.7 / 21.8 | 11.5 / 12.2 / 31.3 | 35.2 / 29.3 / 43.0 | 0.212 / 0.247 / 0.320 | 0.65 | 326 / 403 / 216 |
| 0.002 | **3.3 / 1.9 / 5.0** | 13.7 / 13.4 / 29.2 | 43.8 / 22.9 / 41.8 | 0.236 / 0.195 / 0.290 | 0.72 | 333 / 408 / 250 |
| 0.005 | 5.0 / 2.4 / 7.7 | 13.7 / 11.3 / 26.4 | 44.8 / 17.7 / 46.9 | 0.233 / 0.186 / 0.268 | 0.80 | 337 / 411 / 284 |
| 0.01 | 2.2 / 4.5 / 13.5 | 14.1 / 13.3 / 24.9 | 43.5 / 20.3 / 50.9 | 0.233 / 0.188 / 0.267 | 0.84 | 337 / 407 / 297 |
| 0.02 | 3.1 / 6.1 / 10.1 | 14.2 / 15.2 / 25.8 | 41.0 / 23.2 / 48.3 | 0.226 / 0.197 / 0.250 | 0.90 | 336 / 401 / 311 |
| 0.05 | 3.5 / 4.6 / 11.5 | 14.3 / 16.7 / 27.2 | 40.5 / 22.7 / 51.2 | 0.220 / 0.194 / 0.254 | 0.91 | 337 / 405 / 310 |
| 0.1 | 2.7 / 5.7 / 9.9 | 17.8 / 16.0 / 27.4 | 40.1 / 19.6 / 52.5 | 0.196 / 0.186 / 0.263 | 0.91 | 335 / 403 / 313 |
| 0.3 | 3.5 / 5.9 / 9.4 | 17.4 / 16.7 / 18.1 | 33.8 / 20.0 / 52.7 | 0.173 / 0.195 / 0.263 | 0.92 | 336 / 403 / 313 |

Trp is exempt from the width prior, so every ridge-off arm gives the same Trp:
group CV 46.5 / 22.4 / 19.4%, TRT 25.6 / 6.5% (current 46.8 / 22.0 / 19.9%,
25.5 / 6.2%).

NADH2/H6 group CV moves little across the grid (U01 13.5–16.7 / 16.2–18.0,
Pilot-Sinc 9.2–13.6 / 6.3–12.2, DFS 9.5–13.2 / 16.2–17.6; current 15.0 / 15.1,
8.6 / 8.0, 12.4 / 15.2). A mean z-score composite over TRT H2/H6, group CV H2/H6
(3 cohorts) and MAD (3 cohorts) ranks current first (−0.76), then 0.002 (−0.48),
0.001 (−0.37), 0.005 (−0.27).

The prior becomes a hard constraint quickly: median |wV − prior|/prior is
≤5% by 0.01, ≤3.2% by 0.02 and ≤1.6% by 0.05.

### Why nothing recovers the current config

With the ridge off the width target is ~6 Hz (U01) narrower. The fit follows
it, and the baseline spline rises under all three NAD peaks to take the wings
(visible in nearly every DFS subject; e.g. NADH4 S167 241 → 98, S092 331 → 167,
S169 293 → 139 µM at 0.002). Residuals are marginally *smaller* (DFS median RMS
0.0168 → 0.0131), so fit quality gives no reason to reject that partition — the
data do not decide between them.

Raising `lambdaWidth` does not bring NADH4 back: past 0.05 the widths simply sit
on the narrower target and DFS NADH4 plateaus at 310–313 µM. **It is the target
value, not the prior strength, that matters.** The ridge's extra width has been
acting as a useful regularizer, mostly for NADH4.

0.002 is the best ridge-off value and gives the best Pilot-Sinc TRT in the grid
(the evidence §15 leaned on), but it costs DFS MAD (0.290 vs 0.232), NADH4 in
all three cohorts, and U01 NADH6 TRT. Four pairs in one cohort do not outweigh
that.

As always, none of this says which concentration scale is correct (§10).

### Reproducing

Scratch harness (session scratchpad, not the repo): an isolated copy of the tree
with a harness-only `hsvdCoarse.lambda` option passed to `hsvd_fit` (empty =
heuristic, 0 = plain LS), recorded as `coarse.hsvd.lambda` in the output.
`ridge_run.m` / `ridge_score.m` / `ridge_zoom.m` for the A/B; `lw_run.m` /
`lw_score.m` / `lw_zoom.m` for the sweep. **Start MATLAB with `cd` into the
scratch folder:** launched from `Z:\code\svs_code`, the current folder shadows
the tree copy and the repo code runs instead — the harness asserts on `which`
for exactly this reason.

## 22. Sep 2026: HSVD results saved for comparison with VARPRO

**Outcome:** `dat2mat_svsNAD` now saves a per-peak summary of the coarse HSVD
and an HSVD decomposition of the water reference, so HSVD and VARPRO results can
be compared after processing. **No fit result changed**: concentrations, fit
parameters, `areasConv`, water area/parameters, seeds and applied phase are
bit-identical to the previous code on all 18 brain (U01 + Pilot-Sinc) and 14
calf sessions. Blood (`bru2mat_svsNAD`) was not changed.

### New metabolite fields

`output.met.fit.hsvd` (vectors in `pars.peaks.name` order; NaN where HSVD found
no qualifying component):

| field | meaning | compare with |
|---|---|---|
| `found`, `nComp` | peak found; number of components summed | — |
| `ppm` | \|amp\|-weighted component center | VARPRO `pos` |
| `width`, `widthHz` | \|amp\|-weighted Lorentzian FWHM (ppm, Hz) | VARPRO total width `wV` |
| `ampl` | complex sum of component amplitudes, FID units | — |
| `area` | `N·Δppm·|ampl|/2`, same units as `areasConv` | `areasConv` |
| `phaseDeg` | angle of `ampl`, before coarse phase correction | — |
| `phaseDegCorr` | `phaseDeg` minus the applied coarse phase at the peak | — |
| `phaseDegFit` | `-phaseDegCorr`, the fine fit's `ph` sign convention | VARPRO `ph` |
| `lambda` | `hsvd_fit` ridge actually used (§21) | — |

Also added: `init.found` (HSVD seed vs 30 Hz fallback), `init.widthHz`,
`init.priorWidthHz` (the width prior target the fitter received, NaN where
unused), and in `region(r).coarse.hsvd`: `peakSel` (components × peaks
assignment), `nCompPeak`, `centerPeak`, `widthPeak`, `phasePeakDeg`, `lambda`.

**The `area` scaling.** For a sampled one-sided decaying exponential of
amplitude `a`, the DFT sums to `N·a`; half of that is the Lorentzian line and
half a constant offset from the first point, which the VARPRO spline absorbs.
So the line area is `N·Δppm·|a|/2`, the same factor the `watfit.method = 'fid'`
path uses.

**The phase sign.** VARPRO `ph` runs opposite to the HSVD corrected phase
(brain: corr −0.64 against `phaseDegCorr`, +0.64 against `phaseDegFit`; median
|`phaseDegFit` − `ph`| 7 / 7 / 11° for H2/H6/H4). Same conj(Vc) convention as
§11. Always compare `phaseDegFit` with `ph`.

### What the comparison looks like (current config)

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| brain HSVD area / areasConv (median) | 1.14 | 0.99 | 1.23 | 2.49 |
| brain HSVD − VARPRO center (ppm) | 0.005 | 0.002 | 0.007 | 0.003 |
| brain HSVD width / VARPRO wV (Hz) | 31.4 / 30.0 | 35.7 / 36.1 | 42.5 / 39.2 | 61.4 / 42.4 |

Trp's 2.5x is the broad HSVD Trp component (§13, §15) carrying pedestal. Calf is
noisier: NADH4 found in 71% and Trp in 43% of sessions, and calf HSVD widths run
~10 Hz above VARPRO for NADH2 and NADH6 (NADH4 40.1 vs 42.5 Hz). HSVD metabolite amplitudes carry the
ridge's inflation (§21; +9–18% at 8000 Hz on synthetic lines, larger for
narrower lines).

### Water reference HSVD (`output.wat.fit.hsvd`)

There was no HSVD of the water reference before. It now runs for every
`watfit.method`, on the same spectrum the VARPRO water fit sees, through the
**same code as the metabolite coarse fit**: `hsvdComponents` (the shared
`hsvd_fit` wrapper) and `hsvdAssign`, with the water as a one-row peak table.
Components inside `ppm_range` with linewidth in [`lwMinHz` `lwMaxHz`] are
complex-summed;
fields mirror the metabolite ones (`ampl`, `area`, `ppm`, `widthHz`,
`phaseDeg`, `nComp`) plus the full component list and `assigned`.

Options in `pars.watfit.hsvd`: `enable` (true), `K` (**6**), `npts` (1024),
`lambda` (**[]**, the `hsvd_fit` ridge heuristic), `lwMinHz` (5), `lwMaxHz`
(90), `ppm_range` (`watfit.ppm_range`).

**Why K = 6 with the ridge.** Chosen to match SpecTickle's water model order
(`WaterModelOrderEditField` = 6) and to keep the water decomposition on the same
footing as the metabolite one, which keeps the ridge (§21). An HSVD met/water
ratio therefore carries the same amplitude bias on both sides.

**What that costs.** HSVD describes the strong, not purely Lorentzian water line
with a conjugate **pair** of poles at ±~0.01 ppm whose phases sit near ±83–90°,
so their complex amplitudes very nearly cancel: `|Σa|/Σ|a|` runs 0.01–0.8 and the
small residual sum *is* the line area. That is normal and by itself harmless —
see the failure mode below for when it stops being harmless. Sweep over all 32
sessions, HSVD area / VARPRO `areaTotal` (median [min max]):

| config | brain | calf |
|---|---|---|
| K 60, ridge (metabolite order) | 0.96 [0.05 4.18] | 0.96 [0.67 2.20] |
| K 60, λ = 0 | 0.91 [0.60 1.28] | 0.93 [0.74 1.22] |
| K 20, λ = 0 | 0.92 [0.43 2.19] | 0.94 [0.66 1.28] |
| K 10, λ = 0 | 0.94 [0.66 1.09] | 0.94 [0.86 1.17] |
| K 6, λ = 0 | 0.94 [0.82 1.14] | 0.95 [0.60 1.17] |
| **K 6, ridge (current)** | **0.92 [0.80 1.21]** | **0.88 [0.08 1.24]** |

Brain is well behaved at the current setting (all 18 sessions within
0.80–1.21 of the VARPRO area). **Four calf sessions are not**: S004_04162026
(ratio 0.08), S139_07102026 (0.13), S0000_04172026 (0.45), S0000_04202026
(0.66). The other ten calf sessions sit at 0.87–1.24.

**The failure is `lwMaxHz`, not the cancellation.** The overlays from
`plot_hsvd_varpro_fits` (§23) show the HSVD water model tracking the data in all
four, with a wide pedestal sitting in the *unassigned* trace. In each one HSVD
put a **broad component at the water frequency — 103, 111, 141 and 153 Hz —
which `lwMaxHz = 90` throws out**. That component carries the positive area; the
conjugate pair that remains is left with its near-cancelled integral, so the
assigned area collapses. Recomputing the area from the saved component list
with the cap raised (pure arithmetic on `compAmp`, no reprocessing):

| lwMaxHz | U01 | Pilot-Sinc | DFS | calf |
|---|---|---|---|---|
| 90 (current) | 0.89 [0.80 0.99] | 0.96 [0.85 1.21] | 0.94 [0.89 1.02] | 0.88 [0.08 1.24] |
| 200 | 0.89 [0.80 0.99] | 0.91 [0.83 0.96] | 0.94 [0.89 1.02] | 0.89 [0.87 1.24] |

Every affected session lands near 1.0 (0.08 → 0.95, 0.13 → 0.97, 0.45 → 0.99,
0.66 → 1.00), U01 and DFS are untouched (no >90 Hz component in range at all),
and Pilot-Sinc tightens (its 1.21 outlier → 0.96). `lwMaxHz = 200` in
`pars.watfit.hsvd` is therefore a strict improvement over the default, which was
copied from SpecTickle's *peak-table* water row where a 90 Hz cap is meant to
keep broad baseline out of a metabolite. For the water reference the broad
component **is** the water. Left at 90 pending a decision — it is a `pars`
change, not a code change.

`K = 10` with `lambda = 0` also removes the instability (calf 0.86–1.17) but
changes the decomposition instead of the assignment.

The VARPRO water area itself agrees with the model-free \|complex integral\|
over the window to a median 1.02 (0.87–1.17), so VARPRO remains the stable water
estimate; HSVD water reads ~8% lower on brain and ~12% lower on calf.

---

## 23. Sep 2026: HSVD vs VARPRO across every data set (`compare_hsvd_varpro`)

**Outcome:** with §22's fields in place, a full HSVD-only quantitation can be
compared with the production VARPRO one on all 44 sessions. The two agree on the
peaks that are genuinely Lorentzian and diverge everywhere else. Brain
concentration, HSVD vs VARPRO (median per cohort):

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 (n=10) | +31% | +11% | +37% | +193% |
| Pilot-Sinc (n=8) | +9% | +1% | +23% | +194% |
| DFS (n=12) | +14% | +5% | +44% | +209% |

Calf HSVD is not usable at all. Nothing in the pipeline changed.

### The script

`compare_hsvd_varpro.m` (repo root). Runs every session through its production
driver, then reads `output.met.fit.hsvd` / `output.wat.fit.hsvd` back out and
compares them with the VARPRO fit in the same output.

```matlab
compare_hsvd_varpro                              % process what is missing, then report
compare_hsvd_varpro('cohorts',{'u01','calf'})
compare_hsvd_varpro('run',false)                 % report from existing outputs only
```

Writes to `_hsvd_varpro\`: the processed `.mat` per session, `*_peaks.tsv`
(session × peak), `*_water.tsv`, `*_report.txt`, and one PNG per cohort.
Cohorts: `u01` (10 brain), `pilotsinc` (8), `dfs` (12), `calf` (14).

`plot_hsvd_varpro_fits.m` draws the two fits on top of the data from those same
outputs — per fit region and for the water reference, data / VARPRO (fit +
spline baseline) / HSVD (all components, unassigned dashed), residuals below.
Per-session PNGs in `_hsvd_varpro\fits\<cohort>\`, plus `fits_<cohort>_metabolite.png`
and `fits_<cohort>_water.png` overview grids. The metabolite HSVD traces are the
ones the coarse step already stores (`region.coarse.fitPhased` /
`.baselinePhased`, in the data's phase frame); the water model is rebuilt from
the saved component list with `utils/hsvdSpectra.m` — which is the same function
`dat2mat_svsNAD` uses, moved out of it for this purpose (verified bit-identical
on one brain and one calf session).

What the overlays show for calf: both models track the data closely in 12 of 14
sessions and differ mainly in where they put the **baseline** — the HSVD
baseline (unassigned components) climbs steeply under the NAD cluster where the
VARPRO spline stays flat, which is the ~+46% NADH2/NADH6 area difference. Two
panels are worth knowing about: S092_04232026, where the *VARPRO* fit visibly
misses 9.3–9.5 ppm (NADH6 = 129 µM, the low outlier) and HSVD follows the data,
so its 2.19 area ratio is VARPRO under-fitting; and S139_06052026, where the
HSVD baseline dives to −1e-4 on 1e-4 data and its ratio of 3.68 is meaningless.

**HSVD concentrations are computed in the script, not in `dat2mat_svsNAD`.**
`absoluteQuant_svs` is linear in (metabolite area / water area) and the
segmentation, T1, T2 and timing are shared, so

```
conc_HSVD = conc_VARPRO * (areaHSVD_met/areaVARPRO_met) / (areaHSVD_wat/areaVARPRO_wat)
```

Both sides of the ratio are HSVD — the metabolite area **and** the water area.

Pilot-Sinc is driven through a generated copy of `processU01_Brain_1H` with the
.dat pair pinned (`*2X9_7ppm*` / `*WATERREF_2x4_7ppm*`): the production picker
reads the water reference as water-suppressed, because its name contains "ppm",
and then takes the highest MID, which is the 3X10 voxel. Fit settings are the
driver's own.

### Brain: the area ratio tracks the Gaussian fraction

HSVD area / `areasConv`, median per cohort, with the VARPRO Voigt η (`wG/wT`)
that peak was fitted with:

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 (n=10) | 1.14 | 1.00 | 1.21 | 2.48 |
| Pilot-Sinc (n=8) | 1.02 | 0.96 | 1.23 | 3.11 |
| DFS (n=12) | 1.06 | 0.99 | 1.31 | 2.82 |
| median η (VARPRO) | 0.00–0.14 | 0.00–0.01 | 0.01–0.02 | 0.87–0.99 |
| median width HSVD / VARPRO | 1.04–1.11 | 0.99–1.02 | 1.07 | 1.37–1.43 |

NADH6 is essentially pure Lorentzian (η ≈ 0) and the two areas agree to 1–4% in
all three cohorts — an independent check that §22's `N·Δppm·|A|/2` scaling is
right. Trp is fitted near-Gaussian (η 0.87–0.99) and HSVD is forced to describe
it with a Lorentzian ~1.4x wider: a Lorentzian carries ~1.5x the area of a
Gaussian of the same height and FWHM, so ~2x is expected from lineshape alone
and the rest is the broad Trp pedestal (§13, §15). NADH4's +21–31% is the
window-clipping peak (§15) where HSVD and the spline baseline divide the wing
differently.

The concentration differences above are these area ratios divided by the water
area ratio (0.89–0.96, §22), which is why every peak shifts up by a further
~5–13%. Centers agree to 0.002–0.014 ppm; median |`phaseDegFit` − VARPRO `ph`|
is 5–22°.

### Brain: precision and self-consistency

Group CV% and test–retest CV% of concentration, VARPRO → HSVD:

| | NADH2 | NADH6 | NADH4 | Trp |
|---|---|---|---|---|
| U01 group CV | 15.0 → 15.8 | 15.1 → 18.8 | 29.4 → 30.1 | 46.8 → 38.8 |
| U01 test–retest (2 pairs) | 13.6 → 10.1 | 11.4 → 9.2 | 19.8 → 21.4 | 25.5 → 64.0 |
| Pilot-Sinc group CV | 8.6 → 26.8 | 8.0 → 28.5 | 17.9 → 30.4 | 22.0 → 18.2 |
| Pilot-Sinc test–retest (4 pairs) | 2.9 → 10.4 | 2.8 → 15.9 | 9.1 → 17.6 | 6.2 → 17.1 |
| DFS group CV | 12.4 → 14.2 | 15.2 → 20.7 | 29.8 → 68.9 | 19.9 → 15.2 |

Pilot-Sinc is the clearest case: VARPRO test–retest 2.8–9.1% becomes 10–18%
under HSVD, on the same spectra. H6/H2 and H4/H2 (MAD from 1.0, VARPRO → HSVD):
U01 0.246 → 0.153, Pilot-Sinc 0.218 → 0.262, DFS 0.232 → 0.378. HSVD improves
U01's ratio balance (H6/H2 1.25 → 1.08) and worsens DFS's (H4/H2 0.97 → 1.46
with SD 0.96) — no consistent direction, and this is self-consistency, not
accuracy (there is no gold standard here).

### Calf: HSVD does not survive

Per-peak availability is already poor — NADH4 found in 71%, NR 79%, Trp 43% of
the 14 sessions — and group CV% of concentration goes 26–52 (VARPRO) to 84–186
(HSVD). Two sessions blow up through the water denominator (S004_04162026 ratio
0.08, S139_07102026 0.13; §22): NADH2 reads 25 870 and 10 536 µM against
VARPRO's 1194 and 909. Median area ratios are 1.46 (NADH2), 1.10 (NADH6), 0.47
(NADH4) with group CVs of the ratio of 48–177%. Calf lines are broader, less
resolved and fitted with larger η (0.07–0.44 on the NAD peaks, 0.93 on
Peak10_3), which is exactly the regime where the Lorentzian-only model fails.

### Calf NR row widened to [9.60 9.80] (2026-09-15)

The overlays raised a second assignment problem, this one in the peak table.
`processU01_Calf_1H` had NR at `[9.60 9.75]`, but across the 14 calf sessions the
HSVD component near 9.7 sits at **9.717–9.765 ppm, median 9.746** — on the top
edge. Three sessions fell outside it (9.761, 9.761, 9.765) and HSVD called the
NR peak baseline; the VARPRO center, bounded by the same row, railed at 9.7500
in two sessions and sat within 3 mppm of it in three more. NADH2/NADH6 come out
at 9.35/9.17 in those scans, so the axis is right and the line really is ~0.045
above the nominal 9.70.

Widening the row exposed a second, bigger bound problem. `processU01_Calf_1H` set
no `widthBoundsHz`, so the **default [10 60] Hz** applied to all six peaks while
the calf coarse HSVD is allowed 90 Hz: in 15 of 84 peak-sessions the HSVD seed
was already above the fit's own upper bound, and Trp sat pinned on the cap in 3–4
sessions. Both settings are now changed in the driver:

```matlab
pars.peaks.range(4,:)              = [9.6 9.8];   % NR
pars.fit.baselineOpt.widthBoundsHz = [10 90];
```

Four configurations, all 14 calf sessions (relative residual RMS = ‖resid‖/‖data‖
over the fit window):

| | A 9.75/60 (old) | B 9.80/60 | **C 9.80/90 (current)** | D 9.75/90 |
|---|---|---|---|---|
| HSVD finds NR | 11/14 | 14/14 | **14/14** | 11/14 |
| NR centers railed | 2 | 0 | **0** | 4 |
| widths railed (all peaks) | 6 | 7 | **0** | 1 |
| seeds above the cap | 15 | 15 | **0** | 0 |
| median residual | 0.0136 | 0.0140 | **0.0134** | 0.0132 |
| worst residual | 0.0933 | 0.0935 | **0.0929** | 0.0925 |

The width cap is the clear win: **S092_04232026 was simply broken at 60 Hz** —
residual 0.0880, the fit visibly missing 9.3–9.5 ppm, NADH6 reading 129 µM — and
at 90 Hz it drops to 0.0143 with NADH6 at ~500 µM. Group CV improves for NADH6
(30.7 → 22.0) and NR (49.7 → 44.4); NADH4 and Trp grow instead (medians 984 →
1383 and 105 → 161 µM, widths 42.5 → 56.5 and 40.6 → 49.1 Hz), which is the
peak-versus-baseline trade again: give a peak room and it takes area off the
spline.

C vs D (the NR row is the only difference) is a wash on fit quality — C better in
8 sessions, D in 6, medians 0.0134 vs 0.0132 — so the row was decided on the
artifacts: D rails 4 NR centers on the 9.75 bound, C rails none and lets the
centers land where the data put them (9.72–9.77). **C is what is in the driver.**

The one cost is **S139_03182026**, which the NR row (not the width cap) flips into
a worse optimum: A 0.0114 → B 0.0167 → C 0.0157, with NR 255 → 553 µM and Trp
83 → 357 µM as the downfield trio widens. Running that scan at NR 9.75 with the
90 Hz cap (D) returns it to 0.0119, NR 274, Trp 123. Treat its NR/Trp/Peak10_3
numbers as unreliable, and note what it says generally: the single-region calf
layout couples all six peaks through one baseline and one coarse phase, so a seed
change anywhere can move everything. Whether a resonance at 9.745 is NR/NMN is a
separate question none of this settles.

### Two-region calf layout: tried and rejected (2026-09-16)

The obvious structural answer to that coupling is the `perRegion` layout still
sitting commented out in `processU01_Calf_1H` — NAD cluster in one window, the
NR/Trp/10.3 group in another, each with its own coarse fit, phase and baseline.
Both variants were run on all 14 sessions (shadow drivers, nothing in the repo
touched), keeping NR `[9.6 9.8]` and the 90 Hz cap:

* **E1** `[8.95 9.5]` + `[9.5 10.7]` — the split alone, same total span as now
* **E2** `[8.70 9.5]` + `[9.5 10.4]` — the commented block exactly

| | C 1 region (current) | E1 2 regions | E2 2 regions (as commented) |
|---|---|---|---|
| median residual | **0.0134** | 0.0167 | 0.0258 |
| worst residual | 0.0929 | 0.0315 | 0.0322 |
| better / worse than C | — | 2 / 12 | 1 / 13 |
| NADH2 median, group CV | 997 µM, 26% | 464 µM, **68%** | 1467 µM, 50% |
| NADH4 median, group CV | 1383 µM, 56% | 815 µM, 43% | **164 µM, 105%** |
| Trp group CV | 53% | 81% | 79% |
| NR group CV | 44% | **27%** | 41% |
| widths railed | 0 | 2 | 8 |
| median η, NADH4 / NR | 0.10 / 0.08 | **1.00 / 1.00** | 0.87 / 1.00 |

The overlays say why. In a 0.55 ppm window the spline baseline **arches inside
the region** — in E1's NAD window it rises to ~60% of the peak height mid-window
instead of being the gentle envelope it is across the full range — so the peaks
lose that area (NADH2/NADH6 medians roughly halve) and η rails at 1.0, the Voigt
going pure Gaussian to follow what is left. The downfield region under-fits
instead: the VARPRO trace sits visibly below both NR and the 10.3 peak with
structured residual. E2, which gives the NAD window the extra 0.25 ppm that E1
clips off NADH4, is worse still on every global measure.

NR itself is the one thing the split improves (group CV 44% → 27%) — expected,
since it is no longer competing with the NAD cluster — but not at the price of
NADH2 at 68% CV. **Keep the single region.** If this is ever revisited, the lever
is per-region baseline settings (`regions(r).baselineOpt.knotSpacing`, as the
brain driver sets 0.43 / 2), not the split by itself; both arms here left
`knotSpacing` at the global 2 ppm. Neither arm rescues S139_03182026 either
(E1 0.0177, E2 0.0242 against C's 0.0157), so the coupling was not what flipped
that session.

### S139_06052026: the coarse phase points the wrong way (2026-09-16)

This calf scan has the worst fit in every configuration tried (relative residual
0.093 against a cohort median of 0.013, unchanged by the NR row, the width cap or
the region split). The cause is the applied coarse phase, and it is worth
recording because nothing in the pipeline flags it.

`shiftSpectrumPhase(spec, [pc0 pc1], ppm, pivot)` applies **exp(−i·pc0)**, so the
stored `phaseCorr.pc0 = +27.7°` is a −27.7° rotation of the spectrum. Sweeping
the rotation on this scan's own data: the fit window stays non-negative from
about −10° through +110°, and the real-part integral peaks at +55°. Every other
calf session lands in that range (applied rotations −8° to +79°); this one lands
at −27.7°, past the edge. The consequences:

* **63% of the phased fit window is below zero** — the only calf session with any
  negative points at all; the unphased spectrum is smooth and positive across the
  whole window.
* The spline baseline (`knotSpacing` 2 ppm, λ = 1) cannot follow a trough, so the
  fit misses 9.5–10.4 ppm entirely, equally badly on both halves of the window.
* The peaks absorb the mismatch: **every** peak goes pure Gaussian (η ≈ 1 against
  0.0–0.4 elsewhere) and three phases pin near the ±60° bounds (+59, −54, −56°,
  a 113° spread where a good scan spans ~45°).
* Concentrations come out lowest in the cohort (NADH2 384 µM vs a median of 997)
  and the HSVD/VARPRO area ratio is 3.6 — HSVD builds the trough from unassigned
  components and so keeps far more area in its peaks.

Re-running the scan with the VARPRO coarse estimator (`hsvdCoarse.enable =
false`) gives a different phase and drops the residual to 0.040 with
concentrations nearer the cohort (NADH2 470, NADH6 916, NADH4 243).
`coarseFirstOrder = true` reaches 0.027 but only by fitting **pc1 = −1695°** —
the linear phase fit is degenerate on this spectrum, not a usable setting. The
coarse HSVD phase estimate is unreliable here because the NAD peaks are small
against a large sloping background and its components straddle both.

Treat this scan's 1H numbers as invalid.

#### `pars.fit.coarseQC` — the guard, added 2026-09-16

`dat2mat_svsNAD` can now check this automatically, **off by default**:

| field | default | meaning |
|---|---|---|
| `enable` | `false` | run the check |
| `maxNegFrac` | `0.25` | fail when more than this fraction of the phased fit window is negative |

On failure the coarse step is redone for that region with the VARPRO estimator
(`hsvdCoarse.enable = false` internally) and the fit proceeds from its seeds and
phase. Both flags are saved in `output.met.fit.coarseQC`:

| field | meaning |
|---|---|
| `enabled` | whether the check ran |
| `failed` | whether any region failed it |
| `negFrac` | negative fraction per region, after any fallback |
| `negFracRaw` | the same for the unphased spectrum, as a reference |
| `maxNegFrac` | the threshold used |
| `fallback` | `'none'` or `'varpro'`, per region |

`negFrac` is recorded even when the check is disabled, so an existing pipeline
gets the diagnostic without any change in behavior. The threshold is safe by a
wide margin: across all 44 sessions, 43 have **exactly 0%** of the phased window
negative and S139_06052026 has 63%.

In union mode a fallback loses the component decomposition, so the region
baselines carry no HSVD prior and `output.met.fit.hsvd` stays NaN — the fit still
runs, and the per-region `coarse` struct simply holds the VARPRO coarse result.

Verified: with the check off, brain (union) and calf (perRegion) outputs are
bit-identical to the previous code; with it on, S139_06052026 falls back
(negative fraction 0.63 → 0.13, relative residual 0.093 → 0.040) while a
well-behaved session is untouched; and the union fallback path was exercised by
forcing the threshold negative on a brain session.

**Enabled for calf 2026-09-17** (`processU01_Calf_1H`: `pars.fit.coarseQC.enable
= true`) and all 14 sessions reprocessed. It fires on S139_06052026 and nothing
else — the other 13 are bit-identical with the check on or off. That scan's
concentrations rise 10–88% (NADH2 384 → 470, NADH6 831 → 916, NADH4 158 → 243,
NR 65 → 103, Trp 62 → 83, Peak10_3 61 → 115 µM) and its peak phases come back
inside ±50°. It remains the weakest scan in the set. The brain drivers do not
enable it; no brain session has ever produced a negative phased window.

One knock-on: a session that falls back has **no HSVD summary** —
`output.met.fit.hsvd` is all NaN and `region(r).coarse` carries no
`baselinePhased`/`fitPhased`, since `coarseVarpro` produces no components.
`plot_hsvd_varpro_fits` now draws the VARPRO trace alone for such a region and
labels the panel "coarse-phase QC fallback"; any other consumer of
`region(r).coarse` must tolerate those fields being absent.

#### Coarse phase estimate reworked + phase cap (2026-10-09)

Prompted by calf 859875_S01_10062026. Its raw window is 0% negative, and the fine
fit at 0° applied phase gives a relative residual of 0.019, but production gave 0.079.
Both coarse estimators had rotated it in the wrong direction:

* **HSVD:** NADH4 was assigned an inverted component (−150°) at 8.802 ppm, below
  the 8.95 window edge. The sort-by-ppm `unwrap` then moved every later peak by
  360°, giving pc0 = −333.4° (= +26.6°). The window came out 64% negative, so QC failed.
* **VARPRO fallback:** `pc0 = mean(ph0)` was unweighted, and NADH4 (−66.6°,
  centered at 8.832) dominated it. The result was a +18.3° rotation, 25.4% negative,
  which is over the limit, but the fallback was never re-checked.

Changes in `dat2mat_svsNAD`:

1. **QC picks among three options** (`coarseQCPick`). On failure it scores the HSVD
   phase, the VARPRO phase and no correction by negative fraction and keeps the
   least negative. Ties go to HSVD, then VARPRO. `fallback` is now `'none'`,
   `'varpro'` or `'zeroPhase'`. With no correction, the VARPRO seeds are kept.
   New fields: `coarseQC.negFracHsvd` and `.negFracVarpro`.
2. **Edge peaks don't vote.** A peak whose estimated center lies outside the fit
   window still seeds the fine fit but is left out of the phase estimate. This
   applies to both estimators (`coarseFit.phaseUse`).
3. **Weighted circular mean** (`peakLinearPhase`). Both estimators now weight
   peaks by signal² (HSVD |A|², VARPRO `areasConv`²) and average angles, so the
   ±180° wrap no longer matters. For first order, each phase is first put on the
   branch nearest that mean, then a weighted line is fitted. `phResidDeg` is
   wrapped to ±180°.
4. **`pars.fit.coarsePhaseMaxDeg`** (default `Inf`; set to **25** in
   `processU01_Brain_1H` and `processU01_Calf_1H`). `capLinearPhase` clamps the
   applied phase line to ±max at both fit-window edges, for both estimators.
   Output: `phaseCorr.maxDeg`, `phaseCorr.capped`, `region(r).capped`, and the
   uncapped estimate in `coarse.phaseEst`.

Fixes 2–3 move the brain HSVD phase by at most 2.3° (U01, Pilot-Sinc, DFS, 859875;
replayed offline from `spec_full`). In calf they matter on 10062026 (+26.6 → +16.1)
and S092_04232026 (303.3 → −57.1; the old value was an unwrap artifact).

**Choosing the cap.** For 32 brain and 16 calf sessions, the fine fit was replayed
from `spec_full` at applied zero-order phases −45…+45° in 5° steps (scratch
`sweep_phase.m`). What the sweep shows:

* **Positive phase errors break the fit; negative ones barely matter.** Residuals
  are almost flat from −45° up to a few degrees below 0, then climb steeply past
  about +5 to +10° in every cohort (U01 brain S092_03032026: 0.022 at −20, 0.026 at
  0, 0.031 at +10, 0.045 at +20).
* **The data are not centered on 0.** The lowest-residual phase is negative in all
  16 calf sessions (−5 to −35°, median ≈ −17.5°). In brain it ranges from −35 to
  +15°, with 859875 brain at −20 and −35.
* **Effect of a symmetric cap** (each session's own estimator plus QC, residual
  divided by its sweep minimum, all 48 sessions): cap 25 gives a median of 1.07 and
  3 sessions over 1.5 (the fewest of any value). 30, 45 and Inf are the same or
  slightly worse. 5–20 are worse, and on S139_06052026 caps of 10–20 are actively
  harmful: clamping HSVD's +28.5° to +20 drops the window below the QC limit, so
  the fallback never runs (ratio 3.67, against 1.14 with a 25° cap).
* Cap 25 only changes four sessions whose estimates are ≥ 25° in size: DFS S167
  (−44 → −25, improved), calf S092_04232026 (−57, improved), S0000_04172026 (−47,
  neutral) and S004_04162026 (−82, worse: 1.71 vs 1.22).

**10062026 is not fixed by any symmetric cap.** Its new HSVD phase is +16.1°,
leaving the window only 10.4% negative, so QC passes and the fit runs at +16.1°.
The residual is 0.036, against 0.079 before, 0.019 at 0° and 0.013 at −25°.
Possible next steps, not done yet: tighten calf `maxNegFrac` (only 10062026 and
S139_06052026 have any negative window under the new HSVD phase); pick the phase by
fine-fit residual over [−cap, cap]; or apply a fixed calf offset near −15°.

#### Calf: fixed coarse phase of −20° (2026-10-09)

Calf NADH4 is not used in quantitation (user, 2026-10-09), so calf phase
strategies were re-scored on NADH2, NADH6, NR, Trp and Peak10_3 over the 16
swept sessions (scratch `calf_strategies.m`). Residuals are relative to each
session's sweep minimum; repeat pairs are median |log ratio| of area/water over 6
pairs.

| strategy | median residual vs best | NADH2/NADH6 median (IQR) | pairs NADH2 / NADH6 / NR / Trp | NADH2 group CV |
|---|---|---|---|---|
| current (HSVD + QC) | 1.14 | 0.85 (0.72–0.97) | 22 / 30 / 25 / 62 % | 23 % |
| best of HSVD/VARPRO/zero by residual | 1.07 | 0.87 (0.75–1.02) | 19 / 27 / 12 / 11 % | 18 % |
| always zero | 1.30 | 0.77 (0.70–0.85) | 11 / 15 / 6 / 30 % | 19 % |
| fixed −20° | 1.07 | 0.99 (0.86–1.09) | 12 / 22 / 11 / 29 % | 14 % |

NADH2/NADH6 falls steadily with applied phase: 1.27 at −45°, 1.05 at −30°, 0.99 at
−20°, 0.87 at −10° and 0.77 at 0°. The median residual is lowest at −20° (1.07).
Those two all-session measures agree; the repeat pairs are too few to separate −20°
from −25° (NADH6 jumps 8 → 22% between them). A fixed phase also removes the 10–30°
session-to-session jitter of any per-session estimate on the flat negative side of
the residual curve, and that jitter feeds straight into NADH2/NADH6.
NADH2/NADH6 ≈ 1 is self-consistency, not truth.

**`pars.fit.coarsePhaseFixedDeg`** (default empty, meaning not used): when finite,
the coarse step still runs its estimator (HSVD seeds; the estimate is kept in
`coarse.phaseEst`), but the applied zero-order phase is this value, in
`shiftSpectrumPhase`'s convention (applied as exp(−i·φ)). The cap and the
coarse-phase QC do not run (`coarseQC.enabled` is saved as false). Saved in
`phaseCorr.fixedDeg`. **Set to −20 in `processU01_Calf_1H`.** `coarseQC.enable`
and `coarsePhaseMaxDeg` were removed from that driver because they no longer do
anything there.

Data: 14 sessions on `%CustomerSeq%\svssel_20260101` and 2 early sessions
(S139_03122026, S139_03182026) on `svssel_latest`, all TR 500 / TE 13 ms. A new calf
protocol or coil needs this value re-derived.

**Rerun, 2026-10-09:** all 17 calf sessions went through `processU01_Calf_1H` into
their own processed folders: 15 U01 sessions, including S077_07232026 (written to
`processed_incomplete_31Pfailed`), plus both 859875 scans. The previous 1H outputs
are in `prev_1H_20261009\`. Compared with the previous production outputs (Sep 16
for 14 sessions, Oct 8 for 859875, Jul 26 for S077):

* median residual 0.0143 → 0.0126; worst 0.097 → 0.028;
* group CV: NADH2 36 → 14%, NR 48 → 30%, Trp 107 → 52%; NADH6 22 → 24%,
  Peak10_3 34 → 40%;
* NADH2/NADH6 median 0.85 → 0.94 (IQR 0.71–1.16 → 0.85–1.06);
* repeat pairs (6, median |log ratio|): NADH2 8 → 13%, NADH6 24 → 9%, NR 29 → 20%,
  Trp 49 → 42%, Peak10_3 43 → 28%.

**The fine fit is multimodal, and a fixed phase doesn't change that.** Sweeping
the new outputs:

* **859875 10062026:** at exactly −20° it lands in a bad mode (Trp area 7.0 →
  475 µM, Peak10_3 dispersion tail under a broad Trp, NADH2 −20%). At −25° and −15°
  it gives Trp ≈ 1.3–1.4 (≈ 90 µM) with equal or lower residual.
* **S004_04162026** jumps modes between −35° and −30° (residual 0.011 → 0.020). At
  −20° NADH2 is 679 µM, against 1100–1300 elsewhere.
* **S077_07232026** (out of sample, not in the sweep) keeps improving to −45°
  (residual 0.0175, against 0.028 at −20° and 0.099 at 0°). Its NADH6 is low (286 µM)
  at every phase: a sharp notch near 9.1 ppm bends NADH6, so this is not a phase effect.

Next lever: a multi-start fine fit (best of N starts), as for 31P.

#### Calf multi-start, phase range, kernel and shared NAD lineshape (2026-10-09)

All tests used the 17 calf sessions at the fixed −20° phase, replayed from the
saved outputs. The 859875 scans got 49 starts and the others 33. Each set was one
default start, 75% narrow starts (center ±30% of the peak window, width ×0.67–1.5,
phase N(0, X/4)) and 25% wide starts (uniform within the bounds, random Gaussian
fraction), with the same random draws for every configuration. Scratch:
`multistart_phr*.m`, `multistart_tie.m`, `compare_configs.m`.

* **One start is not enough.** The default start is within 1% of the best cost in
  only 3/17 sessions (median gain about 8%; S004_04162026 0.0204 → 0.0119,
  859875 10062026 0.0175 → 0.0121).
* **Phase range.** ±15° is too tight: median residual +7%, up to +44%. ±30° and ±60°
  fit about equally.
* **Kernel off** costs about 1.3–2.5% in residual and makes the search much easier:
  10 random starts reach the best fit in a median 96–99% of sessions, against
  51–66% with the kernel on.
* **Unshared, 859875 NADH6 is not identifiable.** Within 2% of the best cost,
  10062026 NADH6 spans 1340–1683 (kernel on) or 1055–1590 µM (kernel off). NADH4 is
  centered on the window edge (all fitted centers 8.931–8.955), so its width and
  amplitude trade freely against NADH6 and the baseline. A tighter NADH4 center row
  would not bind.

New `baselineOpt.widthGroup` / `.etaGroup` in `curvefitAuto_varproBaseline`: peaks
with the same positive group number share one total width or one Gaussian
fraction. The solver sees one parameter per group (`applyParTies`/`expandPars`).
When neither is set, the solver call is unchanged and the result bit-identical
(checked on 859875 09012026). `dat2mat_svsNAD` takes both in `pars.peaks.name` order
and subsets them per region.

All rows below are ±60°. "Shared width" is the total width only; "shared
width+shape" also shares the Gaussian fraction. Repeat pairs and NADH6 spread are in
percent.

| config | median residual vs best unshared | H2/H6 median (IQR) | repeat pairs H2 / H6 / NR / Trp / P10 | NADH6 spread within 2% of best (median, max) |
|---|---|---|---|---|
| kernel on, unshared | 1.000 | 1.00 (0.87–1.05) | 9 / 12 / 23 / 45 / 15 | 6, 23 |
| kernel on, shared width+shape | 1.066 | 1.06 (0.97–1.12) | 9 / 8 / 31 / 20 / 15 | 4, 8 |
| kernel off, unshared | 1.013 | 1.07 (0.89–1.12) | 9 / 17 / 22 / 20 / 4 | 13, 52 |
| kernel off, shared width | 1.101 | 1.05 (0.98–1.13) | 15 / 13 / 27 / 15 / 4 | 5, 26 |
| kernel off, shared width+shape | 1.098 | 1.08 (1.01–1.13) | 8 / 9 / 20 / 17 / 9 | 4, 29 |

859875 NADH6 with kernel off and shared width+shape (best fit, with the range
among fits within 2% of it): 1082 (1062–1106) for 09012026 and 1052 (1008–1102) for
10062026. NADH2 is 1268 and 1322. Sharing costs 7–10% in residual, but it is
physically motivated: the same molecule in the same voxel gives the same lineshape.
Width-only sharing leaves the two scans' NADH6 11–12% apart.

**Production calf configuration (2026-10-09):**

* `coarsePhaseFixedDeg = -20`
* `ph_range = [-60 60]`
* `lineShapeOpt.enable = false`
* `baselineOpt.widthGroup = etaGroup = [1 1 1 0 0 0]` (NADH2/NADH6/NADH4)
* `nStarts = 32`

New `pars.fit.nStarts` (default 1, unchanged behavior): `nadFineFit` refits from
that many starts (`fineStarts`, fixed seed, so it is reproducible; start 1 is the
single-start fit) and keeps the lowest `output.cost`, which is the solver's own
objective (new `curvefitAuto_varproBaseline` output). Diagnostics are in
`region(r).fit.multiStart`. The starts run on an open parallel pool if one exists,
serially otherwise, with identical results either way. A bug surfaced on the first
run: the water fit inherits `pars.fit.baselineOpt` and so got the 6-peak groups; the
sharing fields are now stripped when it inherits.

Rerun of all 17 calf sessions; the previous outputs are in `prev_1H_20261009_fix20\`.
Compared with the fixed −20° single-start run:

* repeat pairs NADH2 13 → 6%, NADH6 9 → 8%;
* NADH6 group CV 24 → 19%; NADH2 14% (unchanged);
* NADH2/NADH6 0.94 (0.85–1.06) → 1.08 (1.02–1.13);
* **859875: 09012026 NADH2 1277 / NADH6 1100, 10062026 1332 / 1077 µM** (4% and 2%
  apart; residuals 0.0141 and 0.0130, at their sweep minima).

Still open:

* **NR / Trp / Peak10_3.** Median repeat pairs improve (Trp 42 → 18%, Peak10_3
  28 → 8%), but there are outliers. In S139_03122026 (Trp 561, Peak10_3 737 µM)
  Trp and Peak10_3 carry large, opposite dispersion lobes that cancel, a phase
  trade-off between neighbors allowed by the ±60° per-peak phases; S167 and
  S004_04162026 look similar. Candidate fix: a shared phase, or narrower per-peak
  phase bounds, for the downfield group.
* **S077_07232026** gets worse (residual 0.028 → 0.040; NADH2 793, NADH6 569 µM): a
  shared lineshape cannot follow its sharp notch near 9.05 ppm. Looks like a data
  problem.
* **Some sessions are still under-searched.** S0000_04202026 found its best at
  start 32 of 32, and 859875 09012026 had only one start within 1%, although its
  value matches the 49-start test (1277 vs 1268). Consider 48–64 starts.

#### Downfield group: shared phase vs narrower phase limits (2026-10-09)

New `baselineOpt.phaseGroup` (shared phase, modes 3/4/5), same mechanism as
`widthGroup`. Tested on the production calf fit over 17 sessions (49/33 starts;
scratch `multistart_df.m`, `analyze_df.m`), changing only NR/Trp/Peak10_3:

| variant | median residual vs A (max) | pairs NR / Trp / P10 | CV NR / Trp / P10 | near-tie spread (within 2%) NR / Trp / P10 | phase at bound |
|---|---|---|---|---|---|
| A independent ±60 (production) | 1.000 | 20 / 17 / 9 % | 37 / 43 / 37 % | 11 / 36 / 22 % | 3/17 |
| B shared ±60 | 1.028 (1.128) | 20 / 34 / 23 % | 34 / 56 / 29 % | 13 / 43 / 14 % | 0/17 |
| C30 independent ±30 | 1.008 (1.095) | 23 / 21 / 15 % | 37 / 48 / 43 % | 12 / 50 / 24 % | 10/17 |
| C15 independent ±15 | 1.029 (1.187) | 21 / 58 / 16 % | 34 / 58 / 38 % | 9 / 37 / 16 % | 14/17 |

None beats A. The shared phase also disturbs the NAD peaks through the common
baseline: NADH2/NADH6 repeat pairs go 8/9 → 11/14%, and 859875 10062026 NADH6
1052 → 963 µM. The narrower limits rail (10/17 and 14/17 sessions).

**The outliers are not phase-driven.** Rerun as A with a different start set,
S139_03122026 gives Trp 354 / Peak10_3 312 / NR 753 µM at residual 0.0116.
Production found Trp 561 / Peak10_3 737 / NR 665 at 0.0115. These are near-equal
solutions, and Trp's spread among near-ties is 36–50% in every variant. Trp and
Peak10_3 (small and broad, beside each other on a sloping baseline) are not
well determined in this window. **Production stays A**; `phaseGroup` is
available but unused.

### Calf fit range down to 8.8 ppm: tried and rejected (2026-09-16)

The 2x2 over the NR row and the width cap left NADH4 unexplained: the cap moves
its median 40% (984 → 1383 µM) with the CV stuck at 52–56%. The suspicion was
window clipping — NADH4's center row is `[8.80 9.05]` while the calf fit range
starts at **8.95**, so a third of its own center window holds no fitted data.
Extending the region to `[8.8 10.7]` (NR `[9.6 9.8]`, cap 90, `block_align`
left alone) does not fix it — it makes everything worse:

| | C range 8.95 (current) | F range 8.80 |
|---|---|---|
| NADH4 median, CV | 1383 µM, 56% | **1959 µM, 64%** |
| NADH4 median width | 56.5 Hz | **89.5 Hz** (railed at the 90 Hz cap in 7/14) |
| NADH4 centers on the 8.80 bound | 0 | **4/14** |
| NADH6 median, CV | 1016 µM, 22% | 655 µM, **49%** |
| NR median, CV | 315 µM, 44% | 219 µM, 59% |
| median η, NADH4 / NR | 0.10 / 0.08 | 0.74 / 0.76 |
| residual on the common window [8.95 10.7] | **0.0134** | 0.0184 (C better in 12/14) |

Opening the window lets NADH4 run into the corner — center on its lower bound,
width on the cap — and swallow the steep background rising below 8.95, which the
model has no component for. Everything else then re-partitions around it
(NADH6 loses a third of its area, η climbs everywhere as peaks go Gaussian to
follow the envelope), and S092_04232026 breaks again (0.0143 → 0.0690).

So the 8.95 edge is load-bearing: it excludes background the model cannot
represent. NADH4 stays ill-determined either way — clipped and edge-pinned at
8.95 (its fitted center sits within 0.01 ppm of the edge in most sessions), or
divergent at 8.80. The next thing to try is not a wider window but a **broad
component** below 9 ppm to carry that background (the `isBroad` machinery in
`dat2mat_svsNAD` is generic), or a tighter NADH4 center row such as
`[8.90 9.00]`. Until then, treat calf NADH4 as the least trustworthy of the six.

### What this says about keeping HSVD where it is

It is a good seeder and a poor quantifier. Where the line really is Lorentzian
and isolated it reproduces VARPRO to a few percent; where lineshape, baseline or
SNR matter it does not, and it loses test–retest precision on the one cohort
with four repeat pairs. Nothing here argues for changing the coarse HSVD (§21)
or for quantifying from it.
