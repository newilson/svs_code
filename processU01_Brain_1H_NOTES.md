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
H6/H2, NADH2 CV, and NADH2 agreement — not on H4/H2.** Keep NADH4 as a MODELLED
peak in region 1 regardless (removing it leaks its signal into H6 and the
baseline); just don't report it. This also settles the 8.6-vs-8.7 debate in
favour of 8.7.

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
is centred BELOW the window — NAD H2/H6/H4 are not modelled in region 2, and
SpecTickle's component census puts ~105% of the Trp-window baseline in the tails
of components centred below 9.8. A sum of absorption-mode Lorentzian tails
approaching from below is **monotone decreasing** across the window. An interior
minimum requires a negative-amplitude component or residual dispersive content,
neither of which would track the fitted Trp position in every subject. The notch
is peak/baseline degeneracy: the spline digs under Trp so Trp can be larger.

**Lesson worth generalising: a residual gain bought by an inadmissible baseline
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

Region 2's baseline is almost entirely the tail of components centred below
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
- `lambdaWidth` — normalised so 1.0 = "a 100% relative width error costs as much
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

### `coarseMode = 'perRegion'`

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

Replaces the coarse VARPRO; old code kept commented out in
`dat2mat_svsNAD/nadCoarsePhase`.

`K` and `npts` match SpecTickle's defaults (`hsvd_order` 60, `hsvd_np` 1024).

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

**Trp: indole NH proton (1H) at 10.1 ppm — T1/T2 NOT FINALISED.** Jul 2026,
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

Two naming conventions are recognised, because the U01 sessions and the
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
in the session folder. Previous behaviour (required `datfiles\`):

```matlab
datDir = fullfile(SessionDir, 'datfiles');
assert(isfolder(datDir), 'datfiles\\ subfolder not found in %s', SessionDir);
```

---

## 10. Open issues and do-nots

**Do not:**

- Re-select region-2 spline DOF on Trp-window residual, or on agreement with
  SpecTickle (§2).
- Narrow `ph_range` before the coarse phase step is fixed (§4).
- Reconsider the phase-invariant `areasConv` (§8) — settled.
- Report H4/H2 as a tuning target (§1, scope note).
- Treat agreement with SpecTickle as accuracy. There is no gold standard; the
  H2/H6/H4 ratio test is *self-consistency*, not truth.

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
