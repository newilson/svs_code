# processU01 31P — fitting notes

Record of 31P fitting tests on the U01 TE 0.1 ms scans (8 brain, 6 calf).
Settings in the drivers get a bare line; the reasoning lives here.

## 1. Pre-fit phase and the NAD+/NADH phase constraint (Oct 12–13 2026)

### 1.1 Starting point

- Calf NAD+ phi0 sat at or near the +20 deg bound in production fits. The
  question was whether the NAD+/NADH phase and shift bounds were too tight.
- Test: `_p31drop/test_bounds_calf.m` (report `_p31drop/bounds/test_bounds_calf.txt`),
  6 calf scans, best of 8, only the NAD bounds changed.
  - Widening the NAD phase bound to +-45/+-90: NAD+ moved to +25..+49 deg,
    NADH phase became undetermined (-77..+45), cost gain <= 0.7%, and 12 of
    24 widened (nested) fits ended above the production cost (missed minima).
    NAD+/PCr moved -40..+33% per scan.
  - Shift: +-12 Hz is about half the 23 Hz NADH–NAD+ gap; at +-24 Hz NADH
    runs onto NAD+. The shift bound is what keeps them apart.
  - The rail was not the NAD bound itself: the whole tissue spectrum sat at
    phi0 +8..+18 deg, leaving NAD+ about 5 deg of room.

### 1.2 Why the tissue sat at +8..+18 deg

- Sign: the model is basis * exp(-i*phi0), so phi0 = +12 means the data sit
  at -12 deg.
- The peak-based pre-fit phase read `angle()` at the sampled magnitude
  maximum. With 1024 points over 9.6 kHz (9.4 Hz/point) and PCr 11 Hz wide in
  calf (15 Hz brain), the line center is up to half a point off the sample,
  where the line already has dispersion: the reading was off by 18–30 deg.
  gATP and aATP are doublets, and the least-squares line through the three
  readings spread their errors into PCr.
- Check: the 16x zero-filled PCr phase agreed with the fit within 4 deg
  (calf -9.7/-17.7 vs fit -11.6/-17.5).
- The "after peak-based correction" figure was also unphased: the step keeps
  the delay phase 2*pi*(ppm - adcshift)*f0*tFirst in the data on purpose (the
  fit models it), and the figure skipped the rotation into the absorptive
  frame. Fixed: `delayPhase` is defined once after tFirst, and the figure is
  drawn after the peak-based and manual phase steps in the fit's frame.

### 1.3 Changes kept in dat2mat_svs31p

- Phase read on the 16x zero-filled, line-broadened spectrum after frequency
  correction (not before it); the first-order line is fitted in ppm.
- Frequency correction shifts with `t + tFirst`, so a shifted line carries
  the delay phase of its new position (0.5–1.5 deg here).
- Drivers keep the 3-peak (PCr/gATP/aATP) first-order pre-phase ("ZF3").

### 1.4 Arms tested (`_p31drop/test_prephase.m`, reports in `_p31drop/prephase/`)

All 14 scans, driver pars captured with portshim, best of 8. OLD = the code
before 1.3 (`_p31drop/dat2mat_svs31p_oldphase.m`).

| Arm | What | Outcome |
|---|---|---|
| ZF3 | 1.3 as kept | brain PCr phi0 -3..+5 (was -19..+16); calf still +10..+14 (doublet readings still biased); cost ~ OLD |
| NEW | zero order from PCr only, slope left to the fit's phi1 | PCr centered (-5.5..+1.6) but calf cost +1.8..+8.7% in all 6: phi1 rotates the peaks, not the causal baseline, so ~3 deg/ppm stays in the data the baseline sees |
| H3/H4 | HSVD coarse phase (1H-style; K 30; delay phase removed per component; PCr/gATP/aATP, +Pi for H4) | HSVD splits each 31P line into opposite-phase components; the summed phase moves with K (calf PCr -63/-85/-80 deg at K 30/60/100 vs fit -69). Sweep (`sweep_hsvd_phase.txt`): best setting RMS vs the fitted phase line brain 5/16, calf 9/14 deg (at 0 / -8.3 ppm) vs ZF 3/24, 12/14, with >100 deg outliers. Full fits: PCr scattered -16..+15, calf cost worse than ZF3 in all 6 (up to +9%). Rejected |
| G_ZF3/G_NEW | in-fit global phase: phiG0 + phi1 rotate the data, so the baseline too; per-metabolite phi0 relative to it | Frame largely independent of the pre-step (within 3 deg in 11/14), cost below ZF3 in 10/14, but the frame follows the baseline: calf PCr +4..+14, brain -17..+15 relative to it, so bounds still not centered on the tissue |
| *_A | every phi0 anchored to PCr (offsets +-20) | With G: degenerate — the baseline takes its own phase (phiG0 -83..+74, PCr compensates), data shown 50–70 deg rotated, brain gATP/PCr 0.94 -> 0.58. Without G: NAD+ moves to +16..+20 from PCr (~+30 abs), NAD+/PCr down to -24%. Rejected. Report: `test_prephase_v4_anchor.txt` |
| ZF3_T | ZF3 + NADplus/NADH phi0 tied to aATP | **Adopted.** NAD phase rail gone in all 14 (shared phase calf +12..+17, brain -14..+5); calf cost below ZF3 in all 6 with 2 fewer parameters (ZF3's search was missing minima); brain within +-0.6 points |
| G_NEW_T | G_NEW + the same tie | Lower cost than ZF3_T in brain 6/8, higher in calf 4/6, more NAD+ movement; not worth the model change |

Common thread: whenever NAD+ gets more phase room it moves to the edge, its
amplitude drops, and the cost barely changes. Its phase is not determined by
the data, so it is tied to its well-determined neighbor aATP rather than left
to a bound. This tie is narrower than the 17-metabolite tissue phase tie that
was rejected on Oct 12 (separate phi0 per metabolite otherwise kept).

NAD+/PCr x100, ZF3 -> ZF3_T: calf 1.192->0.878, 1.070->0.899, 1.471->1.378,
1.200->1.211, 1.275->1.408, 1.253->1.135; brain -19..+5%.

### 1.5 UDP-sugar phases (`_p31drop/test_udp_phase.m`, reports in `_p31drop/udpphase/`)

With ZF3_T in the drivers (PROD; reproduced the ZF3_T arm exactly), the four
UDP-sugar phases were free: 24 of 56 sat at +-20 deg, in 11 of 14 scans, with
opposite signs for neighboring sugars. The UDP alpha-P lines (-8.0..-8.3 ppm)
sit under NAD+/NADH, so free UDP phases act directly on the NAD region.

| Arm | What | Cost vs PROD | Notes |
|---|---|---|---|
| UDP_A | UDP phi0 in the aATP/NADplus/NADH phase group | +0.04..+0.97% (mean 0.44) | Calf forced off its preferred phase; NADH shift unpinned in 3 calf scans but NADH roughly halved in two (S139 0710 0.124 -> 0.048) |
| UDP_S | one shared phase for the four UDP (free in +-20) | -0.00..+1.20% (mean 0.33) | **Adopted.** Cheaper than UDP_A in 9/14 (calf 5/6) |

- With its own phase, UDP sits 11–25 deg below aATP in every calf scan; in
  brain within 9 deg of aATP in 7/8. A consistent 20 deg offset in calf
  suggests the UDP basis may be taking a feature that is not UDP (open).
- NAD+/PCr x100 PROD -> UDP_S mostly within +-10%; S007 brain 3.621 -> 4.530,
  calf 859875 09/01 0.878 -> 1.084. Repeat pairs did not improve with either
  tie (n = 2, not a basis for choosing).

Drivers now: widths of aATP/NAD+/NADH/UDP tied (lbL, lbG); phi0 of NAD+/NADH
tied to aATP; one shared phi0 for the four UDP.

### 1.6 Removed options (Oct 13 2026)

Tested and taken out of the code; the arms' saved results remain in
`_p31drop/prephase/`.

- `pars.hsvdPhase` in dat2mat_svs31p (HSVD coarse phase). The shared HSVD
  helpers moved for it (`utils/hsvdComponents.m`, `hsvdAssign.m`,
  `peakLinearPhase.m`) stay in utils for dat2mat_svsNAD; the move was checked
  identical on 859875 brain 1H (`isequaln` of the full output).
- `fitOpt.globalPhase` / `phaseG0Init` / `phaseG0Bounds` in
  curvefitAuto_basisVarpro (global phase applied by rotating the data).
- `fitOpt.phaseRef` / `phaseRefBounds` (phi0 anchored to a reference).

### 1.7 Open

- UDP's own phase 11–25 deg below aATP in calf (1.5).
- NADH shift sits at -10..-12 Hz (toward NAD+) in nearly every scan and arm,
  pinned at the -12 Hz bound in 6/14 with ZF3_T — NADH basis position, or
  NADH taking a different feature.
- Fit fragility: the same ZF3 fit run in the MATLAB client (multithreaded)
  instead of a pool worker gave calf 859875 10/06 NAD+ 0.908 vs 1.070 at a
  0.5% lower cost. Compare options only on the same execution path, best of N.
  Production rerun of the four 859875 sessions with ZF3_T in the drivers
  (ps1 pipeline, client MATLAB) vs the ZF3_T test arm (pool workers): NAD+/PCr
  x100 brain 7.452 vs 7.279 and 7.077 vs 7.089, calf 0.948 vs 0.878 and
  0.919 vs 0.899.
- Brain 859875 ext/int: 5.26 (09/01) vs 12.05 (10/06), phantom SNR 12.9 vs 7.4.
