# Theoretical Framework for CH₄ Clumped Isotopologue Thermal Calibration

*A working note summarising the experimental setup, theoretical choices, and their practical implementation. All numbers below are recomputed from `calibration_data_table.csv` by `make_AMT_outputs.py` (Good-only revision); nothing is carried over from stale output files.*

---

## 1. What we are trying to measure

Methane contains rare doubly-substituted isotopologues, ¹³CH₃D and ¹²CH₂D₂, whose abundances relative to a random (stochastic) isotope distribution are reported as Δ¹³CH₃D and Δ¹²CH₂D₂ (in per mille, ‰). For an isotopically equilibrated system, these Δ values decrease monotonically with temperature and can therefore serve as a geothermometer.

Calibrating the thermometer requires producing CH₄ in known isotopic equilibrium across a wide temperature range and anchoring those measurements to a physically rigorous theoretical prediction. Three separate things must therefore agree: the experiments, the mass-spectrometric measurements, and the theory.

**Data policy (Good-only).** Rows flagged `Bad` in the CSV are dropped at load and appear in no fit, figure, or table. Of the 36 CSV rows, 30 Good rows are used; the 6 dropped Bad rows are: CAL1549 at 450 °C, 50 °C, and 2 × 25 °C, and AP613 2025 at 4 °C and −20 °C (~3 months). One Good row — CAL1549 at −20 °C — is additionally *withheld from every fit* (it is the only row carrying a raw value, and its corrected−raw, ≈ −6.4 / −8.4 ‰, is not the working-gas shift); it stays in the tables and figures as an open symbol.

## 2. The experimental heating setup

AP613 and CAL1549 — commercially pure CH₄ reference cylinders with known bulk isotopic compositions — are used as working gases. For each calibration point, approximately 10 catalyst pellets are loaded into a 20 mL borosilicate glass tube fitted with a yellow Teflon stopper (or a flame-sealed borosilicate nozzle for multi-month runs). The tube is evacuated to 10⁻³ mbar, oxygen-burned at 550 °C for 5 h to activate the catalyst, re-evacuated overnight, and then filled with pure CH₄ before being set to the target equilibration temperature.

Two catalysts cover the full temperature range:

- γ-Al₂O₃ for temperatures below 200 °C (as low as −20 °C in our protocol)
- Pt/Al₂O₃ for 200–550 °C

Once equilibration is complete, the tube is connected to the High-Concentration Extraction System (HCES). The extracted CH₄ is passed through two GC columns (5 Å molecular sieve + HayeSep D), cryogenically collected on silica gel, and the purity verified by the TCD chromatogram — negligible O₂ and N₂ peaks confirm no atmospheric contamination. The purified gas is then transferred to the 253 Ultra isotope-ratio mass spectrometer for measurement of Δ¹³CH₃D and Δ¹²CH₂D₂.

Calibration points span −20 to 550 °C. Points below 270 K (−3.15 °C) have not been shown to reach equilibrium; they are treated as a **kinetic constraint**, not as calibration points (Section 7).

## 3. The two theoretical models

### 3.1 BMU — Bigeleisen-Mayer-Urey (1947)

BMU computes isotopologue partition function ratios (PFRs) analytically using a harmonic oscillator potential (V = ½kr²) and classical rotational statistics. Individual PFRs deviate substantially from the exact answer: up to 95 ‰ for D/H-related exchange and 5–6 ‰ for ¹³C-related exchange relative to PIMC. Crucially, however, these errors largely cancel when computing the equilibrium constant K, leaving Δ values that agree with PIMC to within ≤ 0.1–0.4 ‰ for these specific symmetric exchange reactions.

The cancellation is a coincidence of reaction symmetry, not a fundamental property of BMU, and it breaks down for cross-molecular fractionation (e.g., CH₄ vs H₂) and worsens at low temperatures where quantum effects grow. BMU is nonetheless practically useful because it is valid at all temperatures and can be evaluated quickly.

### 3.2 PIMC — Path Integral Monte Carlo

PIMC (Eldridge et al. 2019, building on Webb & Miller 2014) maps the quantum mechanical partition function onto a classical ring-polymer and samples it by Monte Carlo integration using an anharmonic ab initio potential energy surface computed at the CCSD(T)/cc-pVTZ+QZ level. It naturally captures vibrational anharmonicity, quantised rotation, and the full zero-point energy without approximation. Errors on individual Δ values are ≤ 0.03 ‰ for Δ¹³CH₃D and ≤ 0.33 ‰ for Δ¹²CH₂D₂.

PIMC is the ground-truth model, but it has been computed only over **270–800 K** (−3 to 527 °C). Its predictions are provided as seventh- and sixth-order polynomials in 1/T forced through 0 ‰ at infinite temperature (Eldridge et al. 2019, equations 10 and 11):

- Δ¹³CH₃D: 7th-order polynomial in 1/T (coefficients 1.47348e19 … 8.89370e2)
- Δ¹²CH₂D₂: 6th-order polynomial in 1/T (coefficients −9.67634e15 … 1.86258e3)

### 3.3 Why not improve BMU with DFT corrections?

Liu & Liu (2016) applied ex post facto anharmonic corrections to harmonically computed BMU RPFRs (using MP2/aug-cc-pVTZ). Individual PFR errors were substantially reduced, yet the final Δ values became *worse* by 0.2–0.4 ‰ relative to PIMC. The corrections disrupted the precise error cancellation without fully correcting its source. BMU corrections remain reaction-specific and structurally unreliable; PIMC is the only approach that is correct at the level of the partition function itself.

## 4. The reference curve: linking PIMC and BMU below 270 K

Our calibration covers −20 to 550 °C, but the PIMC polynomials are valid only down to 270 K. The reference curve used throughout is therefore **hybrid**:

- **T ≥ 270 K:** the PIMC polynomial itself.
- **T < 270 K:** a 5th-order polynomial in 1/T with **no constant term** (the anomaly is 0 at infinite temperature), least-squares fitted to the Eldridge et al. (2019) **Table 4 BMU values** (32 points, −3.1 to 526.9 °C; fit RMS 0.0013 ‰ for Δ¹³CH₃D and 0.0035 ‰ for Δ¹²CH₂D₂), shifted by a constant bridge
  **δ = mean(PIMC − BMU fit) over every integer kelvin from 283 to 423 K**: δ = **+0.0494 ‰ (Δ¹³CH₃D)**, **+0.0547 ‰ (Δ¹²CH₂D₂)**.

> Note: an earlier version of this note bridged PIMC ↔ BMU with a ln–ln linear map fitted on 283–423 K. That mapping is superseded. The current implementation fits the BMU Table-4 values directly and shifts by the grid-mean offset, which keeps the low-T branch anchored to the published BMU numbers rather than to a regression between two models.

**Reference values at −20 °C (253.15 K):** Δ¹³CH₃D = **7.527 ‰**, Δ¹²CH₂D₂ = **27.807 ‰** (the 283–423 K bridge). The extrapolated PIMC polynomial would give 7.766 / 27.243 ‰.

**Known systematic — the 270 K join is not continuous.** Because PIMC − BMU is temperature-dependent but δ is a constant, the reference curve steps at 270 K by **+0.019 ‰ (Δ¹³CH₃D)** and **−0.374 ‰ (Δ¹²CH₂D₂)** (above minus below). The Δ¹²CH₂D₂ step is comparable to the offset standard errors; residuals, derivatives, and curve inversions across the join must be read with this in mind. The step is printed at runtime and recorded in the workbook manifest.

**Theory envelope at 253.15 K** (bracketing the bridge choice — grid-mean bridges on 283–423, 270–400, and 300–450 K windows, a value-match at 270 K, and the extrapolated PIMC polynomial):

- Δ¹³CH₃D: **7.527 – 7.766 ‰** (spread 0.24 ‰)
- Δ¹²CH₂D₂: **27.243 – 27.855 ‰** (spread 0.61 ‰)

The local slope of the reference curve at 253.15 K is **−0.0460 ‰/K (Δ¹³CH₃D)** and **−0.2255 ‰/K (Δ¹²CH₂D₂)**, so a ±2 °C temperature error moves the reference by ±0.09 ‰ / ±0.45 ‰ there (Table S12) — much larger than the ≤ 0.23 ‰ sometimes quoted for Δ¹²CH₂D₂ at warmer temperatures.

## 5. Frames and the working-gas offset correction

Measurements are made relative to an unheated cylinder gas, so raw Δ values carry the cylinder's own anomaly as a systematic offset. Two frames are used and must not be mixed:

- **Working-gas (raw) frame** — AP613 2021 (Sivan et al. 2024) and CAL1549 are stored and fitted in this frame.
- **Absolute frame** — AP613 2025 and the Eldridge means are stored in this frame; AP613 extended (both campaigns pooled) is absolute.

The 2021 raw→absolute shift is **measured from the CSV**: corrected − raw = **+2.2633 ‰ (SD 0.005, n = 6)** for Δ¹³CH₃D and **+2.7450 ‰ (SD 0.010, n = 6)** for Δ¹²CH₂D₂. It is treated as a fixed frame definition, not a random error. The AP613 2025 working-gas offset is its absolute offset **minus this shift**; an alternative that uses −(2021 Method-1 offset) as the shift is reported but not used — for Δ¹²CH₂D₂ the two shifts differ by **0.51 ‰** (for Δ¹³CH₃D they agree to 0.0001 ‰).

**Method-1 offsets (weighted mean of Δ − reference, slope fixed at 1), Good points at T ≥ 270 K only, SE inflated by √(χ²ᵣ) when χ²ᵣ > 1:**

| Data set | Frame | Δ¹³CH₃D (‰) | Δ¹²CH₂D₂ (‰) | n | χ²ᵣ (13 / 12) |
|---|---|---|---|---|---|
| AP613 2021 | working gas | **−2.263 ± 0.115** | **−3.260 ± 0.817** | 6 | 0.32 / 0.80 |
| AP613 2025 | absolute | **−0.312 ± 0.274** | **−0.084 ± 1.250** | 5 | 2.66 / 1.17 |
| AP613 2025 | working gas (abs − shift) | **−2.575 ± 0.274** | **−2.829 ± 1.250** | 5 | — |
| ↳ alternative shift (−2021 offset) | working gas | −2.575 | **−3.344 ± 1.495** (reported, not used) | 5 | — |
| CAL1549 | working gas | **−5.876 ± 0.197** | **−8.306 ± 0.742** | 6 | 2.94 / 0.62 |
| AP613 extended (consistency) | absolute | −0.100 ± 0.115 | −0.370 ± 0.667 | 11 | 1.46 / 0.88 |
| Eldridge et al. 2019 means (consistency) | absolute | −0.018 ± 0.071 | +0.001 ± 0.440 | 11 | 7.28 / 7.13 |

The AP613 offsets in 2021 and 2025 are not systematically different: on the absolute scale the 2025−2021 campaign difference at T ≥ 270 K is **−0.311 ± 0.297 ‰ for Δ¹³CH₃D (z = −1.05, p = 0.30)** and **+0.430 ± 1.494 ‰ for Δ¹²CH₂D₂ (z = 0.29, p = 0.77)**; the ANCOVA intercept difference is −0.430 ± 0.239 ‰ (p = 0.11). A drift of a few tenths of a permil is the central value and is unresolved — calibration stability over four years is confirmed at that level.

### 5a. Measurement scales: what is actually measured

- **Eldridge et al. (2019)** and **Sivan et al. (2024)** are independent studies that reported anomalies on the **absolute (thermodynamic) scale**. Eldridge's house gas was assigned from the PIMC curve itself, so its offset is ≈ 0 by construction; the Eldridge rows are a consistency check, never a test of PIMC.
- **AP613 2021 (Sivan et al. 2024)**: the CSV carries both frames — `D13_raw`/`D12_raw` (instrument working-gas frame) and `D13`/`D12` = raw + 2.2633 / +2.7450 ‰ (absolute; the shift measured on these six rows). The raw-frame Method-1 offset (−2.263 / −3.260 ‰) is the AP613 cylinder's clumped anomaly on the absolute scale.
- **This study (2025)**: the new AP613 tubes (Good points: 50, 100, 150 °C plus two 550 °C tubes) and CAL1549 (Good points: 50–450 °C, the full range) were measured with **Sivan's AP613 as the instrument reference gas**, so the 2025 values are relative to AP613 — absolute values were not measured. Their "absolute" column is derived by adding the 2021 shift, so the 2025 working-gas offset is conditional on that shift: for Δ¹³CH₃D the CSV shift (+2.2633) equals −(2021 offset) to 0.0001 ‰, but for Δ¹²CH₂D₂ it differs by 0.51 ‰ (both variants are reported in the manifest; only the CSV shift is used).

### 5b. The three AP613 interpretations

| version | Good points | offset Δ¹³CH₃D (frame) | slope Δ¹³CH₃D | role |
|---|---|---|---|---|
| AP613 2021 only (Sivan) | 6, 50–400 °C, raw frame | −2.263 ± 0.115 (wg) | 1.062 ± 0.102 | primary — the 2021 anchor |
| AP613 2025 only (this study) | 5, 50–550 °C (two 550 °C tubes), absolute frame | −0.312 ± 0.274 abs = −2.575 ± 0.274 wg | 1.173 ± 0.171 | primary — weakest alone (550 °C bulk-composition outlier: dropping both moves m by −0.28; Δ¹² CI 0.26–2.11) |
| AP613 extended (pooled) | 11, −20–550 °C, absolute frame | −0.100 ± 0.115 (abs) | 1.067 ± 0.085 | consistency / continuity statement — nests the 2021 anchor, so not an independent slope test |

Continuity of the Utrecht reference gas is shown by the campaign comparison, not by pooling: 2025 − 2021 at T ≥ 270 K is −0.311 ± 0.297 ‰ (Δ¹³CH₃D, z = −1.05, 95 % CI −0.89…+0.27) and +0.430 ± 1.494 ‰ (Δ¹²CH₂D₂, z = 0.29); zero lies inside both intervals, so pooling is statistically permissible, but a drift of a few tenths of a permil is not excluded. **Recommendation:** keep the three gases separate in the six-test primary family, quote the campaign difference as the continuity evidence, and use AP613 extended only as the pooled consistency fit.

## 6. Regression against theory and key results

Offset-corrected values are regressed against the reference curve with two weighted least-squares methods (measurement SEs only — the ±2 °C temperature error is **not** in the fit weights; its effect is quantified separately in Table S12):

- **Method 1:** offset only, slope fixed at 1 (Section 5).
- **Method 2:** weighted line Δ_meas = m·Δ_ref + c on the same points, with t-based 95 % CIs and two-sided p for m = 1. SEs are inflated by √(χ²ᵣ) when χ²ᵣ > 1. The point-estimate identity c = offset − (m−1)·x̄_w is checked at runtime (residual < 1e-8).

**Method-2 slopes at T ≥ 270 K (Good points only):**

| Data set | m ± SE(Δ¹³CH₃D) | z, p | m ± SE(Δ¹²CH₂D₂) | z, p |
|---|---|---|---|---|
| AP613 2021 (primary) | 1.062 ± 0.102 | +0.61, 0.58 | 0.794 ± 0.147 | −1.40, 0.24 |
| AP613 2025 (primary) | 1.173 ± 0.171 | +1.01, 0.39 | 1.185 ± 0.291 | +0.64, 0.57 |
| CAL1549 (primary) | 0.894 ± 0.145 | −0.73, 0.50 | 1.009 ± 0.154 | +0.06, 0.96 |
| AP613 extended (consistency) | 1.067 ± 0.085 | +0.79, 0.45 | 0.899 ± 0.127 | −0.80, 0.45 |
| Eldridge means (consistency) | 1.037 ± 0.031 | +1.23, 0.25 | 0.872 ± 0.053 | −2.41, 0.04 |

The confirmatory family is the 6 primary tests (3 working gases × 2 isotopologues); adjusted p = min(1, 6p). The smallest unadjusted primary p is 0.235, so the smallest adjusted p is **1.00** — **none of the six primary slopes is distinguishable from 1 after the pre-declared correction; several intervals are wide (AP613 2021 Δ¹²CH₂D₂ 95 % CI 0.385–1.203), and a CI that contains 1 is not a demonstration that the slope is 1**. All-temperature (all-Good) and ln-space refits are sensitivities, not members of the confirmatory family. Inverse-variance pools of the primary slopes: Δ¹³CH₃D **1.038 ± 0.075**, Δ¹²CH₂D₂ **0.931 ± 0.100** (both consistent with 1). The Eldridge Table-2 mean-fit slopes (refit here) do not reproduce the published Δ¹²CH₂D₂ slope 0.98 ± 0.05, which was fit to individual analyses, and because Eldridge's house-gas anomalies were assigned from this same PIMC curve, the Eldridge rows are a consistency check, not a test of the curve (Table S9: working-gas offsets recomputed here 2.61 ± 0.07 / 5.87 ± 0.45 ‰ vs published 2.59 ± 0.14 / 5.86 ± 0.6 ‰).

**Sensitivities (Table S4, S11):** refits using all Good points including T < 270 K, the extrapolated-PIMC low-T curve, and 1000·ln(1 + Δ/1000) space; leave-one-out of the two Good AP613 −20 °C points (dropping both moves the extended Δ¹³CH₃D slope from 0.868 to 1.067; dropping one moves it ≤ 0.06 — the two single-point changes do not add); and dropping the 550 °C tubes (different bulk composition) moves the AP613 2025 Δ¹³CH₃D slope by −0.28 and its offset by +0.21 ‰.

## 7. Low-temperature results (−20 °C): a kinetic constraint, not a calibration point

Good AP613 2025 aliquots at −20 °C (253.15 K, absolute frame): Δ¹³CH₃D = 6.25 ± 0.23 and 6.65 ± 0.26 ‰; Δ¹²CH₂D₂ = 24.55 ± 1.2 and 26.88 ± 2.1 ‰ (the two rows share one duration string, "~5 months", and may be replicates of one tube). Against the reference at 253.15 K plus the 2025 offset, the mean residuals are **−0.79 ± 0.32 ‰ (Δ¹³CH₃D, z = −2.4)** and **−2.60 ± 1.63 ‰ (Δ¹²CH₂D₂, z = −1.6)**. Under a monotonic approach from below, an incomplete experiment is only a **lower bound** on the equilibrium anomaly: the recommended bound is **6.08 ‰** for Δ¹³CH₃D (higher-aliquot bound 6.08, replicate-mean bound 6.09), which inverts to an apparent temperature of **+16.5 °C** on the reference curve; the reference curve equals the extrapolated PIMC(253.15 K) at −25.1 °C. The CAL1549 −20 °C row is withheld from all fits: raw 6.56 / 28.23 ‰, fitted-absolute 6.06 / 28.11 ‰.

The Eq. (5) fraction of the cylinder-to-reference distance covered (working-gas frame, cylinder = 0) at −20 °C: AP613 **0.805 ± 0.064** and **0.886 ± 0.072** (Δ¹³CH₃D), **0.873** and **0.966** (Δ¹²CH₂D₂); CAL1549 (withheld) 0.109 ± 0.297 (Δ¹³CH₃D), 1.015 (Δ¹²CH₂D₂).

## 8. Temperature uncertainty of the thermometer (Table S8)

The 68 % temperature interval is obtained by Monte Carlo inversion of the reference curve: draw Δ ~ N(ref(T), √(prec² + SE_offset²)) 20 000 times, invert the curve on a 200–1000 K grid, and take the 15.865/84.135 percentiles. Precision variants: long-term reproducibility 0.15 ‰ (Δ¹³CH₃D) / 1.2 ‰ (Δ¹²CH₂D₂, Sivan et al. 2024) and internal precision 0.3 / 2.4 ‰. Example (AP613 extended, reproducibility, Δ¹³CH₃D): at true T = 50, 150, 300 °C the 68 % half-widths are **−7.0/+7.5, −12.6/+14.0, −26.6/+31.1 °C**; Δ¹²CH₂D₂ at 50 °C: −11.1/+12.6 °C. The intervals are asymmetric because the curve is nonlinear and flattens at high T.

## 9. Outputs and how to run

```
python3 make_AMT_outputs.py [--csv calibration_data_table.csv] [--out AMT/] [--json results_summary.json]
```

Requires numpy, matplotlib, openpyxl. The optional JSON is used **only** for the published Eldridge benchmark numbers in Table S9; the reference curve, offsets, slopes, and low-temperature bounds are always recomputed from the CSV, so a stale JSON cannot move the curve.

Figures (PDF + PNG, 300 dpi): `Fig1_Delta_vs_T`, `Fig2_residuals_vs_T`, `FigS1_method2_regression`, `FigS2_shifted_by_offset`, `FigS3_low_temperature`. Tables: `AMT_tables.xlsx` with a Manifest sheet plus Table 1 (calibration summary), Table 2 (experiments ≤ 25 °C with residuals, offset-aware z, Eq. 5 fractions), and Tables S1–S12 (protocol/durations, Method-1 offsets, all-Good Method-2, four slope sensitivities, campaign comparison, low-T inferences, withheld rows, Monte-Carlo temperature intervals, Eldridge benchmark, pooled slopes, leverage sensitivities, reference-curve slopes).

In figures, filled symbols are Good points in the T ≥ 270 K offset fit; open symbols are Good points below 270 K or the withheld CAL1549 −20 °C row; the grey band marks T < 270 K.

---

### References

- Bigeleisen, Mayer & Urey (1947) J. Chem. Phys. 15:261
- Eldridge et al. (2019) ACS Earth Space Chem. 3:2747–2764
- Webb & Miller (2014) J. Phys. Chem. A 120:467–474
- Liu & Liu (2016) GCA 175:252–270
- Sivan et al. (2024) AMT 17:2687–2705
