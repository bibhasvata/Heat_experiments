Theoretical Framework for CH4 Clumped Isotopologue Thermal Calibration
A working note summarising the experimental setup, theoretical choices, and their practical implementation
Dasgupta, Eckhardt, Sivan, Gilbert & Röckmann — IMAU, Utrecht

1. What we are trying to measure
Methane contains rare doubly-substituted isotopologues, ¹³CH₃D and ¹²CH₂D₂, whose abundances relative to a random (stochastic) isotope distribution are reported as Δ¹³CH₃D and Δ¹²CH₂D₂ (in per mille). For an isotopically equilibrated system, these Δ values decrease monotonically with temperature and can therefore serve as a geothermometer.
Calibrating the thermometer requires producing CH₄ in known isotopic equilibrium across a wide temperature range and anchoring those measurements to a physically rigorous theoretical prediction. Three separate things must therefore agree: the experiments, the mass-spectrometric measurements, and the theory.
2. The experimental heating setup
AP613 and CAL1549 — commercially pure CH₄ reference cylinders with known bulk isotopic compositions — are used as working gases. For each calibration point, approximately 10 catalyst pellets are loaded into a 20 mL borosilicate glass tube fitted with a yellow Teflon stopper (or a flame-sealed borosilicate nozzle for multi-month runs). The tube is evacuated to 10⁻³ mbar, oxygen-burned at 550 °C for 5 h to activate the catalyst, re-evacuated overnight, and then filled with pure CH₄ before being set to the target equilibration temperature.
Two catalysts cover the full temperature range:
•	γ-Al₂O₃ for temperatures below 200 °C (as low as −20 °C in our protocol)
•	Pt/Al₂O₃ for 200–550 °C
Once equilibration is complete, the tube is connected to the High-Concentration Extraction System (HCES). The extracted CH₄ is passed through two GC columns (5 Å molecular sieve + HayeSep D), cryogenically collected on silica gel, and the purity verified by the TCD chromatogram — negligible O₂ and N₂ peaks confirm no atmospheric contamination. The purified gas is then transferred to the 253 Ultra isotope-ratio mass spectrometer for measurement of Δ¹³CH₃D and Δ¹²CH₂D₂.
3. The two theoretical models
3.1 BMU — Bigeleisen-Mayer-Urey (1947)
BMU computes isotopologue partition function ratios (PFRs) analytically using a harmonic oscillator potential (V = ½kr²) and classical rotational statistics. Individual PFRs deviate substantially from the exact answer: up to 95‰ for D/H-related exchange and 5–6‰ for ¹³C-related exchange relative to PIMC. Crucially, however, these errors largely cancel when computing the equilibrium constant K, leaving Δ values that agree with PIMC to within ≤0.1–0.4‰ for these specific symmetric exchange reactions.
The cancelation is a coincidence of reaction symmetry, not a fundamental property of BMU, and it breaks down for cross-molecular fractionation (e.g., CH₄ vs H₂) and worsens at low temperatures where quantum effects grow. BMU is nonetheless practically useful because it is valid at all temperatures and can be evaluated quickly.
3.2 PIMC — Path Integral Monte Carlo
PIMC (Eldridge et al. 2019, building on Webb & Miller 2014) maps the quantum mechanical partition function onto a classical ring-polymer and samples it by Monte Carlo integration using an anharmonic ab initio potential energy surface computed at the CCSD(T)/cc-pVTZ+QZ level. It naturally captures vibrational anharmonicity, quantised rotation, and the full zero-point energy without approximation. Errors on individual Δ values are ≤0.03‰ for Δ¹³CH₃D and ≤0.33‰ for Δ¹²CH₂D₂.
PIMC is the ground-truth model, but it has been computed only over 270–800 K (−3 to 527 °C). Its predictions are provided as seventh- and sixth-order polynomials in 1/T forced through 0‰ at infinite temperature (Eldridge et al. 2019, equations 10 and 11).
3.3 Why not improve BMU with DFT corrections?
Liu & Liu (2016) applied ex post facto anharmonic corrections to harmonically computed BMU RPFRs (using MP2/aug-cc-pVTZ). Individual PFR errors were substantially reduced, yet the final Δ values became worse by 0.2–0.4‰ relative to PIMC. The corrections disrupted the precise error cancelation without fully correcting its source. BMU corrections remain reaction-specific and structurally unreliable; PIMC is the only approach that is correct at the level of the partition function itself.
4. How we link PIMC and BMU in practice
Our calibration covers −20 to 550 °C. The PIMC polynomials are valid only down to −3 °C (270 K). For temperatures at or below −3 °C, we established a linear mapping between PIMC and BMU (B3LYP/6-311++G**) using data in the overlapping range 283–423 K:
ln(K¹³CH₃D)ᴘᴵᴹᴶ = 1.03372 × ln(K¹³CH₃D)ᴮᴹᵁ + 0.07769 × 10⁻³
ln(⁸⁄₃ K¹²CH₂D₂)ᴘᴵᴹᴶ = 0.94297 × ln(⁸⁄₃ K¹²CH₂D₂)ᴮᴹᵁ + 0.02904 × 10⁻³
BMU values below −3 °C are passed through these relationships to produce PIMC-equivalent extrapolations. This hybrid approach (PIMC where valid, PIMC-corrected BMU elsewhere) mirrors the strategy in Image 1 of the calibration plot, where the “Hybrid” curve joins BMU(T) for T < 270 K to PIMC(T) for T ≥ 270 K.
5. Working gas offset correction
All measurements are made relative to an unheated cylinder gas (AP613 or CAL1549). Because neither working gas is in stochastic isotopic distribution (Δ = 0‰), raw experimental Δ values carry a systematic offset. This offset is determined by fitting Monte Carlo resampled polynomials to the full temperature series of raw data, then subtracting it.
Offsets from our fits:
•	AP613 (2021): Δ¹³CH₃D = −2.22 ± 0.13‰; Δ¹²CH₂D₂ = −3.2 ± 0.9‰
•	AP613 (2025): Δ¹³CH₃D = −2.48 ± 0.20‰; Δ¹²CH₂D₂ = −2.2 ± 1.3‰
•	CAL1549 (2025): Δ¹³CH₃D = −5.75 ± 0.13‰; Δ¹²CH₂D₂ = −8.2 ± 0.8‰
AP613 offsets in 2021 and 2025 are not systematically different, confirming calibration stability over four years.
6. Regression against theory and key results
Offset-corrected experimental values are plotted against the PIMC polynomial predictions (and their PIMC-corrected BMU extrapolations below −3 °C). Linear regression using YorkFit (York et al. 2004) — which propagates uncertainties in both axes including an assumed ±2 °C temperature error — yields slopes and intercepts for each working gas and year.
Key observations:
•	Slopes are within 3σ of unity for all individual datasets, indicating no significant deviation from the theoretical temperature dependence.
•	Intercepts reflect the working gas offset and are consistent between 2021 and 2025 for AP613.
•	No systematic shift between the two experimental campaigns is detected, confirming long-term reproducibility of the equilibration protocol.
Together, these results demonstrate that our revised protocol for −20 to 550 °C produces equilibrated gases fully consistent with PIMC-based theory, and that the calibration is stable across laboratories and years.

References: Eldridge et al. (2019) ACS Earth Space Chem. 3:2747–2764 • Webb & Miller (2014) J. Phys. Chem. A 120:467–474 • Liu & Liu (2016) GCA 175:252–270 • Sivan et al. (2024) AMT 17:2687–2705 • York et al. (2004) Am. J. Phys. 72:367
