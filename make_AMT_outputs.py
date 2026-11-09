#!/usr/bin/env python3
"""
make_AMT_outputs.py - static AMT figures and tables from calibration_data_table.csv.
GOOD-ONLY REVISION: rows flagged "Bad" are dropped at load and appear in no
fit, no figure and no table.  The only non-filled symbols left are points
below 270 K and the withheld CAL1549 -20 C row.

Self-contained: the ONLY required input is the CSV next to this script
(override with --csv).  No HTML and no Plotly figures are produced.
results_summary.json is optional and is used only for the published Eldridge
et al. (2019) benchmark numbers in Table S9.  The reference curve, offsets,
slopes and low-temperature bounds are always recomputed from the CSV and the
polynomials below.  A stale JSON cannot move the curve.

Reference curve
  PIMC polynomials of Eldridge et al. (2019) at T >= 270 K.
  Below 270 K, a 5th-order polynomial in 1/T fitted to their Table 4 BMU
  values with no constant term (the anomaly is 0 at infinite T), shifted by
  delta = mean(PIMC - BMU fit) on every integer kelvin from 283 to 423.

What is a calibration, and what is not
  GOOD-ONLY POLICY: Bad (QC-failed) rows are removed immediately after the CSV
  is read.  Every fit, figure and table below uses Good rows only.
  Primary offsets and slopes use Good points at T >= 270 K only.
  Primary slope family (6 tests, Bonferroni): AP613 2021, AP613 2025 and
  CAL1549, both isotopologues.  AP613 extended nests 2021, and the Eldridge
  means were assigned from this PIMC curve, so neither is in that family.
  Points below 270 K are a kinetic constraint.  They are not in the offsets
  or the primary slopes.
  CAL1549 at -20 C is withheld from every fit.  It is the only row with a
  raw value, and corrected-raw is not the working-gas shift of the other
  rows.  The row is Good, so it stays in the tables and figures, but it is
  withheld from every fit.

Frames
  AP613 2021 is fitted in the working-gas frame (raw columns).
  AP613 2025 is stored and fitted in the absolute frame.  Its working-gas
  offset is absolute offset minus the 2021 raw-to-absolute shift measured
  from the CSV.  Figure 1 plots that converted value.  An alternative that
  uses the 2021 Method-1 offset as the shift is reported, not silently used:
  for Delta12 the two shifts differ by several tenths of a permil.
  CAL1549 is working-gas.  Eldridge and AP613 extended are absolute.

Outputs (default: the script's folder; --out to change)
  Fig1_Delta_vs_T, Fig2_residuals_vs_T, FigS1_method2_regression,
  FigS2_shifted_by_offset, FigS3_low_temperature (.pdf and .png)
  AMT_tables.xlsx

Usage
  python3 make_AMT_outputs.py
  python3 make_AMT_outputs.py --csv /path/calibration_data_table.csv --out /path/AMT

Requires: numpy, matplotlib, openpyxl.

Author: Bibhasvata Dasgupta | ICIW 2026
"""

import os, sys, json, csv, math, argparse, re
import numpy as np
import matplotlib
matplotlib.use("Agg")                    # headless: no HTML / interactive output
import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill

# ─────────────────────────────────────────────────────────────────────────────
# 0. ARGUMENTS
# ─────────────────────────────────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser(description="AMT figures+tables from the calibration CSV (no HTML)")
parser.add_argument("--csv", default=os.path.join(HERE, "calibration_data_table.csv"))
parser.add_argument("--out", default=HERE, help="output folder (default: script's folder)")
parser.add_argument("--json", default=None,
                    help="optional results_summary.json (published Eldridge benchmark only)")
parser.add_argument("--dpi", type=int, default=300)
parser.add_argument("--repro13", type=float, default=0.15, help="long-term reproducibility Delta13 (permil), Sivan et al. 2024")
parser.add_argument("--repro12", type=float, default=1.2, help="long-term reproducibility Delta12 (permil), Sivan et al. 2024")
parser.add_argument("--internal13", type=float, default=0.3, help="internal precision Delta13 (permil)")
parser.add_argument("--internal12", type=float, default=2.4, help="internal precision Delta12 (permil)")
parser.add_argument("--mc-draws", type=int, default=20000)
parser.add_argument("--t-grid", default="50,100,150,200,300", help="true temperatures (C) for Table S8")
args = parser.parse_args()

OUT = args.out.rstrip("/")
os.makedirs(OUT, exist_ok=True)

# optional JSON: published Eldridge benchmark only (Table S9).
# Theory stored in that file is ignored, so a legacy curve cannot reset 253 K.
R = None
json_path = args.json or os.path.join(HERE, "results_summary.json")
if os.path.isfile(json_path):
    R = json.load(open(json_path, encoding="utf-8"))
    print(f"[json] {json_path} read for the published Eldridge benchmark only")
else:
    print("[json] no results_summary.json - Table S9 published slopes left blank")

N_PRIMARY = 6   # confirmatory slope tests: 3 working gases x 2 isotopologues

# ─────────────────────────────────────────────────────────────────────────────
# 1. STATISTICS HELPERS (pure math, no scipy)
# ─────────────────────────────────────────────────────────────────────────────
def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 200, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0 / d; h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d
        delta = d * c; h *= delta
        if abs(delta - 1.0) < EPS: break
    return h

def _betai(a, b, x):
    """regularized incomplete beta I_x(a,b)"""
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b

def t_two_sided_p(t, nu):
    """two-sided p of |t| under Student t with nu degrees of freedom"""
    if nu <= 0: return float("nan")
    return _betai(nu / 2.0, 0.5, nu / (nu + t * t))

def t_critical_two_sided(nu, alpha=0.05):
    """Student-t critical value for a two-sided test with the given alpha."""
    lo, hi = 0.0, 1.0
    while t_two_sided_p(hi, nu) > alpha:
        hi *= 2.0
    for _ in range(64):
        mid = (lo + hi) / 2.0
        if t_two_sided_p(mid, nu) > alpha:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0

def inv_norm(p):
    """inverse standard normal CDF (Acklam's rational approximation)"""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1.0 - 0.02425
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)
    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)
    q = p - 0.5; r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0)

# ─────────────────────────────────────────────────────────────────────────────
# 2. INPUT CSV AND FRAME MAP
#    siv : working-gas (raw) frame for its own fits/plots (D13_raw/D12_raw);
#          absolute frame = corrected columns D13/D12 (used for AP613 extended)
#    das : AP613 2025, absolute frame (D13/D12)
#    cal : CAL1549, working-gas frame (D13/D12)
#    eld : Eldridge 2019 Table 2 means, absolute (thermodynamic) frame (D13/D12)
# ─────────────────────────────────────────────────────────────────────────────
if not os.path.isfile(args.csv):
    sys.exit(f"Missing input CSV: {args.csv}")
with open(args.csv, encoding="utf-8") as f:
    ROWS = list(csv.DictReader(f))
print(f"[data] {len(ROWS)} rows read from {args.csv}")
print(f"[out ] figures and tables -> {OUT}")

# GOOD-ONLY FILTER: QC-failed (Bad) rows are dropped here and used nowhere
# below - not in any fit, sensitivity, figure, or table.
BAD_ROWS = [r for r in ROWS if r["flag"].strip().lower() != "good"]
ROWS = [r for r in ROWS if r["flag"].strip().lower() == "good"]

def key_of(dataset):
    if dataset.startswith("CAL1549"): return "cal"
    if dataset.startswith("AP613"): return "das"
    if dataset.startswith("Sivan"): return "siv"
    return "eld"

def num(s):
    s = (s or "").strip()
    return float(s) if s else None

def val_raw(row, iso):    # working-gas frame value
    return num(row["D13_raw"] if iso == "13" else row["D12_raw"]) if key_of(row["dataset"]) == "siv" \
        else num(row["D13"] if iso == "13" else row["D12"])

def val_abs(row, iso):    # absolute frame value
    return num(row["D13"] if iso == "13" else row["D12"])

def val(row, iso):        # value in the frame the data set is fitted/plotted in
    return val_raw(row, iso)

def sev(row, iso):
    return num(row["SE13_used"] if iso == "13" else row["SE12_used"])

def is_good(row):
    return row["flag"].strip().lower() == "good"

def is_withheld(row):
    """Kept in the tables, used in no fit.

    CAL1549 at -20 C is the only row that carries a raw value.  corrected-raw
    is about -6.4 and -8.4 permil, which is not the Sivan working-gas shift
    (+2.26 / +2.74).  The bulk composition of that tube is also an outlier,
    but bulk isotopes are not in this CSV, so the rule is the row itself.
    """
    return key_of(row["dataset"]) == "cal" and abs(float(row["T_C"]) + 20.0) < 1e-9

def in_offset_fit(row):
    """True when the row entered the T >= 270 K Method-1 offset."""
    if not is_good(row) or is_withheld(row):
        return False
    if val(row, "13") is None and val(row, "12") is None:
        return False
    return float(row["T_C"]) + 273.15 >= 270.0

KEYS = ["siv", "das", "cal", "eld"]                 # plotted data sets
ALLKEYS = KEYS + ["apex"]                           # apex = AP613 extended (abs frame)
NAME = {"siv": "AP613 2021 (Sivan et al. 2024)", "das": "AP613 2025 (this study)",
        "cal": "CAL1549 (this study)", "eld": "Eldridge et al. 2019",
        "apex": "AP613 extended (both campaigns, absolute)"}
# Standard display names for the raw CSV dataset strings (e.g. Sivan23 -> AP613 2021)
DS_LABEL = {"CAL1549": "CAL1549 (this study)",
            "AP613 (Dasgupta 2025)": "AP613 2025 (this study)",
            "Sivan23 (2021 campaign)": "AP613 2021 (Sivan et al. 2024)",
            "Eldridge19 (Table 2 means)": "Eldridge et al. 2019 (Table 2 means)"}
def ds_label(r):
    return DS_LABEL.get(r["dataset"], r["dataset"])

print(f"[qc  ] Good-only: kept {len(ROWS)} rows, dropped {len(BAD_ROWS)} Bad row(s):")
for r in BAD_ROWS:
    print(f"[qc  ]   dropped {ds_label(r)}  T={r['T_C']} C  catalyst={r['catalyst']}  flag={r['flag']}")

COLOR = {"siv": "#1f77b4", "das": "#d62728", "cal": "#2ca02c", "eld": "#9467bd"}
MARK = {"siv": "o", "das": "s", "cal": "^", "eld": "D"}
FRAME = {"siv": "working gas (raw)", "das": "absolute", "cal": "working gas",
         "eld": "absolute", "apex": "absolute"}
ISO_LAB = {"13": "Δ¹³CH₃D", "12": "Δ¹²CH₂D₂"}

def rows_of(k):
    if k == "apex":
        return [r for r in ROWS if key_of(r["dataset"]) in ("das", "siv")]
    return [r for r in ROWS if key_of(r["dataset"]) == k]

def get_pts(k, iso):
    """points of data set k in its own frame (apex: absolute frame)"""
    pts = []
    for r in rows_of(k):
        v = val_abs(r, iso) if k == "apex" else val(r, iso)
        e = sev(r, iso)
        if v is None or e is None:
            continue
        pts.append({"T": float(r["T_C"]), "y": v, "se": e, "good": is_good(r),
                    "withheld": is_withheld(r),
                    "src": key_of(r["dataset"]), "row": r})
    return pts

# raw -> absolute shift of the Sivan data, measured from the CSV.
# The shift is a frame definition, not a random error: Table 1 does not fold
# its scatter into the 2025 working-gas SE.  Both moments are reported.
SHIFT, SHIFT_SD, SHIFT_N = {}, {}, {}
for iso, rawc, corc in (("13", "D13_raw", "D13"), ("12", "D12_raw", "D12")):
    d = [num(r[corc]) - num(r[rawc]) for r in rows_of("siv")
         if num(r[rawc]) is not None and num(r[corc]) is not None]
    SHIFT[iso] = float(np.mean(d)) if d else 0.0
    SHIFT_N[iso] = len(d)
    SHIFT_SD[iso] = float(np.std(d, ddof=1)) if len(d) > 1 else 0.0

# ─────────────────────────────────────────────────────────────────────────────
# 3. THEORY CURVES
#    PIMC polynomials of Eldridge et al. 2019 (valid 270-800 K, T in K).
#    BMU = sum_{k=1}^{5} c_k / T^k fitted to their Table 4 (no constant term).
#    delta = mean of PIMC minus that fit on integer kelvin 283 through 423.
#    reference = PIMC at T >= 270 K, else BMU + delta.
#    Nothing in results_summary.json is allowed to replace this.
# ─────────────────────────────────────────────────────────────────────────────
def pimc13(T):
    T = np.asarray(T, dtype=float)
    return (1.47348e19/T**7 - 2.08648e17/T**6 + 1.19810e15/T**5 - 3.54757e12/T**4
            + 5.54476e9/T**3 - 3.49294e6/T**2 + 8.89370e2/T)

def pimc12(T):
    T = np.asarray(T, dtype=float)
    return (-9.67634e15/T**6 + 1.71917e14/T**5 - 1.24819e12/T**4
            + 4.30283e9/T**3 - 4.48660e6/T**2 + 1.86258e3/T)

BMU_T_C = np.array([-3.1, 1.2, 6.9, 16.9, 26.9, 36.9, 46.9, 50.5, 56.9, 66.9, 75.7,
                    76.9, 86.9, 96.9, 126.9, 127.8, 151.9, 165.4, 176.9, 201.9, 226.9,
                    251.9, 276.9, 301.9, 326.9, 351.9, 376.9, 401.9, 426.9, 451.9,
                    476.9, 501.9, 526.9])
BMU_13 = np.array([6.757, 6.589, 6.380, 6.032, 5.709, 5.409, 5.130, 5.033, 4.870, 4.628,
                   4.426, 4.401, 4.189, 3.989, 3.462, 3.446, 3.088, 2.907, 2.765, 2.482,
                   2.235, 2.016, 1.824, 1.652, 1.500, 1.364, 1.242, 1.133, 1.035, 0.947,
                   0.868, 0.796, 0.731])
BMU_12 = np.array([24.24, 23.42, 22.41, 20.73, 19.20, 17.80, 16.52, 16.07, 15.34, 14.25,
                   13.37, 13.26, 12.34, 11.50, 9.35, 9.29, 7.91, 7.24, 6.72, 5.74,
                   4.93, 4.24, 3.67, 3.19, 2.78, 2.43, 2.14, 1.88, 1.67, 1.48,
                   1.31, 1.17, 1.05])

def _bmu_design(TK):
    x = 1.0 / np.asarray(TK, dtype=float)
    return np.column_stack([x ** k for k in range(1, 6)])

_bmu_TK = BMU_T_C + 273.15
C_BMU13, *_ = np.linalg.lstsq(_bmu_design(_bmu_TK), BMU_13, rcond=None)
C_BMU12, *_ = np.linalg.lstsq(_bmu_design(_bmu_TK), BMU_12, rcond=None)

def _bmu_eval(coef, TK):
    raw = np.asarray(TK, dtype=float)
    flat = np.atleast_1d(raw)
    out = (_bmu_design(flat) @ coef).reshape(flat.shape)
    if raw.ndim == 0:
        return float(out[0])
    return out

def bmu13(TK): return _bmu_eval(C_BMU13, TK)
def bmu12(TK): return _bmu_eval(C_BMU12, TK)

_BRIDGE_K = np.arange(283.0, 424.0)          # 283..423 K inclusive
BD13 = float(np.mean(pimc13(_BRIDGE_K) - bmu13(_BRIDGE_K)))
BD12 = float(np.mean(pimc12(_BRIDGE_K) - bmu12(_BRIDGE_K)))
BMU_RMS = {
    "13": float(np.sqrt(np.mean((bmu13(_bmu_TK) - BMU_13) ** 2))),
    "12": float(np.sqrt(np.mean((bmu12(_bmu_TK) - BMU_12) ** 2))),
}

def _scalarize(TK, values):
    raw = np.asarray(TK, dtype=float)
    out = np.asarray(values, dtype=float)
    if raw.ndim == 0:
        return float(out.reshape(-1)[0])
    return out

def ref13(TK):
    TK = np.asarray(TK, dtype=float)
    return _scalarize(TK, np.where(TK >= 270.0, pimc13(TK), np.asarray(bmu13(TK), dtype=float) + BD13))

def ref12(TK):
    TK = np.asarray(TK, dtype=float)
    return _scalarize(TK, np.where(TK >= 270.0, pimc12(TK), np.asarray(bmu12(TK), dtype=float) + BD12))

REF = {"13": ref13, "12": ref12}     # reference curve (bridged BMU below 270 K)
ALT = {"13": pimc13, "12": pimc12}   # the other low-T curve: PIMC extrapolated
BMU = {"13": bmu13, "12": bmu12}
PIMC = ALT
BRIDGE = {"13": BD13, "12": BD12}

# The constant-delta bridge is NOT continuous at 270 K (PIMC - BMU is T-dependent,
# delta is a grid mean).  Report the step so residuals across the join are read correctly.
JOIN_STEP = {"13": float(pimc13(270.0) - (bmu13(270.0) + BD13)),
             "12": float(pimc12(270.0) - (bmu12(270.0) + BD12))}
print(f"[warn] reference-curve step at the 270 K join (above minus below): "
      f"Δ13 {JOIN_STEP['13']:+.4f}‰, Δ12 {JOIN_STEP['12']:+.4f}‰ (systematic; "
      "bracketed by the theory envelope in Table S6)")

# ─────────────────────────────────────────────────────────────────────────────
# 4. FITS (all recomputed from the CSV)
# ─────────────────────────────────────────────────────────────────────────────
def reff_for(variant, iso):
    return ALT[iso] if variant == "alt" else REF[iso]

def xy(q, reff, space):
    x = float(reff(q["T"] + 273.15))
    y = q["y"]
    if space == "ln":
        x = 1000.0 * math.log(1.0 + x / 1000.0)
        y = 1000.0 * math.log(1.0 + y / 1000.0)
    return x, y

def subset(pts, variant):
    """Withheld rows never enter a fit.  Bad rows were already removed at load,
    so every point here is Good: "all" = all Good temperatures, "valid" = T >= 270 K."""
    p = [q for q in pts if not q.get("withheld")]
    p = [q for q in p if q["good"]]
    if variant == "valid":
        p = [q for q in p if q["T"] + 273.15 >= 270.0]
    return p

def fit_m1(pts, reff):
    """Method 1: weighted offset (linear Delta space), SE inflated by sqrt(chi2/(n-1))."""
    sw = swr = 0.0; n = 0
    for q in pts:
        x = float(reff(q["T"] + 273.15))
        w = 1.0 / q["se"] ** 2
        sw += w; swr += w * (q["y"] - x); n += 1
    off = swr / sw
    chi2 = 0.0                      # scatter around the fitted offset, not around 0
    for q in pts:
        x = float(reff(q["T"] + 273.15))
        r = q["y"] - x - off
        chi2 += r * r / q["se"] ** 2
    chi2r = chi2 / (n - 1)
    se = math.sqrt(1.0 / sw) * math.sqrt(max(1.0, chi2r))
    return {"off": off, "se": se, "se_infl": se, "n": n, "chi2r": chi2r}

def fit_m2(pts, reff, space="lin"):
    """Method 2: weighted line y = m x + c, SEs inflated by sqrt(chi2r)."""
    sw = sx = sy = sxx = sxy = 0.0; n = len(pts)
    for q in pts:
        x, y = xy(q, reff, space)
        w = 1.0 / q["se"] ** 2
        sw += w; sx += w * x; sy += w * y; sxx += w * x * x; sxy += w * x * y
    D = sw * sxx - sx * sx
    m = (sw * sxy - sx * sy) / D
    c = (sxx * sy - sx * sxy) / D
    chi2 = 0.0
    for q in pts:
        x, y = xy(q, reff, space)
        chi2 += (y - m * x - c) ** 2 / q["se"] ** 2
    chi2r = chi2 / (n - 2)
    infl = math.sqrt(max(1.0, chi2r))
    se_m = math.sqrt(sw / D) * infl
    se_c = math.sqrt(sxx / D) * infl
    z = (m - 1.0) / se_m
    ci_half = t_critical_two_sided(n - 2) * se_m
    return {"m": m, "se_m": se_m, "c": c, "se_c": se_c, "z": z,
            "p": t_two_sided_p(abs(z), n - 2), "chi2r": chi2r, "n": n,
            "ci95_lo": m - ci_half, "ci95_hi": m + ci_half}

VARIANTS = ["all", "valid", "alt", "ln"]   # no "+Bad" variant: Bad rows are excluded at load
PTS = {(k, iso): get_pts(k, iso) for k in ALLKEYS for iso in ("13", "12")}
M1 = {}   # (k, iso, variant) -> m1 dict   (linear space)
M2 = {}   # (k, iso, variant) -> m2 dict
for k in ALLKEYS:
    for iso in ("13", "12"):
        for v in VARIANTS:
            pts = subset(PTS[(k, iso)], v)
            reff = reff_for(v, iso)
            space = "ln" if v == "ln" else "lin"
            M2[(k, iso, v)] = fit_m2(pts, reff, space)
            if v in ("all", "valid", "alt"):
                M1[(k, iso, v)] = fit_m1(pts, reff)

# offsets at T >= 270 K in each plotted data set's own frame
OFF = {k: {iso: (M1[(k, iso, "valid")]["off"], M1[(k, iso, "valid")]["se_infl"])
           for iso in ("13", "12")} for k in ALLKEYS}
# AP613 2025 working-gas offset via the Sivan raw->absolute shift.
# ALT uses the 2021 Method-1 working-gas offset as the shift instead
# (absolute = raw - o_2021).  For Delta13 the two agree.  For Delta12 the
# CSV shift is not equal to -o_2021, so both are reported and only SHIFT
# is used in the tables and figures.
OFF_WG = {iso: OFF["das"][iso][0] - SHIFT[iso] for iso in ("13", "12")}
SE_WG = {iso: OFF["das"][iso][1] for iso in ("13", "12")}
OFF_WG_ALT = {iso: OFF["das"][iso][0] + OFF["siv"][iso][0] for iso in ("13", "12")}
SE_WG_ALT = {iso: math.hypot(OFF["das"][iso][1], OFF["siv"][iso][1]) for iso in ("13", "12")}

# campaign comparison on the absolute scale: separate M1 offsets + ANCOVA
def ancova(pts, reff):
    X, Y, W = [], [], []
    for q in pts:
        x = float(reff(q["T"] + 273.15))
        X.append([x, 1.0 if q["src"] == "das" else 0.0, 1.0 if q["src"] == "siv" else 0.0])
        Y.append(q["y"]); W.append(1.0 / q["se"] ** 2)
    X = np.array(X); Y = np.array(Y); W = np.array(W)
    XtW = X.T * W
    A = XtW @ X; b = XtW @ Y
    beta = np.linalg.solve(A, b)
    cov = np.linalg.inv(A)
    resid = Y - X @ beta
    chi2 = float(np.sum(W * resid ** 2)); n = len(Y)
    infl = math.sqrt(max(1.0, chi2 / (n - 3)))
    cov = cov * infl * infl
    dc = beta[1] - beta[2]
    se_dc = math.sqrt(cov[1, 1] + cov[2, 2] - 2.0 * cov[1, 2])
    z_dc = dc / se_dc
    return {"m": float(beta[0]), "se_m": math.sqrt(cov[0, 0]),
            "dc": float(dc), "se_dc": se_dc, "z_dc": z_dc,
            "p_dc": t_two_sided_p(abs(z_dc), n - 3), "chi2r": chi2 / (n - 3), "n": n}

CAMP = {}
for iso in ("13", "12"):
    for v in ("all", "valid", "alt"):
        pts = subset(PTS[("apex", iso)], v)
        reff = reff_for(v, iso)
        das_pts = [q for q in pts if q["src"] == "das"]
        siv_pts = [q for q in pts if q["src"] == "siv"]
        od, os_ = fit_m1(das_pts, reff), fit_m1(siv_pts, reff)
        diff = od["off"] - os_["off"]
        se_diff = math.hypot(od["se_infl"], os_["se_infl"])
        CAMP[(iso, v)] = {
            "das": od, "siv": os_, "diff": diff, "se_diff": se_diff, "z": diff / se_diff,
            "ancova": ancova(pts, reff)}

def refit(k, iso, base, drop):
    """Method 1 and Method 2 on subset(base) after dropping points with drop(q)."""
    pts = [q for q in subset(PTS[(k, iso)], base) if not drop(q)]
    if len(pts) < 4:
        return None
    return fit_m1(pts, REF[iso]), fit_m2(pts, REF[iso], "lin")

def _is_550(q):
    return abs(q["T"] - 550.0) < 1e-6

def _is_m20(q):
    return abs(q["T"] + 20.0) < 1e-6

# Leave-one-out of the AP613 -20 C Good points on the all-Good extended fit.
# Either point alone still anchors 253 K, so dropping one barely moves the
# slope and dropping both releases it.  The two changes do not add.
LOO = {}
for iso in ("13", "12"):
    cold = [q for q in subset(PTS[("apex", iso)], "all") if _is_m20(q)]
    cold = sorted(cold, key=lambda q: q["y"])
    base = M2[("apex", iso, "all")]["m"]
    specs = [("drop neither", [])]
    if len(cold) >= 1:
        specs.append(("drop lower -20", [cold[0]]))
    if len(cold) >= 2:
        specs.append(("drop higher -20", [cold[-1]]))
    if cold:
        specs.append(("drop all -20", cold))
    LOO[iso] = []
    for label, gone in specs:
        gone_ids = {id(q) for q in gone}
        fit = refit("apex", iso, "all", lambda q, ids=gone_ids: id(q) in ids)
        if fit is None:
            continue
        LOO[iso].append({"label": label, "m": fit[1]["m"], "se": fit[1]["se_m"],
                         "n": fit[1]["n"], "p": fit[1]["p"], "delta_m": fit[1]["m"] - base})

# Hot-end leverage: AP613 550 C tubes have a different bulk composition.
NO550 = {}
for k in ("das", "apex"):
    for iso in ("13", "12"):
        fit = refit(k, iso, "valid", _is_550)
        if fit is None:
            continue
        base_m, base_o = M2[(k, iso, "valid")], M1[(k, iso, "valid")]
        NO550[(k, iso)] = {
            "m": fit[1]["m"], "se_m": fit[1]["se_m"], "n": fit[1]["n"], "p": fit[1]["p"],
            "off": fit[0]["off"], "se_off": fit[0]["se_infl"],
            "delta_m": fit[1]["m"] - base_m["m"],
            "delta_off": fit[0]["off"] - base_o["off"],
        }

# ─────────────────────────────────────────────────────────────────────────────
# 5. RESIDUALS AND EQ. (5) FRACTIONS (offset fitted at T >= 270 K, own frame)
# ─────────────────────────────────────────────────────────────────────────────
def resid(row, iso):
    v = val(row, iso)
    if v is None:
        return None
    off = OFF[key_of(row["dataset"])][iso][0]
    return v - float(REF[iso](float(row["T_C"]) + 273.15)) - off

def frac_covered(row, iso):
    """Eq. (5): fraction of the distance from the starting (cylinder) composition
    to the reference value covered, evaluated in the working-gas frame (start = 0).
    AP613 2025 is stored in the absolute frame, so it is converted with the 2021
    shift and the working-gas offset OFF_WG.  Undefined (None) for the Eldridge
    means and where the reference value lies at or below the starting
    composition.  (Bad rows no longer exist here: they were removed at load.)"""
    v = val(row, iso)
    if v is None:
        return None
    k = key_of(row["dataset"])
    if k == "eld":
        return None

    if k == "das":
        v = v - SHIFT[iso]
        off = OFF_WG[iso]
    else:
        off = OFF[k][iso][0]
    den = float(REF[iso](float(row["T_C"]) + 273.15)) + off
    if den <= 1e-9:
        return None
    return v / den

def frac_se(row, iso):
    """SE of Eq. (5) from the measurement and the T >= 270 K offset.

    The 2021 shift is treated as a fixed frame definition (its SD is reported
    separately).  None if the fraction is undefined or the point itself is
    inside the offset fit, in which case the offset is not an independent error.
    """
    fr = frac_covered(row, iso)
    if fr is None or in_offset_fit(row):
        return None
    e = sev(row, iso)
    if not e:
        return None
    k = key_of(row["dataset"])
    if k == "das":
        v = val(row, iso) - SHIFT[iso]
        off, se_off = OFF_WG[iso], SE_WG[iso]
    else:
        v = val(row, iso)
        off, se_off = OFF[k][iso]
    den = float(REF[iso](float(row["T_C"]) + 273.15)) + off
    return math.hypot(e / den, v * se_off / den ** 2)

def z_with_offset(row, iso):
    """Residual / hypot(measurement SE, offset SE).  Measurement-only z
    overstates a cold residual.  None when the point entered the offset fit."""
    rr, e = resid(row, iso), sev(row, iso)
    if rr is None or not e or in_offset_fit(row):
        return None
    return rr / math.hypot(e, OFF[key_of(row["dataset"])][iso][1])

def dDdT(iso, TK, curve=None):
    curve = curve or REF[iso]
    eps = 0.25
    return float(curve(TK + eps) - curve(TK - eps)) / (2.0 * eps)

def invert_curve(delta, iso, curve=None):
    """Temperature (C) at which curve(T) = delta.  The curves decrease with T."""
    curve = curve or REF[iso]
    lo, hi = 200.0, 1000.0
    flo, fhi = float(curve(lo)), float(curve(hi))
    if not (fhi <= delta <= flo):
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if float(curve(mid)) > delta:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) - 273.15

def theory_envelope(iso):
    """253.15 K values under the reference bridge, a value-match at 270 K,
    the extrapolated PIMC polynomial, and the grid-mean bridge on two
    other windows.  This is the theory range; it is systematic."""
    fit = BMU[iso]
    pimc = ALT[iso]
    windows = ((283.0, 423.0), (270.0, 400.0), (300.0, 450.0))
    values = {}
    for lo, hi in windows:
        grid = np.arange(lo, hi + 1.0)
        delta = float(np.mean(pimc(grid) - fit(grid)))
        values[f"grid {lo:.0f}-{hi:.0f} K"] = float(fit(253.15) + delta)
    values["match PIMC at 270 K"] = float(fit(253.15) + (pimc(270.0) - fit(270.0)))
    values["PIMC polynomial extrapolated"] = float(pimc(253.15))
    vals = list(values.values())
    return values, min(vals), max(vals)

# ─────────────────────────────────────────────────────────────────────────────
# 6. TABLE S8: MONTE CARLO INVERSION OF THE REFERENCE CURVE (68 % interval)
# ─────────────────────────────────────────────────────────────────────────────
T_GRID = [float(t) for t in args.t_grid.split(",") if t.strip()]
# Wide inversion grid: np.interp clamps draws outside it, which would bias the
# 68% intervals at the edges (worst for Delta12 at high T, where the curve flattens).
_TK_GRID = np.linspace(200.0, 1000.0, 32001)


def t_confidence(off, se_off, iso, prec, seed=1):
    """68 % (-/+) temperature uncertainty at each true T in T_GRID.
    Draw Delta ~ N(ref(T), sqrt(prec^2 + se_off^2)) and invert the curve."""
    rng = np.random.default_rng(seed)
    d_grid = REF[iso](_TK_GRID)          # decreasing in T
    out = {}
    for T in T_GRID:
        d = rng.normal(float(REF[iso](T + 273.15)), math.hypot(prec, se_off), args.mc_draws)
        t_est = np.interp(d, d_grid[::-1], _TK_GRID[::-1]) - 273.15
        lo, hi = np.percentile(t_est, [15.865, 84.135])
        out[T] = {"minus": T - lo, "plus": hi - T}
    return out

T_CONF = {}   # (lab, iso, variant) -> {T: {minus, plus}}
for lab in ("apex", "das", "cal"):
    for iso in ("13", "12"):
        off, se = OFF[lab][iso]
        T_CONF[(lab, iso, "repro")] = t_confidence(off, se, iso, args.repro13 if iso == "13" else args.repro12)
        T_CONF[(lab, iso, "internal")] = t_confidence(off, se, iso, args.internal13 if iso == "13" else args.internal12, seed=2)

# ─────────────────────────────────────────────────────────────────────────────
# 7. FIGURES (matplotlib only - no HTML, no Plotly)
# ─────────────────────────────────────────────────────────────────────────────
T_GRID_K = np.linspace(245.0, 830.0, 1200)
T_GRID_C = T_GRID_K - 273.15
LOWT_C = 270.0 - 273.15
XMIN, XMAX = -30.0, 560.0

def save(fig, base):
    fig.savefig(os.path.join(OUT, base + ".pdf"))
    fig.savefig(os.path.join(OUT, base + ".png"), dpi=args.dpi)
    plt.close(fig)
    print(f"[fig ] {base}.pdf / .png")

def draw_theory(ax, iso):
    ax.axvspan(XMIN, LOWT_C, color="grey", alpha=0.18, zorder=0)
    hi = T_GRID_K >= 270.0
    ax.plot(T_GRID_C[hi], PIMC[iso](T_GRID_K[hi]), color="k", lw=1.6, ls="-",
            label="PIMC (Eldridge 2019)")
    ax.plot(T_GRID_C[~hi], PIMC[iso](T_GRID_K[~hi]), color="k", lw=1.4, ls=":")
    ax.plot(T_GRID_C, BMU[iso](T_GRID_K), color="0.45", lw=1.4, ls="--",
            label="BMU, Table 4 fit")
    ax.plot(T_GRID_C, REF[iso](T_GRID_K), color="k", lw=2.6, label="reference curve")

def marker_face(row, k):
    """Filled = in the T >= 270 K fit.  Open = withheld or T < 270 K (Bad rows are gone)."""
    use = is_good(row) and not is_withheld(row) and float(row["T_C"]) + 273.15 >= 270.0
    return COLOR[k] if use else "none"

def draw_points(ax, iso):
    for k in KEYS:
        for r in rows_of(k):
            v = val(r, iso)
            if v is None:
                continue
            if k == "das":
                v -= SHIFT[iso]
            edge = "k" if is_withheld(r) else COLOR[k]
            ax.scatter([float(r["T_C"])], [v], marker=MARK[k], s=42,
                       facecolors=marker_face(r, k),
                       edgecolors=edge, linewidths=1.4, zorder=3)

def fig1():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharex=True)
    for ax, iso in zip(axes, ("13", "12")):
        draw_theory(ax, iso)
        for k in KEYS:
            ts = [float(r["T_C"]) for r in rows_of(k)]
            tg = np.linspace(min(ts) + 273.15, max(ts) + 273.15, 300)
            off, se = (OFF_WG[iso], SE_WG[iso]) if k == "das" else OFF[k][iso]
            ax.plot(tg - 273.15, REF[iso](tg) + off, color=COLOR[k], lw=1.8,
                    label=f"{NAME[k]} (offset {off:+.2f}‰)")
            ax.fill_between(tg - 273.15, REF[iso](tg) + off - se, REF[iso](tg) + off + se,
                            color=COLOR[k], alpha=0.15, lw=0)
        draw_points(ax, iso)
        ax.set_xlabel("equilibration temperature (°C)")
        ax.set_ylabel(f"{ISO_LAB[iso]} (‰)")
        ax.set_xlim(XMIN, XMAX)
        ax.axvline(LOWT_C, color="grey", lw=0.8, ls="--")
        h, l = ax.get_legend_handles_labels()
        ax.legend(h, l, fontsize=7.5, loc="upper right", framealpha=0.9)
    fig.suptitle("Figure 1.  Working-gas frame (AP613 2025 minus the 2021 shift).  "
                 "Curves: reference + Method-1 offset at T ≥ 270 K.  "
                 "Filled symbols are in that fit; open symbols are not "
                 "(T < 270 K or CAL1549 at −20 °C withheld; Bad rows are excluded at load).", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "Fig1_Delta_vs_T")

def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharex=True)
    for ax, iso in zip(axes, ("13", "12")):
        ax.axvspan(XMIN, LOWT_C, color="grey", alpha=0.18, zorder=0)
        ax.axhline(0.0, color="k", lw=1.2)
        for k in KEYS:
            for r in rows_of(k):
                v, e = val(r, iso), sev(r, iso)
                if v is None or e is None:
                    continue
                rr = resid(r, iso)
                edge = "k" if is_withheld(r) else COLOR[k]
                ax.errorbar([float(r["T_C"])], [rr], yerr=[e], fmt="none",
                            ecolor=COLOR[k], elinewidth=1.0, alpha=0.8, zorder=2)
                ax.scatter([float(r["T_C"])], [rr], marker=MARK[k], s=42,
                           facecolors=marker_face(r, k),
                           edgecolors=edge, linewidths=1.4, zorder=3)
        ax.set_xlabel("equilibration temperature (°C)")
        ax.set_ylabel(f"residual {ISO_LAB[iso]} (‰)")
        ax.set_xlim(XMIN, XMAX)
    fig.suptitle("Figure 2.  Residual from the reference curve + the T ≥ 270 K offset.  "
                 "Error bars are the measurement SE only (Table 2 also gives the "
                 "offset-aware z).  Filled: in the offset fit.  Open: not in the fit.  "
                 "Grey: T < 270 K.", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "Fig2_residuals_vs_T")

def figS3():
    """Low-temperature results in the absolute frame, with gas starts and reference."""
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.0))
    ap_rows = [r for r in ROWS if key_of(r["dataset"]) == "das" and float(r["T_C"]) == -20.0]
    cal_rows = [r for r in ROWS if key_of(r["dataset"]) == "cal" and float(r["T_C"]) == -20.0]
    plotted_values = {}

    for ax, iso in zip(axes, ("13", "12")):
        ref = float(REF[iso](253.15))
        starts = {"AP613": SHIFT[iso], "CAL1549": -OFF["cal"][iso][0]}
        plotted_values[iso] = {"reference": ref, "AP613 start": starts["AP613"],
                               "CAL1549 start": starts["CAL1549"]}
        ax.axhline(ref, color="black", lw=1.8, label=f"reference {ref:.2f}‰")
        ax.hlines(starts["AP613"], 0.3, 2.0, color=COLOR["das"], ls=":", lw=1.5,
                  label=f"AP613 start {starts['AP613']:.2f}‰")
        ax.hlines(starts["CAL1549"], 2.55, 3.12, color=COLOR["cal"], ls=":", lw=1.5,
                  label=f"CAL1549 start {starts['CAL1549']:.2f}‰")

        for i, r in enumerate(ap_rows):   # Good -20 C AP613 rows only (Bad dropped at load)
            x = 0.65 + 0.89 * i
            value = val_abs(r, iso)
            error = sev(r, iso)
            if value is None or error is None:
                continue
            ax.errorbar(x, value, yerr=error, fmt=MARK["das"], color=COLOR["das"],
                        mfc=COLOR["das"], mec=COLOR["das"], capsize=3, zorder=3)

        for r in cal_rows:
            value = val(r, iso)
            error = sev(r, iso)
            if value is None or error is None:
                continue
            estimate = value - OFF["cal"][iso][0]
            ax.errorbar(2.8, estimate, yerr=error, fmt=MARK["cal"], color=COLOR["cal"],
                        mfc=COLOR["cal"], capsize=3, zorder=3,
                        label="CAL1549: working-gas value minus the fitted offset (absolute frame)")
            raw = num(r["D13_raw"] if iso == "13" else r["D12_raw"])
            plotted_values[iso]["CAL1549 fitted"] = estimate
            if raw is not None:
                ax.plot(2.94, raw, marker=MARK["cal"], color=COLOR["cal"],
                        mfc="none", ms=7, zorder=4, label="CAL1549: raw-frame value")
                plotted_values[iso]["CAL1549 raw"] = raw

        ap_x = [0.65 + 0.89 * i for i in range(len(ap_rows))]
        ax.set_xticks(ap_x + [2.8],
                      [f"AP613\n({r['duration']})" for r in ap_rows] + ["CAL1549\n(5 months)"])
        ax.set_xlim(0.25, 3.15)
        ax.set_ylabel(f"{ISO_LAB[iso]} (‰, absolute frame)")
        ax.grid(axis="y", alpha=0.2)
    handles, _ = axes[0].get_legend_handles_labels()
    legend_keys = ("reference", "AP613 start", "CAL1549 start", "CAL1549 fitted", "CAL1549 raw")
    legend_labels = [
        f"{key}: Δ¹³ {plotted_values['13'][key]:.2f}‰; "
        f"Δ¹² {plotted_values['12'][key]:.2f}‰"
        for key in legend_keys
    ]
    fig.legend(handles, legend_labels, fontsize=7.5, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, 0.01), frameon=False)
    axes[0].set_title("Δ¹³CH₃D")
    axes[1].set_title("Δ¹²CH₂D₂")
    fig.suptitle("Figure S3.  −20 °C in the absolute frame: a kinetic constraint, not a calibration point",
                 fontsize=11)
    fig.text(0.5, 0.93, "Bars are measurement SE.  Only Good rows are shown (Bad rows were dropped at load).  "
             "The two AP613 points share one duration string and may be replicates of one tube.  "
             "CAL1549 is withheld from every fit; the open symbol is the raw column, "
             "the filled symbol is working-gas value minus the fitted offset.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=[0, 0.12, 1, 0.90])
    save(fig, "FigS3_low_temperature")

def figS1():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))
    for ax, iso in zip(axes, ("13", "12")):
        for k in KEYS:
            y = M2[(k, iso, "valid")]
            y_all = M2[(k, iso, "all")]
            xs = []
            for r in rows_of(k):
                v = val(r, iso)
                if v is None:
                    continue
                x = float(REF[iso](float(r["T_C"]) + 273.15))
                xs.append(x)
                low = float(r["T_C"]) < LOWT_C
                ax.scatter([x], [v], marker=MARK[k], s=70 if low else 42,
                           facecolors=marker_face(r, k),
                           edgecolors="k" if low or is_withheld(r) else COLOR[k],
                           linewidths=1.5 if low else 1.2, zorder=4 if low else 3)
            xs = np.array(xs)
            xg = np.linspace(xs.min() - 0.2, xs.max() + 0.2, 100)
            ax.plot(xg, y["m"] * xg + y["c"], color=COLOR[k], lw=1.7,
                    label=(f"{NAME[k]} T≥270 K: m={y['m']:.3f}±{y['se_m']:.3f}"))
            if abs(y_all["m"] - y["m"]) > 0.01:
                ax.plot(xg, y_all["m"] * xg + y_all["c"], color=COLOR[k], lw=1.0, ls="--",
                        label=f"{NAME[k]} all Good: m={y_all['m']:.3f} (sensitivity)")
            band = np.sqrt((y["se_m"] * xg) ** 2 + y["se_c"] ** 2)
            ax.fill_between(xg, y["m"] * xg + y["c"] - band, y["m"] * xg + y["c"] + band,
                            color=COLOR[k], alpha=0.12, lw=0)
        ax.set_xlabel(f"reference-curve {ISO_LAB[iso]} (‰)")
        ax.set_ylabel(f"measured {ISO_LAB[iso]} (‰, own frame)")
        ax.legend(fontsize=6.6, loc="best", framealpha=0.9)
    fig.suptitle("Figure S1.  Method 2 in each series' own frame.  Solid lines and bands are the "
                 "T ≥ 270 K fit (primary).  Dashed lines, where drawn, include Good points below 270 K.  "
                 "Open symbols are not in the primary fit; black edges mark T < 270 K.", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save(fig, "FigS1_method2_regression")

def figS2():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharex=True)
    for ax, iso in zip(axes, ("13", "12")):
        draw_theory(ax, iso)
        for k in KEYS:
            off, se = OFF[k][iso]
            ax.fill_between(T_GRID_C, REF[iso](T_GRID_K) - se, REF[iso](T_GRID_K) + se,
                            color=COLOR[k], alpha=0.12, lw=0)
            for r in rows_of(k):
                v = val(r, iso)
                if v is None:
                    continue
                face = marker_face(r, k)
                ax.scatter([float(r["T_C"])], [v - off], marker=MARK[k], s=42,
                           facecolors=face,
                           edgecolors="k" if is_withheld(r) else COLOR[k], linewidths=1.4,
                           alpha=0.45 if face == "none" else 1.0, zorder=3)
        ax.set_xlabel("equilibration temperature (°C)")
        ax.set_ylabel(f"{ISO_LAB[iso]} shifted by own offset (‰)")
        ax.set_xlim(XMIN, XMAX)
        ax.axvline(LOWT_C, color="grey", lw=0.8, ls="--")
    fig.suptitle("Figure S2.  Each series minus its own Method-1 offset (T ≥ 270 K), own frame.  "
                 "Filled symbols are in that fit.  Faded open symbols are not.  "
                 "Bands: ±1σ of the offset.  Grey: T < 270 K.", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "FigS2_shifted_by_offset")

# ─────────────────────────────────────────────────────────────────────────────
# 8. TABLES (collected into one Excel workbook)
# ─────────────────────────────────────────────────────────────────────────────
TABLES = {}

def write_table(name, header, rows, footnotes=()):
    TABLES[name] = (header, rows, footnotes)

def excel_value(value):
    if not isinstance(value, str):
        return value
    value = value.strip()
    if not value:
        return None
    if re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", value):
        return float(value)
    return value

def write_workbook(manifest_rows):
    workbook = Workbook()
    manifest = workbook.active
    manifest.title = "Manifest"
    sheets = {"Manifest": (["Item", "Value"], manifest_rows, ())}
    table_order = sorted(TABLES, key=lambda name: (0 if name in ("Table1", "Table2") else 1,
                                                    int(name.removeprefix("TableS")) if name.startswith("TableS") else 0))
    sheets.update({name: TABLES[name] for name in table_order})
    for sheet_name, (header, rows, footnotes) in sheets.items():
        worksheet = manifest if sheet_name == "Manifest" else workbook.create_sheet(sheet_name)
        worksheet.append(header)
        for row in rows:
            worksheet.append([excel_value(value) for value in row])
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = f"A1:{worksheet.cell(row=1, column=len(header)).column_letter}{len(rows) + 1}"
        if footnotes:
            worksheet.append([])
            for note in footnotes:
                worksheet.append([note])
                row_number = worksheet.max_row
                worksheet.merge_cells(start_row=row_number, start_column=1,
                                      end_row=row_number, end_column=len(header))
                cell = worksheet.cell(row=row_number, column=1)
                cell.font = Font(italic=True, size=9)
                cell.alignment = Alignment(wrap_text=True, vertical="top")
                worksheet.row_dimensions[row_number].height = 30
        for cell in worksheet[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4E78")
        for column in worksheet.columns:
            width = min(max(len(str(cell.value or "")) for cell in column) + 2, 42)
            worksheet.column_dimensions[column[0].column_letter].width = width

    path = os.path.join(OUT, "AMT_tables.xlsx")
    workbook.save(path)
    return path

def f2(x): return "" if x is None else f"{x:.2f}"
def f3(x): return "" if x is None else f"{x:.3f}"
def f4(x): return "" if x is None else f"{x:.4f}"

def adjusted_p(p, primary):
    if not primary or p is None:
        return None
    return min(1.0, N_PRIMARY * p)

_MONTH_DAYS = 365.2425 / 12.0
_DAY_UNITS = {"hour": 1.0 / 24.0, "day": 1.0, "week": 7.0, "month": _MONTH_DAYS}

def duration_to_days(text):
    """Approximate days.  Ranges become midpoints.  Bare numbers are not guessed."""
    raw = (text or "").strip()
    if raw.lower() in ("", "none", "nan"):
        return None, "not recorded"
    s = (raw.replace("−", "-").replace("–", "-").replace("~", "")
            .replace("≈", "").strip())
    if re.fullmatch(r"[A-Za-z]+", s):
        return None, "unit without a number"
    m = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*-\s*([0-9]*\.?[0-9]+)\s*([A-Za-z]+)?", s)
    if m:
        a, b, unit = float(m.group(1)), float(m.group(2)), m.group(3)
        if not unit:
            return None, "range without a unit; left unconverted"
        days, note = _scale_duration(0.5 * (a + b), unit)
        return days, "midpoint of the recorded range"
    m = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*([A-Za-z]+)", s)
    if not m:
        return None, "unparsed"
    return _scale_duration(float(m.group(1)), m.group(2))

def _scale_duration(number, unit):
    key = unit.lower().rstrip("s")
    if key not in _DAY_UNITS:
        return None, f"unknown unit ({unit})"
    return number * _DAY_UNITS[key], "converted; months use 365.2425/12 days"

def table1():
    rows = []
    role_of = {"siv": "primary", "das": "primary", "cal": "primary",
               "apex": "consistency", "eld": "consistency"}
    for k in ("siv", "das", "cal", "apex", "eld"):
        role = role_of[k]
        for iso in ("13", "12"):
            v = M2[(k, iso, "valid")]
            off, se = (OFF_WG[iso], SE_WG[iso]) if k == "das" else OFF[k][iso]
            p_adj = adjusted_p(v["p"], role == "primary")
            rows.append([NAME[k], role, ISO_LAB[iso], v["n"], f2(off), f2(se),
                         f3(v["m"]), f3(v["se_m"]), f3(v["ci95_lo"]), f3(v["ci95_hi"]),
                         f"{v['z']:+.2f}", f4(v["p"]),
                         "" if p_adj is None else f4(p_adj)])
    write_table("Table1",
                ["data set", "role", "isotopologue", "n",
                 "offset (‰)", "SE offset (‰)", "slope m", "SE(m)",
                 "95% CI low", "95% CI high", "t", "p", "Bonferroni 6p"], rows,
                footnotes=[
                    "Offset is Method 1 (slope fixed at 1) at T >= 270 K.  Slope is Method 2 on the same points.  The two are not the intercept and slope of one fit.  SE is inflated by sqrt(chi2r) when chi2r > 1.",
                    "Primary family, 6 tests: AP613 2021, AP613 2025 and CAL1549, both isotopologues, T >= 270 K.  Adjusted p = min(1, 6p).  A CI that contains 1 is not a demonstration that the slope is 1.",
                    "AP613 2025 offset is the absolute-frame Method-1 offset minus the 2021 raw-to-absolute shift measured on the CSV.  The SE is the absolute-fit SE.  It does not include the shift scatter (see the manifest).  For Delta12 that shift is not equal to minus the 2021 Method-1 offset; the alternative is in the manifest and is not used here.",
                    "Unmarked primary offsets are in the working-gas frame.  AP613 extended and Eldridge are in the absolute frame and are consistency checks, not members of the slope family.  Extended nests the 2021 anchor.  Eldridge's values were assigned from this PIMC calibration, and Table S9 refits published means rather than their individual analyses.",
                    "Points below 270 K and the withheld CAL1549 -20 C row are not in this table.  Bad rows are excluded at load and appear in no table.",
                ])

def table2_and_S6():
    header = ["data set", "T (°C)", "catalyst", "duration", "flag", "in offset fit", "withheld",
              "Δ13 (‰)", "SE13", "resid13 (‰)", "z13 meas", "z13 with offset", "frac13", "SE frac13",
              "Δ12 (‰)", "SE12", "resid12 (‰)", "z12 meas", "z12 with offset", "frac12", "SE frac12"]
    rows = []
    for r in sorted(ROWS, key=lambda r: (r["dataset"], float(r["T_C"]))):
        if float(r["T_C"]) > 25.0:
            continue
        out = [ds_label(r), r["T_C"], r["catalyst"], r["duration"], r["flag"],
               "yes" if in_offset_fit(r) else "no",
               "yes" if is_withheld(r) else "no"]
        for iso in ("13", "12"):
            v, e = val(r, iso), sev(r, iso)
            rr, fr = resid(r, iso), frac_covered(r, iso)
            zo, fs = z_with_offset(r, iso), frac_se(r, iso)
            out += [f2(v), f2(e), f2(rr),
                    ("" if rr is None or not e else f"{rr / e:+.1f}"),
                    ("" if zo is None else f"{zo:+.1f}"),
                    ("" if fr is None else f"{fr:.3f}"),
                    ("" if fs is None else f"{fs:.3f}")]
        rows.append(out)
    write_table("Table2", header, rows, footnotes=[
        "z meas = residual / measurement SE.  It ignores the offset uncertainty and is too small a denominator for the -20 C deficit.",
        "z with offset = residual / hypot(measurement SE, T>=270 K offset SE).  Blank when the point is inside that offset fit, because the offset is then not an independent error.",
        "frac is Eq. (5) in the working-gas frame, where the cylinder is 0.  AP613 2025 is converted with the 2021 shift.  The fraction SE includes the measurement and the offset, not the theory error and not the shift scatter.",
        "Blank fractions: Eldridge means have no cylinder zero in this table; a non-positive distance from the cylinder to the reference is undefined.  Bad rows are excluded at load and do not appear in this table.",
        "CAL1549 at -20 C is withheld from every fit.  Do not read its Delta13 residual as agreement with the curve.",
    ])

def tableS1():
    rows = []
    for r in ROWS:
        days, note = duration_to_days(r["duration"])
        rows.append([ds_label(r), r["T_C"], r["catalyst"], r["duration"],
                     "" if days is None else f"{days:.2f}", note,
                     r["flag"], "yes" if is_withheld(r) else "no"])
    write_table("TableS1",
                ["data set", "T (°C)", "catalyst", "duration as recorded",
                 "duration (days)", "duration note", "flag", "withheld from fits"],
                rows,
                footnotes=[
                    "Days are a convenience for comparing the cold experiments.  A range is its midpoint.  One month = 365.2425/12 days.  Eldridge durations have no unit in the CSV and are not converted.",
                    "The setpoint used for -20 C rows is -20 C (253.15 K), not -18 C.",
                ])

def tableS2():
    rows = []
    for k in ("siv", "apex", "cal", "eld"):
        for iso in ("13", "12"):
            a, v = M1[(k, iso, "all")], M1[(k, iso, "valid")]
            rows.append([NAME[k], FRAME[k], ISO_LAB[iso],
                         f2(a["off"]), f2(a["se_infl"]), a["n"], f2(a["chi2r"]),
                         f2(v["off"]), f2(v["se_infl"]), v["n"], f2(v["chi2r"])])
    write_table("TableS2",
                ["data set", "frame", "isotopologue",
                 "offset all (‰)", "SE infl all", "n all", "chi2r all",
                 "offset T>=270K (‰)", "SE infl T>=270K", "n T>=270K", "chi2r T>=270K"], rows)

def tableS3():
    rows = []
    for k in ("siv", "apex", "cal", "eld"):
        for iso in ("13", "12"):
            y = M2[(k, iso, "all")]
            rows.append([NAME[k], ISO_LAB[iso], y["n"], f3(y["m"]), f3(y["se_m"]),
                         f3(y["ci95_lo"]), f3(y["ci95_hi"]),
                         f"{y['z']:+.2f}", f4(y["p"]),
                         f2(y["c"]), f2(y["se_c"]), f2(y["chi2r"])])
    write_table("TableS3",
                ["data set", "isotopologue", "n", "m", "SE(m)", "95% CI low", "95% CI high",
                 "t", "p",
                 "intercept c (‰)", "SE(c)", "chi2r"], rows,
                footnotes=[
                    "Sensitivity only: all Good points, including Good points below 270 K that have not been shown to be at equilibrium.  These eight rows are not the confirmatory family and are not Bonferroni-adjusted.",
                    "The p value that used to be quoted as adjusted 0.20 is the AP613-extended Delta13 row here (unadjusted p near 0.025) multiplied by 8.  That row is moved by the -20 C points.  Do not cite it as the T >= 270 K result.",
                    "CAL1549 at -20 C is withheld, so it is not in this table.",
                ])

def tableS4():
    lab = {"all": "all Good", "valid": "T>=270K", "alt": "alt low-T curve",
           "ln": "ln space"}
    rows = []
    for k in ("siv", "apex", "cal", "eld"):
        for iso in ("13", "12"):
            for v in VARIANTS:
                y = M2[(k, iso, v)]
                rows.append([NAME[k], ISO_LAB[iso], lab[v], f3(y["m"]), f3(y["se_m"]),
                             f3(y["ci95_lo"]), f3(y["ci95_hi"]), y["n"], f4(y["p"])])
    write_table("TableS4",
                ["data set", "isotopologue", "choice", "m", "SE(m)", "95% CI low",
                 "95% CI high", "n", "p"], rows,
                footnotes=[
                    "Four sensitivities, not confirmatory tests.  'all Good' includes Good points below 270 K.  'alt' extrapolates the PIMC polynomial below 270 K.  'ln' regresses 1000 ln(1+Delta/1000).  There is no '+Bad' variant: Bad rows are excluded at load.",
                    "CAL1549 at -20 C is withheld from every choice.  The primary numbers are the T>=270K rows of AP613 2021, AP613 2025 and CAL1549 in Table 1; AP613 2025 is not repeated in this table because the extended fit already contains it.",
                    "Eldridge slopes in this table refit published means.  They do not reproduce the published Delta12 slope of 0.98 ± 0.05 (Table S9) and they are not a test of PIMC.",
                ])

def _pool_slopes(keys, iso):
    fits = [M2[(k, iso, "valid")] for k in keys]
    weights = [1.0 / fit["se_m"] ** 2 for fit in fits]
    total = sum(weights)
    pooled_m = sum(w * fit["m"] for w, fit in zip(weights, fits)) / total
    pooled_se = math.sqrt(1.0 / total)
    z = (pooled_m - 1.0) / pooled_se
    zcrit = inv_norm(0.975)
    eld_w = weights[keys.index("eld")] / total if "eld" in keys else 0.0
    return pooled_m, pooled_se, z, pooled_m - zcrit * pooled_se, pooled_m + zcrit * pooled_se, eld_w

def tableS10():
    specs = [
        ("primary: AP613 2021, AP613 2025, CAL1549", ("siv", "das", "cal")),
        ("sensitivity: primary plus Eldridge means", ("siv", "das", "cal", "eld")),
    ]
    rows = []
    for label, keys in specs:
        for iso in ("13", "12"):
            m, se, z, lo, hi, eld_w = _pool_slopes(keys, iso)
            rows.append([label, ISO_LAB[iso], f3(m), f3(se), f3(lo), f3(hi),
                         f"{z:+.2f}", "" if "eld" not in keys else f"{100.0 * eld_w:.1f}"])
    write_table("TableS10",
                ["analysis", "isotopologue", "pooled slope", "SE", "95% CI low",
                 "95% CI high", "z vs 1", "Eldridge weight (%)"], rows,
                footnotes=[
                    "Inverse-variance pool of the Method-2 slopes at T >= 270 K.  The primary pool does not include Eldridge.  Adding those means gives them a large weight because the published standard errors are small, and the Delta12 mean-fit does not match the published individual-analysis slope.",
                    "The normal 95% interval is a summary of the slope estimates.  It is not a substitute for the six-test family in Table 1.",
                ])

def ivw(values, ses):
    weights = [1.0 / se ** 2 for se in ses]
    total = sum(weights)
    mean = sum(w * v for w, v in zip(weights, values)) / total
    return mean, math.sqrt(1.0 / total)

def low_temperature_summary():
    ap_rows = [r for r in ROWS if key_of(r["dataset"]) == "das"
               and abs(float(r["T_C"]) + 20.0) < 1e-9 and is_good(r)
               and val_abs(r, "13") is not None]
    cal_row = next(r for r in ROWS if key_of(r["dataset"]) == "cal"
                   and abs(float(r["T_C"]) + 20.0) < 1e-9)
    ref13 = float(REF["13"](253.15))
    ref12 = float(REF["12"](253.15))
    ys = [val_abs(r, "13") for r in ap_rows]
    ses = [sev(r, "13") for r in ap_rows]
    anchor_se = CAMP[("13", "valid")]["siv"]["se_infl"]
    mean, se_mean = ivw(ys, ses)
    i_hi = max(range(len(ys)), key=lambda i: ys[i])
    bound_hi = ys[i_hi] - 2.0 * math.hypot(ses[i_hi], anchor_se)
    # one-sided 95% on the mean, if the two rows are replicates of one tube
    bound_mean = mean - inv_norm(0.95) * math.hypot(se_mean, anchor_se)
    recommended = min(bound_hi, bound_mean)
    env13, env13_lo, env13_hi = theory_envelope("13")
    env12, env12_lo, env12_hi = theory_envelope("12")
    deficits = {}
    for label, key in (("AP613 2025 offset", "das"), ("AP613 extended offset", "apex")):
        for iso, ref in (("13", ref13), ("12", ref12)):
            offset, offset_se = OFF[key][iso]
            rows = [r for r in ap_rows if val_abs(r, iso) is not None and sev(r, iso)]
            residuals = [val_abs(r, iso) - ref - offset for r in rows]
            rse = [sev(r, iso) for r in rows]
            m, se = ivw(residuals, rse)
            se = math.hypot(se, offset_se)
            deficits[(label, iso)] = (m, se, m / se)
    return {
        "ap_rows": ap_rows,
        "cal_row": cal_row,
        "ref13": ref13,
        "ref12": ref12,
        "durations": sorted({r["duration"] for r in ap_rows}),
        "anchor_se": anchor_se,
        "higher": ys[i_hi],
        "higher_se": ses[i_hi],
        "mean": mean,
        "se_mean": se_mean,
        "bound_hi": bound_hi,
        "bound_mean": bound_mean,
        "recommended_bound": recommended,
        "T_bound": invert_curve(recommended, "13"),
        "T_pimc_on_ref": invert_curve(float(ALT["13"](253.15)), "13"),
        "T_env_lo": invert_curve(env13_hi, "13"),
        "T_env_hi": invert_curve(env13_lo, "13"),
        "slope13": dDdT("13", 253.15),
        "slope12": dDdT("12", 253.15),
        "env13": env13, "env13_lo": env13_lo, "env13_hi": env13_hi,
        "env12": env12, "env12_lo": env12_lo, "env12_hi": env12_hi,
        "pimc13": float(ALT["13"](253.15)),
        "pimc12": float(ALT["12"](253.15)),
        "deficits": deficits,
        "cal_fit_absolute": val(cal_row, "13") - OFF["cal"]["13"][0],
        "cal_raw": num(cal_row["D13_raw"]),
        "cal_fit_absolute12": val(cal_row, "12") - OFF["cal"]["12"][0],
        "cal_raw12": num(cal_row["D12_raw"]),
        "cal_fraction": frac_covered(cal_row, "13"),
        "cal_fraction_se": frac_se(cal_row, "13"),
        "cal_fraction12": frac_covered(cal_row, "12"),
        "cal_fraction12_se": frac_se(cal_row, "12"),
    }

def tableS6(low):
    """Low-temperature inferences.  Not a second copy of Table 2."""
    d13 = low["deficits"][("AP613 2025 offset", "13")]
    d13x = low["deficits"][("AP613 extended offset", "13")]
    d12 = low["deficits"][("AP613 2025 offset", "12")]
    rows = [
        ["reference at -20 C", f"{low['ref13']:.3f}", f"{low['ref12']:.3f}",
         "PIMC at and above 270 K; BMU fit + grid-mean delta below"],
        ["PIMC polynomial extrapolated to 253.15 K", f"{low['pimc13']:.3f}", f"{low['pimc12']:.3f}",
         "outside the published 270-800 K range"],
        ["theory envelope min", f"{low['env13_lo']:.3f}", f"{low['env12_lo']:.3f}",
         "grid-mean bridges, match at 270 K, and extrapolated PIMC"],
        ["theory envelope max", f"{low['env13_hi']:.3f}", f"{low['env12_hi']:.3f}",
         "systematic; not in the fit error bars"],
        ["local slope (permil per K)", f"{low['slope13']:.4f}", f"{low['slope12']:.4f}",
         "central difference, 0.25 K"],
        ["change for +2 C / -2 C",
         f"{2*low['slope13']:+.3f} / {-2*low['slope13']:+.3f}",
         f"{2*low['slope12']:+.3f} / {-2*low['slope12']:+.3f}",
         "linear local approximation at -20 C"],
        ["AP613 Good -20 C durations", ", ".join(low["durations"]), "",
         "same string means the CSV does not show two tubes"],
        ["AP613 IVW mean of the Good absolute values", f"{low['mean']:.3f}", f"{low['se_mean']:.3f}",
         "measurement SE only"],
        ["higher Good aliquot", f"{low['higher']:.3f}", f"{low['higher_se']:.3f}",
         "selecting the higher of two rows makes this an optimistic point"],
        ["2021 absolute-scale SE used in the bound", f"{low['anchor_se']:.3f}", "",
         "conversion / anchor uncertainty, Delta13"],
        ["bound from the higher aliquot minus 2 combined SE", f"{low['bound_hi']:.3f}", "",
         "the previously quoted style of bound"],
        ["bound from the IVW mean minus a one-sided 95% SE", f"{low['bound_mean']:.3f}", "",
         "use this if the two rows are replicates"],
        ["recommended lower bound", f"{low['recommended_bound']:.3f}", "",
         "the more conservative of the two bounds; a bound, not a measurement of equilibrium"],
        ["apparent T of the recommended bound on the reference curve",
         "" if low["T_bound"] is None else f"{low['T_bound']:.1f} C", "",
         "nonlinear inversion of the reference curve"],
        ["apparent T where the reference curve equals extrapolated PIMC(253 K)",
         "" if low["T_pimc_on_ref"] is None else f"{low['T_pimc_on_ref']:.1f} C", "",
         "cold end of the theory bracket"],
        ["apparent T of the theory-envelope ends",
         "" if low["T_env_lo"] is None else f"{low['T_env_lo']:.1f} C",
         "" if low["T_env_hi"] is None else f"{low['T_env_hi']:.1f} C",
         "Delta13 column is the cold end (high anomaly); Delta12 column is the warm end"],
        ["AP613 mean Delta13 residual, 2025 offset", f"{d13[0]:+.3f}", f"{d13[1]:.3f}",
         f"IVW, offset SE included, {d13[2]:+.1f} sigma; this is the primary deficit"],
        ["AP613 mean Delta13 residual, extended offset", f"{d13x[0]:+.3f}", f"{d13x[1]:.3f}",
         f"secondary, {d13x[2]:+.1f} sigma; do not quote this in place of the campaign offset"],
        ["AP613 mean Delta12 residual, 2025 offset", f"{d12[0]:+.3f}", f"{d12[1]:.3f}",
         f"offset SE included, {d12[2]:+.1f} sigma; measurement-only z in Table 2 is larger"],
    ]
    for r in low["ap_rows"]:
        rows.append([
            f"AP613 {r['duration']} fraction of the cylinder-to-reference distance",
            "" if frac_covered(r, "13") is None else f"{frac_covered(r, '13'):.3f} ± {frac_se(r, '13'):.3f}",
            "" if frac_covered(r, "12") is None else f"{frac_covered(r, '12'):.3f} ± {frac_se(r, '12'):.3f}",
            "Eq. (5); SE includes measurement and offset, not theory",
        ])
    rows.append([
        "CAL1549 -20 C, withheld",
        f"fit-frame absolute {low['cal_fit_absolute']:.2f}; raw column {low['cal_raw']:.2f}; "
        f"fraction {low['cal_fraction']:.2f} ± {low['cal_fraction_se']:.2f}",
        f"fit-frame absolute {low['cal_fit_absolute12']:.2f}; raw column {low['cal_raw12']:.2f}; "
        f"fraction {low['cal_fraction12']:.2f} ± {low['cal_fraction12_se']:.2f}",
        "not a calibration point and not evidence that Delta13 agrees with theory",
    ])
    write_table("TableS6",
                ["quantity", "Δ13CH3D", "Δ12CH2D2", "note"], rows,
                footnotes=[
                    "Equilibrium at -20 C is not demonstrated.  Under a monotonic approach from below, an incomplete experiment is only a lower bound on the equilibrium anomaly.",
                    "The apparent temperatures invert the reference curve.  They are not -20 plus (7.77-6.1)/0.046.",
                ])

def tableS11():
    rows = []
    for (k, iso), fit in sorted(NO550.items()):
        rows.append([NAME[k], ISO_LAB[iso], "T>=270 K, drop 550 C",
                     f3(fit["m"]), f3(fit["se_m"]), fit["n"], f4(fit["p"]),
                     f2(fit["off"]), f2(fit["se_off"]),
                     f"{fit['delta_m']:+.3f}", f"{fit['delta_off']:+.3f}"])
    for iso, items in LOO.items():
        for item in items:
            rows.append(["AP613 extended", ISO_LAB[iso], item["label"] + " (all Good)",
                         f3(item["m"]), f3(item["se"]), item["n"], f4(item["p"]),
                         "", "", f"{item['delta_m']:+.3f}", ""])
    write_table("TableS11",
                ["data set", "isotopologue", "sensitivity", "slope m", "SE(m)", "n", "p",
                 "offset", "SE offset", "change in slope", "change in offset"], rows,
                footnotes=[
                    "Drop-550 refits the T>=270 K line without the AP613 550 C tubes, whose bulk delta-D differs by about 35 permil from the other AP613 tubes.  Offsets for AP613 2025 are in the absolute frame; the working-gas offset moves by the same amount.",
                    "The -20 C leave-one-out is on the all-Good AP613-extended fit, which is a sensitivity.  Dropping either point alone leaves the other at the same temperature, so the slope hardly moves; dropping both releases it.  The two single-point changes are not supposed to sum to the joint change.",
                ])

def tableS12():
    rows = []
    for T in (-20.0, 1.2, 25.0, 50.0, 100.0, 300.0):
        TK = T + 273.15
        out = [f"{T:g}", f"{TK:.2f}"]
        for iso in ("13", "12"):
            slope = dDdT(iso, TK)
            out += [f"{float(REF[iso](TK)):.3f}", f"{slope:.4f}",
                    f"{float(REF[iso](TK + 2) - REF[iso](TK)):.3f}",
                    f"{float(REF[iso](TK - 2) - REF[iso](TK)):.3f}"]
        rows.append(out)
    write_table("TableS12",
                ["T (C)", "T (K)",
                 "ref Δ13", "dΔ13/dT", "Δ13 at +2 C minus ref", "Δ13 at -2 C minus ref",
                 "ref Δ12", "dΔ12/dT", "Δ12 at +2 C minus ref", "Δ12 at -2 C minus ref"],
                rows,
                footnotes=[
                    "Above 270 K the reference is the PIMC polynomial.  At -20 C it is the bridged BMU curve.  A statement that ±2 C moves Delta12 by at most 0.23 permil is not true at -20 C; use this table.",
                    "Temperature error is not included in the fit weights.",
                ])

def tableS5():
    rows = []
    for iso in ("13", "12"):
        for v, lab in (("all", "all"), ("valid", "T>=270K"), ("alt", "alt low-T curve")):
            c = CAMP[(iso, v)]; an = c["ancova"]
            rows.append([ISO_LAB[iso], lab,
                         f2(c["das"]["off"]), f2(c["das"]["se_infl"]), c["das"]["n"],
                         f2(c["siv"]["off"]), f2(c["siv"]["se_infl"]), c["siv"]["n"],
                         f2(c["diff"]), f2(c["se_diff"]), f"{c['z']:+.2f}",
                         f3(an["m"]), f3(an["se_m"]),
                         f2(an["dc"]), f2(an["se_dc"]), f4(an["p_dc"])])
    write_table("TableS5",
                ["isotopologue", "case", "2025 offset (‰)", "SE", "n",
                 "2021 offset (‰)", "SE", "n", "difference (‰)", "SE diff", "z",
                 "ANCOVA m", "SE(m)", "ANCOVA intercept diff (‰)", "SE", "p"], rows)

def tableS7():
    # Bad rows were removed at load; only rows withheld from fits remain here.
    header = ["data set", "T (°C)", "duration", "status", "Δ13 (‰)", "resid13 (‰)", "z13",
              "Δ12 (‰)", "resid12 (‰)", "z12"]
    rows = []
    for r in ROWS:
        if not is_withheld(r):
            continue
        out = [ds_label(r), r["T_C"], r["duration"], "withheld from fits"]
        for iso in ("13", "12"):
            v, e = val(r, iso), sev(r, iso)
            rr = resid(r, iso)
            out += [f2(v), f2(rr), ("" if rr is None or not e else f"{rr / e:+.1f}")]
        rows.append(out)
    write_table("TableS7", header, rows, footnotes=[
        "Rows withheld from every fit.  Bad (QC-failed) rows are listed nowhere: they are dropped at load, before any fit, figure or table.",
    ])

def tableS8():
    rows = []
    for lab in ("apex", "das", "cal"):
        for iso in ("13", "12"):
            for variant in ("repro", "internal"):
                for T, c in sorted(T_CONF[(lab, iso, variant)].items()):
                    rows.append([NAME[lab], ISO_LAB[iso], variant, f"{T:g}",
                                 f"{c['minus']:.1f}", f"{c['plus']:.1f}"])
    write_table("TableS8",
                ["offset", "isotopologue", "precision variant", "true T (°C)",
                 "minus (°C, 68%)", "plus (°C, 68%)"], rows)

def tableS9():
    # if the JSON is available, reproduce its benchmark verbatim (published values)
    if R and "eld_benchmark" in R:
        rows = []
        for iso in ("13", "12"):
            for mode, lab in (("se", "s.e. weights"), ("sigma", "1σ weights")):
                b = R["eld_benchmark"][f"{iso}_{mode}"]
                rows.append([ISO_LAB[iso], lab, f2(b["wg_offset"]), f2(b["wg_se"]),
                             f2(b["pub_wg"]), f2(b["pub_wg_se"]),
                             f3(b["m_lin"]), f3(b["m_lin_se"]), f3(b["m_ln"]), f3(b["m_ln_se"]),
                             f2(b["pub_m"]), f2(b["pub_m_se"])])
        write_table("TableS9",
                    ["isotopologue", "weights", "working-gas offset, this code (‰)", "SE",
                     "published working-gas Δ (‰)", "SE", "slope linear Δ", "SE",
                     "slope 1000 ln", "SE", "published slope", "SE"], rows,
                    footnotes=[
                        "Consistency check of the arithmetic.  Eldridge's house-gas anomalies were assigned from this PIMC curve, so an offset near zero does not test the curve.",
                        "The Delta12 slope from these Table-2 means is not the published slope (0.98 ± 0.05), which was fit to individual analyses.  The p value on the mean-fit is not evidence against PIMC and the mean-fit is not in the primary slope family.",
                    ])
        return
    # otherwise compute from the CSV with s.e. weights; published working-gas
    # composition of Eldridge et al. (2019) house methane: 2.59 / 5.86 permil
    PUB_WG = {"13": 2.59, "12": 5.86}
    rows = []
    for iso in ("13", "12"):
        pts = []
        for q in PTS[("eld", iso)]:
            pts.append({**q, "y": q["y"] - PUB_WG[iso]})   # working-gas frame
        m1 = fit_m1(subset(pts, "all"), REF[iso])
        m_lin = fit_m2(subset(pts, "all"), REF[iso], "lin")
        m_ln = fit_m2(subset(pts, "all"), REF[iso], "ln")
        rows.append([ISO_LAB[iso], "s.e. weights", f2(-m1["off"]), f2(m1["se_infl"]),
                     f2(PUB_WG[iso]), "", f3(m_lin["m"]), f3(m_lin["se_m"]),
                     f3(m_ln["m"]), f3(m_ln["se_m"]), "", ""])
    write_table("TableS9",
                ["isotopologue", "weights", "working-gas Δ, this code (‰)", "SE",
                 "published working-gas Δ (‰)", "SE", "slope linear Δ", "SE",
                 "slope 1000 ln", "SE", "published slope", "SE"], rows)
    print("[note] Table S9: no results_summary.json - published slopes left blank; "
          "1σ-weight row needs the replicate scatter (not in the CSV)")

# ─────────────────────────────────────────────────────────────────────────────
# 9. RUN
# ─────────────────────────────────────────────────────────────────────────────
def method_identity_error():
    """c = o - (m-1)*x_w must hold for the point estimates.  SEs differ."""
    worst = 0.0
    for k in ALLKEYS:
        for iso in ("13", "12"):
            pts = subset(PTS[(k, iso)], "valid")
            if len(pts) < 3:
                continue
            sw = sx = 0.0
            for q in pts:
                x = float(REF[iso](q["T"] + 273.15))
                w = 1.0 / q["se"] ** 2
                sw += w
                sx += w * x
            pred = M1[(k, iso, "valid")]["off"] - (M2[(k, iso, "valid")]["m"] - 1.0) * (sx / sw)
            worst = max(worst, abs(pred - M2[(k, iso, "valid")]["c"]))
    return worst

_ident = method_identity_error()
print(f"[check] Method-1 / Method-2 identity, max |c - prediction| = {_ident:.2e}")
if _ident > 1e-8:
    sys.exit("Method-1 offset and Method-2 intercept are inconsistent")

if R and "theory" in R:
    for key, got in (("ref253_13", float(REF["13"](253.15))),
                     ("ref253_12", float(REF["12"](253.15)))):
        old = R["theory"].get(key)
        if old is not None and abs(float(old) - got) > 0.02:
            print(f"[json] ignoring {key}={old}; this run uses {got:.4f}")

fig1(); fig2(); figS1(); figS2(); figS3()
LOW_T = low_temperature_summary()
table1(); table2_and_S6(); tableS1(); tableS2(); tableS3(); tableS4(); tableS5()
tableS6(LOW_T); tableS7(); tableS8(); tableS9(); tableS10(); tableS11(); tableS12()

def pool_working_gas(iso):
    return ivw((OFF["siv"][iso][0], OFF_WG[iso]), (OFF["siv"][iso][1], SE_WG[iso]))

primary_ps = [M2[(k, iso, "valid")]["p"] for k in ("siv", "das", "cal") for iso in ("13", "12")]
n_empty_good = sum(1 for r in ROWS if is_good(r) and val(r, "13") is None and val(r, "12") is None)
camp = CAMP[("13", "valid")]
camp_lo = camp["diff"] - inv_norm(0.975) * camp["se_diff"]
camp_hi = camp["diff"] + inv_norm(0.975) * camp["se_diff"]
m12, se12, z12, lo12, hi12, _ = _pool_slopes(("siv", "das", "cal"), "12")
anc = camp["ancova"]
anc_lo = anc["dc"] - inv_norm(0.975) * anc["se_dc"]
anc_hi = anc["dc"] + inv_norm(0.975) * anc["se_dc"]

manifest_rows = [
    ["Generator", "AMT outputs from make_AMT_outputs.py (Good-only revision).  Setpoint of the cold experiments: -20 C (253.15 K)."],
    ["QC filter (Good only)",
     f"{len(ROWS)} Good rows are used in every fit, figure and table.  {len(BAD_ROWS)} Bad rows dropped at load: "
     + ("; ".join(f"{ds_label(r)} T={r['T_C']} C" for r in BAD_ROWS) if BAD_ROWS else "none") + "."],
    ["Input CSV", args.csv],
    ["Optional JSON", (json_path if R else "none") + " (Eldridge published benchmark only; curve recomputed)"],
    ["Primary slope family",
     "6 tests: Method 2 at T>=270 K for AP613 2021, AP613 2025 and CAL1549, both isotopologues.  "
     f"Adjusted p = min(1, 6p).  Smallest unadjusted p = {min(primary_ps):.4f}; "
     f"smallest adjusted p = {min(1.0, N_PRIMARY * min(primary_ps)):.4f}.  "
     "AP613 extended and Eldridge are not in the family."],
    ["Not a calibration",
     "Good points below 270 K are excluded from offsets and primary slopes (kinetic constraint, not equilibrium).  "
     "CAL1549 at -20 C is withheld from every fit but kept in the tables and figures.  "
     "Bad rows are excluded at load.  "
     f"Empty Good rows with no anomaly: {n_empty_good}."],
    ["BMU fit",
     "No constant term.  "
     + ", ".join(f"c{k+1}(Δ13)={C_BMU13[k]:.6e}" for k in range(5))
     + ".  "
     + ", ".join(f"c{k+1}(Δ12)={C_BMU12[k]:.6e}" for k in range(5))
     + f".  rms Δ13 {BMU_RMS['13']:.4f}‰, Δ12 {BMU_RMS['12']:.4f}‰."],
    ["Bridge delta, integer K 283 through 423", f"Δ13 {BD13:+.4f}‰, Δ12 {BD12:+.4f}‰"],
    ["Reference-curve step at the 270 K join",
     f"Δ13 {JOIN_STEP['13']:+.4f}‰, Δ12 {JOIN_STEP['12']:+.4f}‰ (above minus below).  Systematic: the constant-delta bridge is not continuous at 270 K; the theory envelope in Table S6 brackets it."],
    ["Reference at -20 C", f"Δ13 {LOW_T['ref13']:.3f}‰, Δ12 {LOW_T['ref12']:.3f}‰"],
    ["Theory envelope at -20 C",
     f"Δ13 {LOW_T['env13_lo']:.3f} to {LOW_T['env13_hi']:.3f}‰; "
     f"Δ12 {LOW_T['env12_lo']:.3f} to {LOW_T['env12_hi']:.3f}‰.  "
     "Windows and the extrapolated polynomial are listed in Table S6."],
    ["Raw-to-absolute 2021 shift",
     f"Δ13 {SHIFT['13']:+.3f}‰ (SD {SHIFT_SD['13']:.3f}, n={SHIFT_N['13']}); "
     f"Δ12 {SHIFT['12']:+.3f}‰ (SD {SHIFT_SD['12']:.3f}, n={SHIFT_N['12']}).  "
     "Treated as a fixed frame definition.  Not the same object as minus the Method-1 offset."],
    ["2021 Method-1 working-gas offset (the other possible shift)",
     f"Δ13 {OFF['siv']['13'][0]:+.3f} ± {OFF['siv']['13'][1]:.3f}‰; "
     f"Δ12 {OFF['siv']['12'][0]:+.3f} ± {OFF['siv']['12'][1]:.3f}‰.  "
     f"Δ12 shift minus (-offset) = {SHIFT['12'] + OFF['siv']['12'][0]:+.3f}‰."],
    ["AP613 2025 working-gas offset used in Table 1",
     f"Δ13 {OFF_WG['13']:+.3f} ± {SE_WG['13']:.3f}‰; Δ12 {OFF_WG['12']:+.3f} ± {SE_WG['12']:.3f}‰ "
     "(absolute offset minus the CSV shift)."],
    ["AP613 2025 working-gas offset if the shift is minus the 2021 offset",
     f"Δ13 {OFF_WG_ALT['13']:+.3f} ± {SE_WG_ALT['13']:.3f}‰; "
     f"Δ12 {OFF_WG_ALT['12']:+.3f} ± {SE_WG_ALT['12']:.3f}‰.  "
     "Reported, not used.  Delta13 barely moves; Delta12 moves by the shift discrepancy above."],
    ["Pooled AP613 Δ12 working-gas offset (recommended pool)",
     (lambda m, s: f"{m:+.2f} ± {s:.2f}‰.  Delta13 is not pooled: the campaigns differ, and averaging them is not the calibration.")(*pool_working_gas("12"))],
    ["AP613 Δ13 campaign difference at T>=270 K, 2025 minus 2021, absolute",
     f"{camp['diff']:+.2f} ± {camp['se_diff']:.2f}‰, z={camp['z']:+.2f}, "
     f"normal 95% interval {camp_lo:+.2f} to {camp_hi:+.2f}‰.  "
     f"ANCOVA intercept difference {anc['dc']:+.2f} ± {anc['se_dc']:.2f}‰ "
     f"({anc_lo:+.2f} to {anc_hi:+.2f}).  A drift of a few tenths is the central value and is unresolved."],
    ["Primary Δ12 slope pool (no Eldridge)", f"{m12:.3f} ± {se12:.3f}, z vs 1 = {z12:+.2f}, 95% {lo12:.3f} to {hi12:.3f}"],
    ["AP613 -20 C Δ13",
     f"recommended lower bound {LOW_T['recommended_bound']:.2f}‰ "
     f"(higher-aliquot bound {LOW_T['bound_hi']:.2f}, replicate-mean bound {LOW_T['bound_mean']:.2f}).  "
     f"Primary residual vs the 2025 offset: "
     f"{LOW_T['deficits'][('AP613 2025 offset', '13')][0]:+.2f} ± "
     f"{LOW_T['deficits'][('AP613 2025 offset', '13')][1]:.2f}‰.  "
     f"Apparent T of the bound: {LOW_T['T_bound']:.1f} C.  "
     f"Apparent T of extrapolated PIMC on the reference curve: {LOW_T['T_pimc_on_ref']:.1f} C.  "
     "See Table S6.  Not a calibration anchor."],
    ["550 C sensitivity",
     "Table S11.  Bulk delta-D at 550 C is not explained by this CSV; the refit is the quantitative check."],
    ["Temperature error",
     "±2 C is not in the fit.  At -20 C, +2 C moves the reference by "
     f"{LOW_T['slope13']*2:+.3f}‰ (Δ13) and {LOW_T['slope12']*2:+.3f}‰ (Δ12).  Table S12."],
    ["Table S8",
     "68% half-widths, not one interval 'from 7 to 13 C'.  "
     f"Reproducibility {args.repro13}/{args.repro12}‰ can be tighter than the internal SE "
     f"{args.internal13}/{args.internal12}‰ because some chi2r values are below 1."],
    ["Frames",
     "Figure 1 is entirely in the working-gas frame.  AP613 2025 and Eldridge fits are absolute; "
     "Figure S1 keeps each series in its own frame and draws the T>=270 K line as the primary fit."],
    ["Cylinder start vs reference (absolute frame)",
     f"AP613 start = +{SHIFT['13']:.3f}/+{SHIFT['12']:.3f} permil; "
     f"CAL1549 start = +{-OFF['cal']['13'][0]:.3f}/+{-OFF['cal']['12'][0]:.3f} permil.  "
     "Both cylinders lie below the reference curve over the whole experimental range, so the approach "
     "to equilibrium is from below; the -20 C lower-bound argument assumes this monotonic approach from below."],
    ["Figures",
     "Fig1_Delta_vs_T, Fig2_residuals_vs_T, FigS1_method2_regression, "
     "FigS2_shifted_by_offset, FigS3_low_temperature (.pdf + .png)"],
]

workbook_path = write_workbook(manifest_rows)
legacy_text_outputs = [
    "Table1_calibration_Tge270K.txt", "Table2_experiments_le25C.txt",
    "TableS1_protocol.txt", "TableS2_method1_offsets.txt", "TableS3_method2_allGood.txt",
    "TableS4_slopes_five_choices.txt", "TableS5_campaign_comparison.txt",
    "TableS6_experiments_le25C.txt", "TableS7_excluded_points.txt",
    "TableS8_temperature_uncertainty.txt", "TableS9_eldridge_benchmark.txt",
    "TableS10_slope_summary.txt", "manifest.txt",
]
for filename in legacy_text_outputs:
    old_path = os.path.join(OUT, filename)
    if os.path.isfile(old_path):
        os.remove(old_path)
print(f"[done] workbook -> {workbook_path} ({len(TABLES)} tables + Manifest sheet)")
