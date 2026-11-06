#!/usr/bin/env python3
"""
make_AMT_outputs.py - static AMT figures and tables from calibration_data_table.csv.

Self-contained: the ONLY required input is the CSV that sits next to this
script (override with --csv).  No HTML and no Plotly figures are produced,
and HeatExperiments.py is not needed.  If a results_summary.json is found
(--json, or next to the script) it is used ONLY for run config (z_Bonferroni,
reproducibility) and the published Eldridge et al. benchmark values;
everything else is recomputed here from the CSV.

All fits are recomputed from the CSV:
  Method 1  weighted offset, SE inflated by sqrt(chi2/(n-1)) when > 1
  Method 2  weighted line  y = m*x + c  (x = reference-curve anomaly),
            SE inflated by sqrt(chi2r), z = (m-1)/SE(m), p from Student t
  variants  all | valid (T >= 270 K) | alt (extrapolated PIMC low-T curve)
            | ln (1000 ln(1+Delta/1000) space) | incl_bad (Good + Bad)
  campaigns AP613 2025 vs 2021 on the absolute scale: separate Method 1
            offsets + ANCOVA (common slope, one intercept per campaign)
  Table 2/S6/S7 residuals use the offset fitted at T >= 270 K; the Eq. (5)
            fraction is relative to the cylinder-gas starting composition
  Table S8  Monte Carlo inversion of the reference curve (68 % interval),
            with long-term reproducibility and with internal precision

Outputs (written next to the script by default, --out to change):
  Fig1_Delta_vs_T, Fig2_residuals_vs_T, FigS1_method2_regression,
  FigS2_shifted_by_offset            (.pdf vector + .png @ --dpi)
  Table1_calibration_Tge270K.txt, Table2_experiments_le25C.txt,
  TableS1_protocol.txt ... TableS9_eldridge_benchmark.txt, manifest.txt

Usage
  python3 make_AMT_outputs.py                 # CSV and outputs in the script's folder
  python3 make_AMT_outputs.py --csv /path/calibration_data_table.csv --out /path/AMT

Requires: numpy, matplotlib.  (pip install numpy matplotlib)

Author: Bibhasvata Dasgupta | ICIW 2026
"""

import os, sys, json, csv, math, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")                    # headless: no HTML / interactive output
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# 0. ARGUMENTS
# ─────────────────────────────────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser(description="AMT figures+tables from the calibration CSV (no HTML)")
parser.add_argument("--csv", default=os.path.join(HERE, "calibration_data_table.csv"))
parser.add_argument("--out", default=HERE, help="output folder (default: script's folder)")
parser.add_argument("--json", default=None,
                    help="optional results_summary.json (config + published Eldridge values)")
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

# optional JSON: config + published benchmark only
R = None
json_path = args.json or os.path.join(HERE, "results_summary.json")
if os.path.isfile(json_path):
    R = json.load(open(json_path, encoding="utf-8"))
    print(f"[json] optional results_summary.json found: {json_path} (config/published values only)")
else:
    print("[json] no results_summary.json - computing everything from the CSV")

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

# Bonferroni threshold for 8 slope tests at 5 % (two-sided)
Z_BONF = (R or {}).get("config", {}).get("z_bonferroni", inv_norm(1.0 - 0.05 / 16.0))

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

KEYS = ["siv", "das", "cal", "eld"]                 # plotted data sets
ALLKEYS = KEYS + ["apex"]                           # apex = AP613 extended (abs frame)
NAME = {"siv": "AP613 2021 (Sivan et al. 2024)", "das": "AP613 2025 (this study)",
        "cal": "CAL1549 (this study)", "eld": "Eldridge et al. 2019",
        "apex": "AP613 extended (both campaigns, absolute)"}
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
                    "src": key_of(r["dataset"]), "row": r})
    return pts

# raw -> absolute shift of the Sivan data, from the CSV itself (constant to ~0.03 permil)
SHIFT = {}
for iso, rawc, corc in (("13", "D13_raw", "D13"), ("12", "D12_raw", "D12")):
    d = [num(r[corc]) - num(r[rawc]) for r in rows_of("siv")
         if num(r[rawc]) is not None and num(r[corc]) is not None]
    SHIFT[iso] = float(np.mean(d)) if d else 0.0

# ─────────────────────────────────────────────────────────────────────────────
# 3. THEORY CURVES
#    PIMC polynomials of Eldridge et al. 2019 (valid 270-800 K, T in K);
#    BMU = 5th-order polynomial in 1/T fitted to their Table 4;
#    reference = PIMC for T >= 270 K, else BMU + bridge delta (mean PIMC-BMU
#    over 283-423 K; overridden by results_summary.json when available).
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

_bmu_TK = BMU_T_C + 273.15
_x = 1.0 / _bmu_TK
P_BMU13 = np.polyfit(_x, BMU_13, 5)
P_BMU12 = np.polyfit(_x, BMU_12, 5)
_m = (_bmu_TK >= 283) & (_bmu_TK <= 423)
_th = (R or {}).get("theory", {})
BD13 = _th.get("bridge_delta13", float(np.mean(pimc13(_bmu_TK[_m]) - BMU_13[_m])))
BD12 = _th.get("bridge_delta12", float(np.mean(pimc12(_bmu_TK[_m]) - BMU_12[_m])))

def bmu13(TK): return np.polyval(P_BMU13, 1.0 / np.asarray(TK, dtype=float))
def bmu12(TK): return np.polyval(P_BMU12, 1.0 / np.asarray(TK, dtype=float))

def ref13(TK):
    TK = np.asarray(TK, dtype=float)
    return np.where(TK >= 270.0, pimc13(TK), bmu13(TK) + BD13)

def ref12(TK):
    TK = np.asarray(TK, dtype=float)
    return np.where(TK >= 270.0, pimc12(TK), bmu12(TK) + BD12)

REF = {"13": ref13, "12": ref12}     # reference curve (bridged BMU below 270 K)
ALT = {"13": pimc13, "12": pimc12}   # the other low-T curve: PIMC extrapolated
BMU = {"13": bmu13, "12": bmu12}
PIMC = ALT

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
    if variant == "incl_bad":
        return list(pts)
    p = [q for q in pts if q["good"]]
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
    return {"m": m, "se_m": se_m, "c": c, "se_c": se_c, "z": z,
            "p": t_two_sided_p(abs(z), n - 2), "chi2r": chi2r, "n": n}

VARIANTS = ["all", "valid", "alt", "ln", "incl_bad"]
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
# AP613 2025 working-gas offset via the Sivan raw->absolute shift (Sect. S4)
OFF_WG = {iso: OFF["das"][iso][0] - SHIFT[iso] for iso in ("13", "12")}
SE_WG = {iso: OFF["das"][iso][1] for iso in ("13", "12")}

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
    to the reference value covered.  Working-gas frame: start = 0; absolute
    frame: start = -offset.  Undefined where the frame has no dynamic range."""
    v = val(row, iso)
    if v is None:
        return None
    k = key_of(row["dataset"])
    off = OFF[k][iso][0]
    ref = float(REF[iso](float(row["T_C"]) + 273.15))
    den = ref + off
    if abs(den) < 1e-9:
        return None
    start = 0.0 if FRAME[k].startswith("working") else -off
    return (v - start) / den

# ─────────────────────────────────────────────────────────────────────────────
# 6. TABLE S8: MONTE CARLO INVERSION OF THE REFERENCE CURVE (68 % interval)
# ─────────────────────────────────────────────────────────────────────────────
T_GRID = [float(t) for t in args.t_grid.split(",") if t.strip()]
_TK_GRID = np.linspace(250.0, 830.0, 20001)

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

def draw_points(ax, iso):
    for k in KEYS:
        for r in rows_of(k):
            v = val(r, iso)
            if v is None:
                continue
            good = is_good(r)
            ax.scatter([float(r["T_C"])], [v], marker=MARK[k], s=42,
                       facecolors=COLOR[k] if good else "none",
                       edgecolors=COLOR[k], linewidths=1.4, zorder=3)

def fig1():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharex=True)
    for ax, iso in zip(axes, ("13", "12")):
        draw_theory(ax, iso)
        for k in KEYS:
            ts = [float(r["T_C"]) for r in rows_of(k)]
            tg = np.linspace(min(ts) + 273.15, max(ts) + 273.15, 300)
            off, se = OFF[k][iso]
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
    fig.suptitle("Figure 1.  Δ against equilibration temperature; open symbols = excluded points; "
                 "grey region: T < 270 K", fontsize=10)
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
                good = is_good(r)
                ax.errorbar([float(r["T_C"])], [rr], yerr=[e], fmt="none",
                            ecolor=COLOR[k], elinewidth=1.0, alpha=0.8, zorder=2)
                ax.scatter([float(r["T_C"])], [rr], marker=MARK[k], s=42,
                           facecolors=COLOR[k] if good else "none",
                           edgecolors=COLOR[k], linewidths=1.4, zorder=3)
        ax.set_xlabel("equilibration temperature (°C)")
        ax.set_ylabel(f"residual {ISO_LAB[iso]} (‰)")
        ax.set_xlim(XMIN, XMAX)
    fig.suptitle("Figure 2.  Residuals from the reference curve plus the offset fitted at T ≥ 270 K; "
                 "error bars 1 s.e.; open symbols = excluded; grey region: T < 270 K", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "Fig2_residuals_vs_T")

def figS1():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2))
    for ax, iso in zip(axes, ("13", "12")):
        for k in KEYS:
            y = M2[(k, iso, "all")]
            xs, ys, low = [], [], []
            for r in rows_of(k):
                v = val(r, iso)
                if v is None:
                    continue
                x = float(REF[iso](float(r["T_C"]) + 273.15))
                xs.append(x); ys.append(v)
                if float(r["T_C"]) < LOWT_C:
                    low.append((x, v))
            xs = np.array(xs); ys = np.array(ys)
            ax.scatter(xs, ys, marker=MARK[k], s=42, facecolors=COLOR[k],
                       edgecolors=COLOR[k], linewidths=1.2, zorder=3)
            if low:
                ax.scatter([p[0] for p in low], [p[1] for p in low], marker=MARK[k], s=90,
                           facecolors="none", edgecolors="k", linewidths=1.6, zorder=4)
            xg = np.linspace(xs.min() - 0.2, xs.max() + 0.2, 100)
            ax.plot(xg, y["m"] * xg + y["c"], color=COLOR[k], lw=1.6,
                    label=(f"{NAME[k]}: m={y['m']:.3f}±{y['se_m']:.3f}, c={y['c']:+.2f}, "
                           f"χ²r={y['chi2r']:.1f}, z={y['z']:+.2f}, p={y['p']:.3f}"))
            band = np.sqrt((y["se_m"] * xg) ** 2 + y["se_c"] ** 2)
            ax.fill_between(xg, y["m"] * xg + y["c"] - band, y["m"] * xg + y["c"] + band,
                            color=COLOR[k], alpha=0.15, lw=0)
        ax.set_xlabel(f"reference-curve {ISO_LAB[iso]} (‰)")
        ax.set_ylabel(f"measured {ISO_LAB[iso]} (‰)")
        ax.legend(fontsize=7.2, loc="best", framealpha=0.9)
    fig.suptitle("Figure S1.  Method 2: measured against reference-curve anomaly at the same T; "
                 "bands ±1σ of the fitted line; black-ringed points below 270 K", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
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
                good = is_good(r)
                ax.scatter([float(r["T_C"])], [v - off], marker=MARK[k], s=42,
                           facecolors=COLOR[k] if good else "none",
                           edgecolors=COLOR[k], linewidths=1.4,
                           alpha=0.45 if not good else 1.0, zorder=3)
        ax.set_xlabel("equilibration temperature (°C)")
        ax.set_ylabel(f"{ISO_LAB[iso]} shifted by own offset (‰)")
        ax.set_xlim(XMIN, XMAX)
        ax.axvline(LOWT_C, color="grey", lw=0.8, ls="--")
    fig.suptitle("Figure S2.  Each data set shifted by its own Method 1 offset (T ≥ 270 K); "
                 "faded open symbols = excluded; bands ±1σ of the offset; grey region: T < 270 K",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "FigS2_shifted_by_offset")

# ─────────────────────────────────────────────────────────────────────────────
# 8. TABLES (CSV)
# ─────────────────────────────────────────────────────────────────────────────
def write_csv(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=",", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    print(f"[tbl ] {name}")

def f2(x): return "" if x is None else f"{x:.2f}"
def f3(x): return "" if x is None else f"{x:.3f}"
def f4(x): return "" if x is None else f"{x:.4f}"

def table1():
    rows = []
    for k in ("siv", "das", "cal", "apex", "eld"):
        frame = {"siv": "working gas", "das": "working gas (Sect. S4 conversion)",
                 "cal": "working gas", "apex": "absolute", "eld": "absolute"}[k]
        for iso in ("13", "12"):
            v = M2[(k, iso, "valid")]
            off, se = (OFF_WG[iso], SE_WG[iso]) if k == "das" else OFF[k][iso]
            rows.append([NAME[k], frame, ISO_LAB[iso], v["n"], f2(off), f2(se),
                         f3(v["m"]), f3(v["se_m"]), f"{v['z']:+.2f}", f4(v["p"])])
    write_csv("Table1_calibration_Tge270K.txt",
              ["data set", "offset frame", "isotopologue", "n", "offset (‰)", "SE offset (‰)",
               "slope m", "SE(m)", "z", "p"], rows)

def table2_and_S6():
    header = ["data set", "T (°C)", "catalyst", "duration", "flag",
              "Δ13 (‰)", "SE13", "resid13 (‰)", "z13", "frac13",
              "Δ12 (‰)", "SE12", "resid12 (‰)", "z12", "frac12"]
    rows = []
    for r in sorted(ROWS, key=lambda r: (r["dataset"], float(r["T_C"]))):
        if float(r["T_C"]) > 25.0:
            continue
        out = [r["dataset"], r["T_C"], r["catalyst"], r["duration"], r["flag"]]
        for iso in ("13", "12"):
            v, e = val(r, iso), sev(r, iso)
            rr, fr = resid(r, iso), frac_covered(r, iso)
            out += [f2(v), f2(e), f2(rr),
                    ("" if rr is None or not e else f"{rr / e:+.1f}"),
                    ("" if fr is None else f"{fr:.3f}")]
        rows.append(out)
    write_csv("Table2_experiments_le25C.txt", header, rows)
    write_csv("TableS6_experiments_le25C.txt", header, rows)

def tableS1():
    write_csv("TableS1_protocol.txt",
              ["data set", "T (°C)", "catalyst", "duration", "flag"],
              [[r["dataset"], r["T_C"], r["catalyst"], r["duration"], r["flag"]] for r in ROWS])

def tableS2():
    rows = []
    for k in ("siv", "apex", "cal", "eld"):
        for iso in ("13", "12"):
            a, v = M1[(k, iso, "all")], M1[(k, iso, "valid")]
            rows.append([NAME[k], FRAME[k], ISO_LAB[iso],
                         f2(a["off"]), f2(a["se_infl"]), a["n"], f2(a["chi2r"]),
                         f2(v["off"]), f2(v["se_infl"]), v["n"], f2(v["chi2r"])])
    write_csv("TableS2_method1_offsets.txt",
              ["data set", "frame", "isotopologue",
               "offset all (‰)", "SE infl all", "n all", "chi2r all",
               "offset T>=270K (‰)", "SE infl T>=270K", "n T>=270K", "chi2r T>=270K"], rows)

def tableS3():
    rows = []
    for k in ("siv", "apex", "cal", "eld"):
        for iso in ("13", "12"):
            y = M2[(k, iso, "all")]
            rows.append([NAME[k], ISO_LAB[iso], y["n"], f3(y["m"]), f3(y["se_m"]),
                         f"{y['z']:+.2f}", f4(y["p"]), f4(min(1.0, 8.0 * y["p"])),
                         f2(y["c"]), f2(y["se_c"]), f2(y["chi2r"])])
    write_csv("TableS3_method2_allGood.txt",
              ["data set", "isotopologue", "n", "m", "SE(m)", "z", "p", "Bonferroni 8p",
               "intercept c (‰)", "SE(c)", "chi2r"], rows)

def tableS4():
    lab = {"all": "all", "valid": "T>=270K", "alt": "alt low-T curve",
           "ln": "ln space", "incl_bad": "+Bad"}
    rows = []
    for k in ("siv", "apex", "cal", "eld"):
        for iso in ("13", "12"):
            for v in VARIANTS:
                y = M2[(k, iso, v)]
                rows.append([NAME[k], ISO_LAB[iso], lab[v], f3(y["m"]), f3(y["se_m"]),
                             y["n"], f4(y["p"])])
    write_csv("TableS4_slopes_five_choices.txt",
              ["data set", "isotopologue", "choice", "m", "SE(m)", "n", "p"], rows)

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
    write_csv("TableS5_campaign_comparison.txt",
              ["isotopologue", "case", "2025 offset (‰)", "SE", "n",
               "2021 offset (‰)", "SE", "n", "difference (‰)", "SE diff", "z",
               "ANCOVA m", "SE(m)", "ANCOVA intercept diff (‰)", "SE", "p"], rows)

def tableS7():
    header = ["data set", "T (°C)", "duration", "flag", "Δ13 (‰)", "resid13 (‰)", "z13",
              "Δ12 (‰)", "resid12 (‰)", "z12"]
    rows = []
    for r in ROWS:
        if is_good(r):
            continue
        out = [r["dataset"], r["T_C"], r["duration"], r["flag"]]
        for iso in ("13", "12"):
            v, e = val(r, iso), sev(r, iso)
            rr = resid(r, iso)
            out += [f2(v), f2(rr), ("" if rr is None or not e else f"{rr / e:+.1f}")]
        rows.append(out)
    write_csv("TableS7_excluded_points.txt", header, rows)

def tableS8():
    rows = []
    for lab in ("apex", "das", "cal"):
        for iso in ("13", "12"):
            for variant in ("repro", "internal"):
                for T, c in sorted(T_CONF[(lab, iso, variant)].items()):
                    rows.append([NAME[lab], ISO_LAB[iso], variant, f"{T:g}",
                                 f"{c['minus']:.1f}", f"{c['plus']:.1f}"])
    write_csv("TableS8_temperature_uncertainty.txt",
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
        write_csv("TableS9_eldridge_benchmark.txt",
                  ["isotopologue", "weights", "working-gas offset, this code (‰)", "SE",
                   "published working-gas Δ (‰)", "SE", "slope linear Δ", "SE",
                   "slope 1000 ln", "SE", "published slope", "SE"], rows)
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
    write_csv("TableS9_eldridge_benchmark.txt",
              ["isotopologue", "weights", "working-gas Δ, this code (‰)", "SE",
               "published working-gas Δ (‰)", "SE", "slope linear Δ", "SE",
               "slope 1000 ln", "SE", "published slope", "SE"], rows)
    print("[note] Table S9: no results_summary.json - published slopes left blank; "
          "1σ-weight row needs the replicate scatter (not in the CSV)")

# ─────────────────────────────────────────────────────────────────────────────
# 9. RUN
# ─────────────────────────────────────────────────────────────────────────────
fig1(); fig2(); figS1(); figS2()
table1(); table2_and_S6(); tableS1(); tableS2(); tableS3(); tableS4(); tableS5(); tableS7(); tableS8(); tableS9()

with open(os.path.join(OUT, "manifest.txt"), "w", encoding="utf-8") as f:
    f.write("AMT outputs generated by make_AMT_outputs.py from the CSV (no HTML/Plotly)\n")
    f.write(f"input: {args.csv}\noptional json: {json_path if R else 'none'}\n\n")
    f.write(f"z_Bonferroni (8 tests, 5%) = {Z_BONF:.4f}\n")
    f.write(f"raw->absolute shift of Sivan data: Δ13 {SHIFT['13']:+.3f}‰, Δ12 {SHIFT['12']:+.3f}‰\n")
    f.write(f"bridge delta (PIMC-BMU, mean 283-423 K): Δ13 {BD13:+.4f}‰, Δ12 {BD12:+.4f}‰\n")
    f.write(f"reproducibility used for Table S8: repro {args.repro13}/{args.repro12}‰, "
            f"internal {args.internal13}/{args.internal12}‰ (Δ13/Δ12)\n\n")
    f.write("figures: Fig1_Delta_vs_T, Fig2_residuals_vs_T, FigS1_method2_regression, "
            "FigS2_shifted_by_offset (.pdf + .png)\n")
    f.write("tables:  Table1, Table2, TableS1..TableS9 (.txt, comma-separated)\n")
    f.write("frames:  AP613 2021 and CAL1549 in working-gas frame; AP613 2025, AP613 extended\n"
            "         and Eldridge et al. 2019 in the absolute frame (Sect. S4 conversion)\n")
print(f"[done] manifest -> {os.path.join(OUT, 'manifest.txt')}")