"""Generate synthetic datasets + config grid for the cross-language comparison.

data_main.csv and data_weak.csv are the July 2026 datasets (the same draws; since
2026-10 the special FE labels keep their full length, item 38). data_ext.csv
(2026-10, WP8) holds the scenarios added for the fix phase; it is drawn from its
own generator, so the July datasets do not depend on it.
"""
import numpy as np
import os, json

OUT = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(987654321)

# ---------------- main dataset ----------------
n = 900
z1 = rng.normal(size=n)
z2 = rng.normal(size=n)
x1 = rng.normal(size=n)
x2 = 0.4 * x1 + rng.normal(size=n)
x3 = rng.normal(size=n)
reg_a = rng.normal(size=n)
reg_b = rng.normal(size=n)
reg_c = rng.normal(size=n)

Sig = np.array([[1.0, 0.5, 0.2], [0.5, 1.0, 0.3], [0.2, 0.3, 1.0]])
E = rng.multivariate_normal([0, 0, 0], Sig, size=n)
eT, eM, eY = E[:, 0], E[:, 1], E[:, 2]

tvar = 0.8 * z1 + 0.3 * z2 + 0.4 * x1 - 0.2 * x2 + eT
mvar = 0.6 * tvar + 0.5 * x2 + 0.3 * x3 + eM
y = 0.3 * tvar + 0.5 * mvar + 0.4 * x1 + eY

# categorical dummies
dcat = rng.choice(["alpha", "Beta", "beta", "Gamma"], size=n)     # mixed case -> sort-order probe
dnum = rng.choice([1, 2, 3], size=n)                              # integer categories
dflo = rng.choice([1.0, 2.0, 3.0], size=n)                        # float categories -> name-format probe

clu = rng.integers(0, 30, size=n)                                 # 30 clusters
clu_lab = np.array([f"c{int(c):02d}" for c in clu])

# fixed effects (dtype=object so the special labels below keep their full
# length; a '<U3' array cut them to 'gSN', 'gCA', 'hSN', 'hSO' before 2026-10)
fe1 = np.array([f"g{int(g):02d}" for g in rng.integers(0, 45, size=n)], dtype=object)
fe2 = np.array([f"h{int(h):02d}" for h in rng.integers(0, 30, size=n)], dtype=object)
# direct singletons in each dimension
fe1[10] = "gSNG1"                       # singleton in fe1
fe2[20] = "hSNG2"                       # singleton in fe2
# cascade: rows 30,31 alone in fe1 group gCASC; row 30 is alone in fe2 group hSOLO.
# Dropping fe2 singletons removes row 30 -> gCASC becomes a NEW fe1 singleton (row 31).
fe1[30] = "gCASC"; fe1[31] = "gCASC"
fe2[30] = "hSOLO"
# make sure rows 31.. keep multi-member fe2 groups (leave as drawn)

w = rng.uniform(0.5, 2.5, size=n)
w[[100, 200, 300, 400]] = 0.0                                     # zero-weight probe

lat = rng.uniform(36.0, 48.0, size=n)
lon = rng.uniform(-100.0, -90.0, size=n)

# missing values (disjoint from special rows 10,20,30,31,100..400 used above)
x2_m = x2.copy(); x2_m[[50, 51, 52, 53, 54, 55, 56, 57, 58, 59]] = np.nan
mvar_m = mvar.copy(); mvar_m[[60, 61, 62, 63, 64, 65, 66, 67]] = np.nan
z2_m = z2.copy(); z2_m[[70, 71, 72, 73, 74]] = np.nan

cols = dict(y=y, tvar=tvar, mvar=mvar_m, z1=z1, z2=z2_m,
            x1=x1, x2=x2_m, x3=x3, reg_a=reg_a, reg_b=reg_b, reg_c=reg_c,
            dcat=dcat, dnum=dnum, dflo=dflo, clu=clu_lab, fe1=fe1, fe2=fe2,
            w=w, lat=lat, lon=lon)


def write_csv(path, cols):
    """None -> empty field (missing text id), NaN -> empty field, +-inf -> Inf/-Inf."""
    names = list(cols.keys())
    m = len(next(iter(cols.values())))
    with open(path, "w", newline="") as f:
        f.write(",".join(names) + "\n")
        for i in range(m):
            row = []
            for nm in names:
                v = cols[nm][i]
                if v is None:
                    row.append("")
                elif isinstance(v, (str, np.str_)):
                    row.append(str(v))
                elif isinstance(v, (int, np.integer)):
                    row.append(str(int(v)))
                else:
                    fv = float(v)
                    if np.isnan(fv):
                        row.append("")
                    elif np.isinf(fv):
                        row.append("Inf" if fv > 0 else "-Inf")
                    else:
                        row.append(repr(fv))
            f.write(",".join(row) + "\n")


write_csv(os.path.join(OUT, "data_main.csv"), cols)

# ---------------- weak-IV dataset ----------------
# rho_TM = 0.10 and gamma = 0 -> F_{M,Z|T} ~ 0.9 -> AR theta set is 'rays'
# (unbounded weak-IV signature) under both robust and cluster SEs.
n2 = 400
rng2 = np.random.default_rng(24681357)
SigW = np.array([[1.0, 0.10, 0.2], [0.10, 1.0, 0.3], [0.2, 0.3, 1.0]])
z1w = rng2.normal(size=n2)
E2 = rng2.multivariate_normal([0, 0, 0], SigW, size=n2)
tW = 1.0 * z1w + E2[:, 0]
mW = 0.0 * tW + E2[:, 1]            # weak conditional first stage for theta
yW = 0.10 * tW + 0.5 * mW + E2[:, 2]
x1w = rng2.normal(size=n2)
cluw = np.array([f"k{int(c):02d}" for c in rng2.integers(0, 25, size=n2)])
write_csv(os.path.join(OUT, "data_weak.csv"),
          dict(y=yW, tvar=tW, mvar=mW, z1=z1w, x1=x1w, clu=cluw))

# ---------------- extension dataset (2026-10, WP8) ----------------
# One wide file; each config picks its columns (Y, T, M default to y, tvar, mvar).
n3 = 1200
rng3 = np.random.default_rng(20261018)   # chosen so the K = 2 AR sets are bounded
z1e = rng3.normal(size=n3)
z2e = rng3.normal(size=n3)
Xe = rng3.normal(size=(n3, 8))                  # x1..x8
Xe[:, 1] += 0.5 * Xe[:, 0]
Xe[:, 3] += 0.3 * Xe[:, 2]
Ce = rng3.normal(size=(n3, 20))                 # c01..c20 (many small coefficients)
bC = np.array([0.5, -0.4, 0.3, 0.25, -0.2, 0.15, 0.1, 0.08, 0.06, 0.05] + [0.0] * 10)
het = np.exp(0.3 * z1e)                         # heteroskedasticity tied to z1: K_eff != K for K = 2
# Cov(e_T, e_M) = 0.7 gives a strong conditional first stage; Cov(e_T, e_Y) = 0
SigE = np.array([[1.0, 0.7, 0.0], [0.7, 1.0, 0.3], [0.0, 0.3, 1.0]])
Ee = rng3.multivariate_normal([0, 0, 0], SigE, size=n3) * het[:, None]
# T has no C terms, so omitted controls never confound e_T (the K = 2 model stays
# correctly specified under the heteroskedasticity and its AR set is bounded)
tE = 0.5 * z1e + 0.35 * z2e + 0.3 * Xe[:, 0] - 0.2 * Xe[:, 1] + 0.2 * Xe[:, 2] + Ee[:, 0]
mE = 0.5 * tE + 0.4 * Xe[:, 1] + 0.3 * Xe[:, 3] + 0.4 * (Ce @ bC[::-1]) + Ee[:, 1]
yE = 0.25 * tE + 0.45 * mE + 0.4 * Xe[:, 0] + 0.3 * Xe[:, 4] + 0.3 * (Ce @ bC) + Ee[:, 2]
yneg = 0.25 * tE - 0.6 * mE + 0.4 * Xe[:, 0] + 0.3 * Xe[:, 4] + Ee[:, 2]   # theta < 0

we = rng3.uniform(0.5, 2.5, size=n3)

# one-way / two-way FE with long labels
fa_c = rng3.integers(0, 40, size=n3)
fb_c = rng3.integers(0, 25, size=n3)
cl_c = rng3.integers(0, 24, size=n3)
fan_c = cl_c * 3 + rng3.integers(0, 3, size=n3)          # nested in cl (72 groups)
fa = np.array([f"group_a_{c:02d}" for c in fa_c], dtype=object)
fb = np.array([f"group_b_{c:02d}" for c in fb_c], dtype=object)
cl = np.array([f"cluster_{c:02d}" for c in cl_c], dtype=object)
fan = np.array([f"nested_{c:03d}" for c in fan_c], dtype=object)

# disconnected two-way FE: rows 0..599 and 600..1199 share no fd1 or fd2 level (2 components)
half = n3 // 2
fd1_c = np.concatenate([rng3.integers(0, 15, size=half), rng3.integers(15, 30, size=n3 - half)])
fd2_c = np.concatenate([rng3.integers(0, 10, size=half), rng3.integers(10, 20, size=n3 - half)])
fd1 = np.array([f"disc1_{c:02d}" for c in fd1_c], dtype=object)
fd2 = np.array([f"disc2_{c:02d}" for c in fd2_c], dtype=object)

# zero weights with FE: rows (0,1) and (2,3) form two-row groups whose zero-weight
# member must be dropped BEFORE singleton detection (rows 1 and 2 become singletons)
fz = np.array([f"zgroup_{c:02d}" for c in rng3.integers(0, 50, size=n3)], dtype=object)
fz[[0, 1]] = "zgroup_pair_a"
fz[[2, 3]] = "zgroup_pair_b"
wz = we.copy()
wz[[0, 3, 50, 150, 250, 350, 450, 550, 650, 750, 850, 950]] = 0.0

# missing text ids: "" and whitespace-only are missing in every language (WP5.1)
cl_txt = cl.copy()
fa_txt = fa.copy()
dtxt = np.array([f"cat_{c}" for c in rng3.choice(["a", "b", "c"], size=n3)], dtype=object)
for i in range(100, 110):
    cl_txt[i] = None                     # empty field
for i in range(110, 115):
    cl_txt[i] = " "                      # whitespace only
for i in range(120, 130):
    fa_txt[i] = None
for i in range(130, 135):
    fa_txt[i] = "  "
for i in range(135, 140):
    dtxt[i] = None

# bootstrap with missing values (a missing mediator and a missing control)
mmiss = mE.copy(); mmiss[[160, 161, 162, 163, 164, 165, 166, 167]] = np.nan
xmiss = Xe[:, 2].copy(); xmiss[[170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181]] = np.nan

# ill-conditioned and exactly collinear controls
year = rng3.integers(1950, 2020, size=n3)
year2 = year * year
x12 = Xe[:, 0] + Xe[:, 1]                       # exactly x1 + x2 -> dropped (D9)

# shifted copies (LASSO without FE must not depend on the location of any input)
shifts = 1.0 + np.arange(8)

# degenerate inputs (error_configs)
gz = rng3.normal(size=40)
zg = gz[fa_c]                                   # constant within fa
tconst = np.ones(n3)
z1dbl = 2.0 * z1e                               # collinear with z1
xinf = Xe[:, 5].copy(); xinf[5] = np.inf
xsparse = np.full(n3, np.nan); xsparse[:7] = Xe[:7, 6]   # 7 complete rows: n <= k
wneg = we.copy(); wneg[7] = -0.5
cl_one = np.array(["one_cluster"] * n3, dtype=object)
late = rng3.uniform(40.0, 40.5, size=n3)        # ~55 km x 42 km box
lone = rng3.uniform(-90.5, -90.0, size=n3)
# Phase A.1 (drawn last, so every column above is unchanged): firms within four
# 'worlds' of six clusters each, so in (fan, fw) fan is nested in cl, fw is not,
# and the FE graph has four components; and a mediator that is an exact linear
# function of the treatment and x1.
fw = np.array([f"world_firm_{c:02d}" for c in (cl_c // 6) * 5 + rng3.integers(0, 5, size=n3)],
              dtype=object)
mcol = 2.0 * tE + 0.5 * Xe[:, 0]

ext = dict(y=yE, tvar=tE, mvar=mE, yneg=yneg, z1=z1e, z2=z2e)
for j in range(8):
    ext[f"x{j + 1}"] = Xe[:, j]
for j in range(20):
    ext[f"c{j + 1:02d}"] = Ce[:, j]
ext.update(w=we, fa=fa, fb=fb, cl=cl, fan=fan, fd1=fd1, fd2=fd2, fz=fz, wz=wz,
           cl_txt=cl_txt, fa_txt=fa_txt, dtxt=dtxt, mmiss=mmiss, xmiss=xmiss,
           year=year, year2=year2, x12=x12,
           ys=yE + 10.0, ts=tE + 5.0, ms=mE - 3.0, zs=z1e + 2.0)
for j in range(8):
    ext[f"xs{j + 1}"] = Xe[:, j] + shifts[j]
ext.update(zg=zg, tconst=tconst, z1dbl=z1dbl, xinf=xinf, xsparse=xsparse, wneg=wneg,
           cl_one=cl_one, lat=late, lon=lone)
ext.update(fw=fw, mcol=mcol)
write_csv(os.path.join(OUT, "data_ext.csv"), ext)

# ---------------- negative Conley variance (Phase A.1) ----------------
# 12 x 12 lattice, uniform kernel with cutoff 2.5 (not PSD in two dimensions).
# The outcome error is built in the kernel's negative eigenspace and orthogonal
# to the instruments, so the variance of the total effect is negative (the
# Stata reviewer's S-4(a) design): SE NaN, AR set 'unavailable'.
rng4 = np.random.default_rng(31)
gx = np.tile(np.arange(12.0), 12)
gy = np.repeat(np.arange(12.0), 12)
n4 = len(gx)
z4 = rng4.standard_normal(n4)
t4 = z4 + rng4.standard_normal(n4)
m4 = 0.5 * t4 + rng4.standard_normal(n4)
D4 = np.sqrt((gx[:, None] - gx[None, :]) ** 2 + (gy[:, None] - gy[None, :]) ** 2)
lam4, U4 = np.linalg.eigh((D4 <= 2.5).astype(float))
Zm4 = np.column_stack([np.ones(n4), z4])
Xm4 = np.column_stack([np.ones(n4), t4])
Xh4 = Zm4 @ np.linalg.solve(Zm4.T @ Zm4, Zm4.T @ Xm4)
h4 = Xh4 @ np.linalg.inv(Xh4.T @ Xh4)[:, 1]
Un4 = U4[:, lam4 < -1e-8]
c4 = np.linalg.svd(Zm4.T @ (Un4 / h4[:, None]))[2][-1]
e4 = (Un4 @ c4) / h4
e4 = e4 / np.sqrt(e4 @ e4)
write_csv(os.path.join(OUT, "data_negvar.csv"),
          dict(y=1.0 + 0.5 * t4 + e4, tvar=t4, mvar=m4, z1=z4, gx=gx, gy=gy))

# ---------------- weakly connected two-way FE ----------------
# A worker-firm chain: 2,000 workers x 2 periods and 200 firms; worker 10k moves to
# the next firm in period 2, so one mover per ten workers links neighboring firms
# (worker 1, alone in firm 1, forms a second component). Alternating projections
# converge slowly on such a graph; the exact two-way projection handles it.
rng5 = np.random.default_rng(20261017)
n5 = 4000
wkr = np.arange(n5) // 2 + 1                       # ceil(_n/2)
per = np.arange(n5) % 2 + 1
firm = (wkr - 1 + 9) // 10 + 1                     # ceil(200*(wkr-1)/2000) + 1
firm = firm + ((per == 2) & (wkr % 10 == 0))
firm = np.minimum(firm, 200)
z5 = rng5.standard_normal(n5)
e5 = rng5.standard_normal(n5)
t5 = 0.7 * z5 + e5 + firm / 50
m5 = 0.5 * t5 + 0.5 * e5 + rng5.standard_normal(n5)
y5 = 0.6 * m5 + 0.2 * t5 + wkr / 1000 + rng5.standard_normal(n5)
write_csv(os.path.join(OUT, "data_chain.csv"),
          dict(y=y5, tvar=t5, mvar=m5, z1=z5, wkr=wkr, per=per, firm=firm))

# ---------------- configuration grid ----------------
XB = ["x1", "x2", "x3"]
configs = [
    dict(name="base_robust", data="main", Z=["z1"], X=XB, se="robust"),
    dict(name="classical", data="main", Z=["z1"], X=XB, se="classical"),
    dict(name="cluster", data="main", Z=["z1"], X=XB, se="cluster", clusterid="clu"),
    dict(name="spatial_uni_hav", data="main", Z=["z1"], X=XB, se="spatial",
         lat="lat", lon="lon", spatial_cutoff=150.0, spatial_kernel="uniform", spatial_metric="haversine"),
    dict(name="spatial_bart_euc", data="main", Z=["z1"], X=XB, se="spatial",
         lat="lat", lon="lon", spatial_cutoff=1.2, spatial_kernel="bartlett", spatial_metric="euclidean"),
    dict(name="wt_aw_robust", data="main", Z=["z1"], X=XB, se="robust", weights="w", wtype="aweight"),
    dict(name="wt_pw_cluster", data="main", Z=["z1"], X=XB, se="cluster", clusterid="clu",
         weights="w", wtype="pweight"),
    dict(name="wt_pw_classical", data="main", Z=["z1"], X=XB, se="classical", weights="w", wtype="pweight"),
    dict(name="wt_aw_classical", data="main", Z=["z1"], X=XB, se="classical", weights="w", wtype="aweight"),
    dict(name="fe1_cluster", data="main", Z=["z1"], X=XB, se="cluster", clusterid="clu", fe="fe1"),
    dict(name="fe2way_robust", data="main", Z=["z1"], X=XB, se="robust", fe=["fe1", "fe2"]),
    dict(name="dummies_robust", data="main", Z=["z1"], X=XB, se="robust", dummies=["dcat", "dnum", "dflo"]),
    dict(name="glob_string_X", data="main", Z=["z1"], X="x1 x2 reg_*", se="robust"),
    dict(name="k2_cluster", data="main", Z=["z1", "z2"], X=XB, se="cluster", clusterid="clu",
         ar_grid_halfwidth=50.0, ar_grid_points=801),
    dict(name="weak_ar", data="weak", Z=["z1"], X=["x1"], se="robust"),
    dict(name="weak_ar_cluster", data="weak", Z=["z1"], X=["x1"], se="cluster", clusterid="clu"),
    dict(name="alpha10_tau20", data="main", Z=["z1"], X=XB, se="robust", alpha=0.10, op_tau=0.20),
    dict(name="op_offgrid", data="main", Z=["z1"], X=XB, se="robust", op_alpha=0.10),
    dict(name="noX_robust", data="main", Z=["z1"], se="robust"),
    dict(name="fe1_weights", data="main", Z=["z1"], X=XB, se="robust", fe="fe1", weights="w", wtype="aweight"),
    # ---- LASSO ----
    dict(name="lasso_base", data="main", Z=["z1"], X=XB, se="robust", lasso=True,
         dummies=["dnum"], show_selection=True),
    dict(name="lasso_keep_glob", data="main", Z=["z1"], X=XB, se="robust", lasso=True,
         dummies=["dnum"], lasso_keep="x1*",
         expand_x=["square", "cube", "interact_with_X"],
         expand_dummies=["interact_with_X", "interact_with_dummies"]),
    dict(name="lasso_override", data="main", Z=["z1"], X=XB, se="robust", lasso=True,
         lasso_lambda_override=25.0),
    dict(name="lasso_cluster_wt", data="main", Z=["z1"], X=XB, se="cluster", clusterid="clu",
         weights="w", wtype="aweight", lasso=True, dummies=["dnum"]),
    dict(name="lasso_spatial", data="main", Z=["z1"], X=XB, se="spatial",
         lat="lat", lon="lon", spatial_cutoff=150.0, spatial_kernel="uniform",
         spatial_metric="haversine", lasso=True),
    dict(name="lasso_expand_override", data="main", Z=["z1"], X=XB, se="robust", lasso=True,
         dummies=["dnum"], expand_override={"x1": ["square"], "_D_dnum_2": ["square"]}),
    dict(name="lasso_fe", data="main", Z=["z1"], X=XB, se="robust", fe="fe1", lasso=True),
    dict(name="lasso_k2", data="main", Z=["z1", "z2"], X=XB, se="robust", lasso=True,
         ar_grid_halfwidth=50.0, ar_grid_points=801),
    dict(name="lasso_dcat_probe", data="main", Z=["z1"], X=XB, se="robust", lasso=True,
         dummies=["dcat"]),
    dict(name="lasso_dcat_keep", data="main", Z=["z1"], X=["x1"], se="robust", lasso=True,
         dummies=["dcat", "dnum"], lasso_keep="_D_dcat_*",
         expand_dummies=["interact_with_X", "interact_with_dummies"],
         show_selection=True),
    # NOTE: a cutoff >> data extent (e.g., 5000 km uniform) degenerates the
    # Conley meat to (sum of scores) outer-product == 0 by the normal
    # equations. Since 2026-10 this is an error (uniform kernel; see
    # err_spatial_allpairs below).
    dict(name="glob_class_X", data="main", Z=["z1"], X="x[12] x3", se="robust"),
]

# ---- 2026-10 (WP8): new scenarios on data_ext.csv ----
XE = ["x1", "x2", "x3"]
X8 = [f"x{j}" for j in range(1, 9)]
configs += [
    dict(name="neg_theta_robust", data="ext", Y="yneg", Z=["z1"], X=XE, se="robust"),
    dict(name="fe1_wt_aw_robust", data="ext", Z=["z1"], X=XE, se="robust", fe="fa",
         weights="w", wtype="aweight"),
    dict(name="fe1_wt_pw_cluster", data="ext", Z=["z1"], X=XE, se="cluster", clusterid="cl", fe="fa",
         weights="w", wtype="pweight"),
    dict(name="fe2way_wt_aw_robust", data="ext", Z=["z1"], X=XE, se="robust", fe=["fa", "fb"],
         weights="w", wtype="aweight"),
    dict(name="fe2way_wt_pw_cluster", data="ext", Z=["z1"], X=XE, se="cluster", clusterid="cl",
         fe=["fa", "fb"], weights="w", wtype="pweight"),
    dict(name="fe_nested_cluster", data="ext", Z=["z1"], X=XE, se="cluster", clusterid="cl", fe="fan"),
    dict(name="fe2way_nested_cluster", data="ext", Z=["z1"], X=XE, se="cluster", clusterid="cl",
         fe=["fan", "fb"]),
    dict(name="fe2way_disconnected", data="ext", Z=["z1"], X=XE, se="robust", fe=["fd1", "fd2"]),
    dict(name="fe2way_disconnected_cluster", data="ext", Z=["z1"], X=XE, se="cluster", clusterid="cl",
         fe=["fd1", "fd2"]),
    dict(name="k2_het_robust", data="ext", Z=["z1", "z2"], X=XE, se="robust"),
    dict(name="k2_het_cluster", data="ext", Z=["z1", "z2"], X=XE, se="cluster", clusterid="cl"),
    # LASSO without FE must not depend on the location of Y, T, M, Z (same selection and
    # estimates in each pair); shifted controls change their squares, so the second
    # pair has no expansion
    dict(name="lasso_shift_raw", data="ext", Z=["z1"], X=X8, se="robust", lasso=True),
    dict(name="lasso_shift_ytmz", data="ext", Y="ys", T="ts", M="ms", Z=["zs"], X=X8,
         se="robust", lasso=True),
    dict(name="lasso_noexp_raw", data="ext", Z=["z1"], X=X8, se="robust", lasso=True, expand_x=[]),
    dict(name="lasso_noexp_shifted", data="ext", Y="ys", T="ts", M="ms", Z=["zs"],
         X=[f"xs{j}" for j in range(1, 9)], se="robust", lasso=True, expand_x=[]),
    dict(name="lasso_wt_pweight", data="ext", Z=["z1"], X=X8, se="robust", weights="w",
         wtype="pweight", lasso=True),
    dict(name="lasso_cluster_only", data="ext", Z=["z1"], X=X8, se="cluster", clusterid="cl",
         lasso=True),
    dict(name="lasso_postlasso_loadings", data="ext", Z=["z1"],
         X=X8 + [f"c{j:02d}" for j in range(1, 21)], se="robust", lasso=True, show_selection=True),
    dict(name="fe1_zero_weights", data="ext", Z=["z1"], X=XE, se="robust", fe="fz",
         weights="wz", wtype="aweight"),
    dict(name="missing_text_ids", data="ext", Z=["z1"], X=XE, se="cluster", clusterid="cl_txt",
         fe="fa_txt", dummies=["dtxt"]),
    dict(name="illcond_controls", data="ext", Z=["z1"], X=["x1", "year", "year2"], se="robust"),
    dict(name="collinear_control", data="ext", Z=["z1"], X=["x1", "x2", "x12"], se="robust"),
    # Phase A.1
    dict(name="fe2way_nested_disconnected", data="ext", Z=["z1"], X=XE, se="cluster",
         clusterid="cl", fe=["fan", "fw"]),
    dict(name="spatial_negvar", data="negvar", Z=["z1"], se="spatial", lat="gx", lon="gy",
         spatial_cutoff=2.5, spatial_kernel="uniform", spatial_metric="euclidean"),
    # Addendum A.2 (O-2): two-way FE on a weakly connected worker-firm chain
    dict(name="fe2_chain", data="chain", Z=["z1"], se="robust", fe=["wkr", "firm"], ar_ci=False),
]

boot_configs = [
    dict(name="boot_pairs_basic", data="main", Z=["z1"], X=XB, se="robust",
         method="pairs", B=60, seed=555, ci_type="basic"),
    dict(name="boot_block_pct", data="main", Z=["z1"], X=XB, se="cluster", clusterid="clu",
         method="block", B=60, seed=777, ci_type="percentile"),
    dict(name="boot_block_stud", data="main", Z=["z1"], X=XB, se="cluster", clusterid="clu",
         method="block", B=60, seed=888, ci_type="studentized", bootstrap_diagnostics=True),
    # 2026-10: missing values in a role variable, a control and the cluster id
    dict(name="boot_pairs_missing", data="ext", M="mmiss", Z=["z1"], X=["x1", "x2", "xmiss"],
         se="robust", method="pairs", B=60, seed=4242, ci_type="percentile"),
    dict(name="boot_block_missing", data="ext", M="mmiss", Z=["z1"], X=["x1", "x2", "xmiss"],
         se="cluster", clusterid="cl_txt", method="block", B=60, seed=4343, ci_type="studentized"),
    # Phase A.1: FE nested in the clusters get one level per drawn copy
    dict(name="boot_block_fe_nested", data="ext", Z=["z1"], X=XE, fe="fan", se="cluster",
         clusterid="cl", method="block", B=60, seed=4545, ci_type="studentized"),
]

# ---- 2026-10 (WP8): degenerate inputs; every language must stop with an error
# of the class in "expect" (classes: see error_classes below) ----
error_configs = [
    dict(name="err_one_cluster", data="ext", expect="clusters", Z=["z1"], X=XE, se="cluster",
         clusterid="cl_one"),
    dict(name="err_z_absorbed_by_fe", data="ext", expect="identification", Z=["zg"], X=XE,
         se="robust", fe="fa"),
    dict(name="err_z_collinear", data="ext", expect="identification", Z=["z1", "z1dbl"], X=XE,
         se="robust"),
    dict(name="err_const_treatment", data="ext", expect="identification", T="tconst", Z=["z1"],
         X=XE, se="robust"),
    dict(name="err_x_contains_y", data="ext", expect="role", Z=["z1"], X=["x1", "y"], se="robust"),
    dict(name="err_z_equals_t", data="ext", expect="role", Z=["tvar"], X=XE, se="robust"),
    dict(name="err_dup_instrument", data="ext", expect="role", Z=["z1", "z1"], X=XE, se="robust"),
    dict(name="err_inf_value", data="ext", expect="nonfinite", Z=["z1"], X=["x1", "xinf"],
         se="robust"),
    dict(name="err_negative_weight", data="ext", expect="weights", Z=["z1"], X=XE, se="robust",
         weights="wneg", wtype="aweight"),
    dict(name="err_n_le_k", data="ext", expect="too_few_obs", Z=["z1"],
         X=["x1", "x2", "x3", "xsparse"], se="robust"),
    dict(name="err_spatial_allpairs", data="ext", expect="spatial_degenerate", Z=["z1"], X=XE,
         se="spatial", lat="lat", lon="lon", spatial_cutoff=500.0, spatial_kernel="uniform",
         spatial_metric="haversine"),
    dict(name="err_spatial_args_robust", data="ext", expect="option", Z=["z1"], X=XE, se="robust",
         lat="lat", lon="lon"),
    # Phase A.1
    dict(name="err_m_collinear_t", data="ext", expect="identification", M="mcol", Z=["z1"],
         X=["x1"], se="robust"),
]

# Error class = the first class with a phrase contained in the lower-cased error
# message; "other" when none matches. The phrases follow the Python messages.
error_classes = [
    ["clusters", ["at least 2 clusters"]],
    ["spatial_degenerate", ["covers all pairs"]],
    ["nonfinite", ["infinite"]],
    ["weights", ["weights must be non-negative", "all weights are zero"]],
    ["too_few_obs", ["not enough observations", "no observations remaining"]],
    ["identification", ["no variation", "collinear after", "collinear with the treatment",
                        "not identified"]],
    ["role", ["must not contain the outcome", "must differ from the outcome",
              "three different columns", "duplicate instrument"]],
    ["option", ["used only with se=", "must be true or false", "must lie strictly between",
                "must be one of"]],
]

with open(os.path.join(OUT, "configs.json"), "w") as f:
    json.dump(dict(configs=configs, boot_configs=boot_configs, error_configs=error_configs,
                   error_classes=error_classes), f, indent=1)

# ---------------- fnmatch ground truth for glob_check.R ----------------
import fnmatch
gnames = ['x1','x2','x3','x22','reg_a','reg_b','y1','x[','a.b1','aXb1','x1^2','x1^23','x1 Q x2']
gpats = ['x*','x?','x[12]','x[!3]','reg_*','*_a','x[1-2]','[xy]1','x[','a.b*','x[^12]','*','??','x1^2*','x1^2','?[12]']
with open(os.path.join(OUT, "fn_expected.json"), "w") as f:
    json.dump({p: sorted(fnmatch.filter(gnames, p)) for p in gpats}, f)

print("wrote data_main.csv, data_weak.csv, data_ext.csv, data_negvar.csv, data_chain.csv, configs.json, "
      "fn_expected.json")
