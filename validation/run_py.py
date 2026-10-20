"""Run the Python package ivmediate across the config grid; flatten results to JSON.

The package is imported from python/src of this repository (not from an installed
copy), so the harness always tests the working tree.

configs.json holds three lists, each run against data_<data>.csv:
  configs        analysis configs (ivmediate.ivmediate + ivmediate.print_results)
  boot_configs   bootstrap configs (ivmediate.bootstrap + ivmediate.print_bootstrap)
  error_configs  degenerate inputs that must stop with an error of class "expect"
A config may name Y, T and M (default y, tvar, mvar); every other field except
name, data and expect is passed to the function.

Each config's output is a flat dict of strings. Metadata keys start with "_":
  _status        "ok" or "error"
  _error_class   on errors: the first class in configs.json "error_classes" with a
                 phrase contained in the lower-cased message, else "other"
  _ERROR         on errors: "<ExceptionType>: <message>"
  _stdout_lines  number of printed lines (informational)
Bootstrap configs also export the draws so other languages can replay them:
  bootidx_<name>.csv       B rows of 0-based draws: rows of the estimation sample
                           (pairs) or clusters in order of first appearance in the
                           estimation sample (block)
  bootidx_<name>_rows.csv  the estimation sample: 0-based row positions in
                           data_<data>.csv (header excluded), one per line
"""
import io, json, os, sys, contextlib
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # error messages contain Greek letters
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "python", "src")))
import ivmediate
from ivmediate import _core   # private helpers: the static OP table

with open(os.path.join(HERE, "configs.json")) as f:
    CFG = json.load(f)

# round_trip: parse every float exactly, as the other languages do (the default
# pandas parser can be 1 ulp off)
DATASETS = {nm: pd.read_csv(os.path.join(HERE, f"data_{nm}.csv"), float_precision="round_trip")
            for nm in ("main", "weak", "ext", "negvar", "chain")}
ROLE_DEFAULTS = dict(Y="y", T="tvar", M="mvar")
ERROR_CLASSES = CFG.get("error_classes", [])


def error_class(msg):
    low = str(msg).lower()
    for cls, phrases in ERROR_CLASSES:
        if any(ph in low for ph in phrases):
            return cls
    return "other"


def fnum(v):
    v = float(v)
    if np.isnan(v):
        return "NaN"
    return "%.17g" % v


def put(out, path, v):
    if v is None:
        out[path] = "NULL"
    elif isinstance(v, (bool, np.bool_)):
        out[path] = "TRUE" if v else "FALSE"
    elif isinstance(v, (int, np.integer)):
        out[path] = str(int(v))
    elif isinstance(v, (float, np.floating)):
        out[path] = fnum(v)
    else:
        out[path] = str(v)


def put_vec(out, path, v):
    """Scalars as `path`, lists (K >= 2) as `path.1`, `path.2`, ..."""
    if isinstance(v, (list, tuple, np.ndarray)):
        for i, vi in enumerate(v, start=1):
            put(out, f"{path}.{i}", vi)
    else:
        put(out, path, v)


def flat_ar(out, pre, ar):
    if ar is None:
        put(out, pre, "NULL")
        return
    for k in ("shape", "method", "contains_zero", "cv", "df_chi2", "alpha",
              "stat_at_zero", "pval_at_zero", "lower", "upper",
              "left_upper", "right_lower"):
        put(out, f"{pre}.{k}", ar[k])
    segs = ar.get("segments")
    if segs is not None:
        put(out, f"{pre}.n_segments", len(segs))
        for i, (lo, hi, ol, orr) in enumerate(segs, start=1):
            put(out, f"{pre}.seg{i}.lo", lo)
            put(out, f"{pre}.seg{i}.hi", hi)
            put(out, f"{pre}.seg{i}.open_left", ol)
            put(out, f"{pre}.seg{i}.open_right", orr)


def flat_eps(out, pre, e):
    e = np.asarray(e, dtype=float)
    put(out, f"{pre}.mean", float(np.mean(e)))
    put(out, f"{pre}.sd", float(np.sqrt(np.mean((e - e.mean()) ** 2))))
    for i in range(3):
        put(out, f"{pre}.v{i+1}", float(e[i]))


def flat_lasso_detail(out, pre, d):
    if d is None:
        put(out, pre, "NULL")
        return
    put(out, f"{pre}.n_selected", d["n_selected"])
    put(out, f"{pre}.lambda_value", d["lambda_value"])
    put(out, f"{pre}.converged", d["converged"])
    put(out, f"{pre}.n_iter", d["n_iter"])
    put(out, f"{pre}.method", d["method"])
    put(out, f"{pre}.names", "|".join(d["names"]))
    put(out, f"{pre}.coefs", "|".join("%.10g" % c for c in d["coefs"]))
    put(out, f"{pre}.indices", "|".join(str(int(i) + 1) for i in d["indices"]))  # 1-based


def flatten_result(r):
    out = {}
    for k in ("n", "n_instruments", "n_covariates", "n_dummies", "se_type", "alpha",
              "n_fe_groups", "n_dropped_singletons", "df_fe_absorbed", "weighted", "wtype",
              "spatial_cutoff", "spatial_kernel", "spatial_metric",
              "theta_hat", "se_theta", "tau_hat", "se_tau", "gamma_hat", "se_gamma",
              "total_hat", "se_total", "indirect", "se_indirect", "cov_gamma_theta",
              "decomp_check", "decomp_error",
              "F_TZ", "F_MZ_T", "F_TZ_op", "F_MZ_T_op",
              "F_TZ_op_critval", "F_MZ_T_op_critval", "op_tau", "op_alpha",
              "c_TM_hat", "rho_TM_hat", "sigma_eT", "sigma_eY",
              "kappa_star", "kappa_star_CI", "rho_TY_star",
              "f_reduced_form", "f_first_stage_cond",
              "RV_theta", "XRV_theta", "rv_rf", "rv_fs", "xrv_rf", "xrv_fs",
              "lasso", "lasso_penalty",
              # 2026-10 (WP8)
              "n_clusters", "K_eff_TZ", "K_eff_MZ_T",
              "kappa_star_CI_informative", "F_MZ_T_above_crit"):
        put(out, k, r.get(k))
    put_vec(out, "pi1", r.get("pi1"))
    put_vec(out, "se_pi1", r.get("se_pi1"))
    nested = r.get("fe_nested_in_cluster")
    put(out, "fe_nested_in_cluster",
        "NULL" if nested is None else ",".join("TRUE" if b else "FALSE" for b in nested))
    fe = r.get("fe")
    put(out, "fe", ",".join(fe) if fe else "NULL")
    for nm in ("ci_theta", "ci_tau", "ci_gamma", "ci_total", "ci_indirect"):
        ci = r[nm]
        put(out, nm + ".1", ci[0]); put(out, nm + ".2", ci[1])
    for i, v in enumerate(r["F_MZ_T_individual"], start=1):
        put(out, f"F_MZ_T_individual.{i}", v)
    for i, row in enumerate(r["kappa_adj_table"], start=1):
        for k in ("kappa", "kappa_frac", "theta_adj", "ci_lo", "ci_hi"):
            put(out, f"kappa_adj.{i}.{k}", row[k])
    for nm in ("ar_theta", "ar_gamma", "ar_total"):
        flat_ar(out, nm, r.get(nm))
    for nm in ("eps_T_hat", "eps_M_hat", "eps_Y_hat"):
        flat_eps(out, nm, r[nm])
    li = r.get("lasso_info")
    if li is not None:
        put(out, "lasso_info.n_selected_union", li["n_selected_union"])
        put(out, "lasso_info.n_total", li["n_total"])
        put(out, "lasso_info.union_names", "|".join(li["union_names"]))
        det = li.get("details") or {}
        flat_lasso_detail(out, "lasso_info.Y", det.get("Y"))
        flat_lasso_detail(out, "lasso_info.T", det.get("T"))
        flat_lasso_detail(out, "lasso_info.M", det.get("M"))
        for k, dz in enumerate(det.get("Z") or [], start=1):
            flat_lasso_detail(out, f"lasso_info.Z{k}", dz)
    return out


def flatten_boot(bt):
    out = {}
    for k in ("n_failed", "method", "B", "ci_type", "alpha", "n_est"):
        put(out, k, bt[k])
    for p, v in bt["point"].items():
        put(out, f"point.{p}", v)
    for p, v in bt["boot_se"].items():
        put(out, f"boot_se.{p}", v)
    for p, ci in bt["boot_ci"].items():
        put(out, f"boot_ci.{p}.1", ci[0]); put(out, f"boot_ci.{p}.2", ci[1])
    for p, v in bt["analytic_se"].items():
        put(out, f"analytic_se.{p}", v)
    bd = bt.get("boot_dist") or {}
    for p, vec in bd.items():
        vec = np.asarray(vec, dtype=float)
        ok = vec[~np.isnan(vec)]
        put(out, f"boot_dist.{p}.mean", float(np.mean(ok)))
        put(out, f"boot_dist.{p}.sd", float(np.std(ok, ddof=1)))
        for i in range(3):
            put(out, f"boot_dist.{p}.v{i+1}", float(vec[i]))
    return out


def split_cfg(cfg):
    """-> name, data, role names (Y, T, M), function kwargs, expected error class."""
    cfg = dict(cfg)
    name = cfg.pop("name")
    data = DATASETS[cfg.pop("data")]
    expect = cfg.pop("expect", None)
    roles = {k: cfg.pop(k, v) for k, v in ROLE_DEFAULTS.items()}
    if "expand_x" in cfg:
        cfg["expand_x"] = tuple(cfg["expand_x"])
    if "expand_dummies" in cfg:
        cfg["expand_dummies"] = tuple(cfg["expand_dummies"])
    if "expand_override" in cfg:
        cfg["expand_override"] = {k: tuple(v) for k, v in cfg["expand_override"].items()}
    return name, data, roles, cfg, expect


def error_entry(e):
    return {"_status": "error", "_error_class": error_class(e), "_ERROR": f"{type(e).__name__}: {e}"}


ALL = {}
n_unexpected = 0

# ---- static OP table ----
st = {}
put(st, "version", ivmediate.__version__)
for K, row in _core._OP_CRIT_5PCT_WALD.items():
    for tau, v in row.items():
        put(st, f"op.{K}.{tau:.2f}", v)
for K, tau, al in [(1, 0.10, 0.05), (5, 0.30, 0.05), (30, 0.05, 0.05),
                   (31, 0.10, 0.05), (1, 0.15, 0.05), (1, 0.10, 0.01)]:
    put(st, f"opcv.{K}.{tau}.{al}", _core._op_critical_value(K, tau, al))
ALL["_static"] = st

# ---- analysis configs and degenerate inputs ----
for cfg in CFG["configs"] + CFG.get("error_configs", []):
    name, data, roles, kw, expect = split_cfg(cfg)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            r = ivmediate.ivmediate(data, **roles, **kw)
            # smoke test the printer, optional sections included (Addendum A.2)
            ivmediate.print_results(r, title=name, kappa_table=True, r2_sensitivity=True)
        out = flatten_result(r)
        out["_status"] = "ok"
        out["_stdout_lines"] = str(len(buf.getvalue().splitlines()))
        msg = "ok"
    except Exception as e:
        out = error_entry(e)
        msg = f"ERROR [{out['_error_class']}] {out['_ERROR']}"
    got = out.get("_error_class", "ok")
    if got != (expect or "ok"):
        n_unexpected += 1
        msg += f"   <-- UNEXPECTED (expected {expect or 'ok'})"
    ALL[name] = out
    print(f"[py] {name}: {msg}")

# ---- bootstrap configs (also export the draws for the other languages) ----
for bcfg in CFG["boot_configs"]:
    name, data, roles, kw, expect = split_cfg(bcfg)
    data = data.reset_index(drop=True)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            bt = ivmediate.bootstrap(data, **roles, return_boot_dist=True, progress=False,
                                     _return_indices=True, **kw)
            ivmediate.print_bootstrap(bt, title=name)
        np.savetxt(os.path.join(HERE, f"bootidx_{name}.csv"), bt["_indices"], fmt="%d", delimiter=",")
        np.savetxt(os.path.join(HERE, f"bootidx_{name}_rows.csv"), bt["_rows"], fmt="%d")
        out = flatten_boot(bt)
        out["_status"] = "ok"
        msg = f"ok (n_est={bt['n_est']}, n_failed={bt['n_failed']})"
    except Exception as e:
        out = error_entry(e)
        msg = f"ERROR [{out['_error_class']}] {out['_ERROR']}"
    got = out.get("_error_class", "ok")
    if got != (expect or "ok"):
        n_unexpected += 1
        msg += f"   <-- UNEXPECTED (expected {expect or 'ok'})"
    ALL[name] = out
    print(f"[py] {name}: {msg}")

with open(os.path.join(HERE, "out_py.json"), "w") as f:
    json.dump(ALL, f, indent=0, sort_keys=True)
print(f"wrote out_py.json ({len(ALL) - 1} configs; {n_unexpected} with an unexpected status)")
sys.exit(1 if n_unexpected else 0)
