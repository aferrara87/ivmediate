// =============================================================================
//  Stata runner for the cross-language harness (2026-10 harness, ivmediate 2.0.0)
//
//  Runs stata/ivmediate.ado on every config of validation/configs.json and
//  writes validation/out_stata.json in the harness format (validation/README.md):
//  one flat object of strings per config, numbers %.17g, NaN, NULL, TRUE/FALSE,
//  metadata _status, _error_class, _ERROR, _stdout_lines.
//
//  Run from validation/stata:
//      "C:\Program Files\StataNow19\StataSE-64.exe" /e do run_stata.do
//  then compare:  python -B ../compare.py stata
//
//  Notes
//  - The configs are translated by hand into Stata syntax (below); the Python
//    names are kept as config names. Y/T/M overrides, data_ext.csv and
//    data_negvar.csv are used as in configs.json.
//  - Error classes come from the phrase table of configs.json (error_classes),
//    matched against the error message the ado prints (read back from a log).
//    One Stata phrase is added: "used only with vce(conley)" -> option.
//  - err_inf_value is not expressible in Stata (no infinity); it is written as
//    _status "n/a". lasso_expand_override uses expand_override, which the Stata
//    command does not have (expandx()/expanddum() apply to every control); it is
//    also written as "n/a". glob_class_X uses a character class (x[12]), which a
//    Stata varlist cannot express; it runs with the resolved list x1 x2 x3.
//  - LASSO term names are translated to the Python names: a#b -> "a × b",
//    k.v -> _D_v_k (value labels for the encoded string variables dcat, dtxt).
//  - Every run turns on kappatable and r2sensitivity (printer smoke test, A.2).
//  - Bootstrap configs replay the Python draws: bootidx_<name>.csv (+1, since
//    ivmediate's bootidx() is 1-based) from validation/reference_2026-10/. The
//    estimation sample is checked against bootidx_<name>_rows.csv.
// =============================================================================
clear all
set more off
set linesize 255
version 16
adopath ++ "../../stata"
capture program drop ivmediate
// ivmediate must resolve to stata/ivmediate.ado of this repository
which ivmediate
quietly findfile ivmediate.ado
if `"`r(fn)'"' != "../../stata/ivmediate.ado" {
    display as error "ivmediate must resolve to ../../stata/ivmediate.ado"
    exit 601
}
local REF "../reference_2026-10"

// ---- datasets ----------------------------------------------------------------
import delimited using "../data_main.csv", varnames(1) asdouble encoding(utf-8) clear
encode dcat, gen(dcatn)
gen long _row0 = _n - 1
save _main.dta, replace
import delimited using "../data_weak.csv", varnames(1) asdouble encoding(utf-8) clear
gen long _row0 = _n - 1
save _weak.dta, replace
import delimited using "../data_ext.csv", varnames(1) asdouble encoding(utf-8) clear
replace dtxt = strtrim(dtxt)
encode dtxt, gen(dtxtn)
gen long _row0 = _n - 1
save _ext.dta, replace
import delimited using "../data_negvar.csv", varnames(1) asdouble encoding(utf-8) clear
gen long _row0 = _n - 1
save _negvar.dta, replace
import delimited using "../data_chain.csv", varnames(1) asdouble encoding(utf-8) clear
gen long _row0 = _n - 1
save _chain.dta, replace

// ---- Mata: JSON writer and result dumpers --------------------------------------
mata:
mata set matastrict off

void rj_open(string scalar path)
{
    external real scalar _RJ_fh, _RJ_n
    unlink(path)
    _RJ_fh = fopen(path, "w")
    fput(_RJ_fh, "{")
    _RJ_n = 0
}

void rj_close()
{
    external real scalar _RJ_fh
    fput(_RJ_fh, "}")
    fclose(_RJ_fh)
}

string scalar rj_esc(string scalar s)
{
    string scalar o
    o = subinstr(s, char(92), char(92) + char(92))
    o = subinstr(o, char(34), char(92) + char(34))
    o = subinstr(o, char(9), " ")
    o = subinstr(o, char(10), " ")
    o = subinstr(o, char(13), " ")
    return(o)
}

void rj_begin(string scalar cfg)
{
    external string scalar _RJ_cfg
    external string matrix _RJ_buf
    _RJ_cfg = cfg
    _RJ_buf = J(0, 2, "")
}

void rj_put(string scalar k, string scalar v)
{
    external string matrix _RJ_buf
    _RJ_buf = _RJ_buf \ (k, v)
}

void rj_end()
{
    external real scalar   _RJ_fh, _RJ_n
    external string scalar _RJ_cfg
    external string matrix _RJ_buf
    real scalar i
    string scalar line
    fput(_RJ_fh, (_RJ_n > 0 ? "," : "") + char(34) + rj_esc(_RJ_cfg) + char(34) + ": {")
    for (i = 1; i <= rows(_RJ_buf); i++) {
        line = char(34) + rj_esc(_RJ_buf[i, 1]) + char(34) + ": " + char(34) + rj_esc(_RJ_buf[i, 2]) + char(34)
        if (i < rows(_RJ_buf)) line = line + ","
        fput(_RJ_fh, line)
    }
    fput(_RJ_fh, "}")
    _RJ_n = _RJ_n + 1
}

string scalar rj_num(real scalar x)
{
    if (x >= .) return("NaN")
    return(strtrim(sprintf("%23.17g", x)))
}

string scalar rj_numN(real scalar x)
{
    if (x >= .) return("NULL")
    return(strtrim(sprintf("%23.17g", x)))
}

string scalar rj_bool(real scalar x)
{
    if (x >= .) return("NULL")
    return(x != 0 ? "TRUE" : "FALSE")
}

real scalar en(string scalar nm)
{
    real matrix v
    v = st_numscalar("e(" + nm + ")")
    if (rows(v) == 0) return(.)
    return(v[1, 1])
}

// Python's "%.10g"
string scalar rj_g10(real scalar x)
{
    string scalar s, f, ex
    real scalar   e, pos
    if (x >= .) return("nan")
    if (x == 0) return("0")
    s = strtrim(sprintf("%30.9e", x))
    pos = strpos(s, "e")
    e = strtoreal(substr(s, pos + 1, .))
    if (e >= -4 & e < 10) {
        f = strtrim(sprintf("%40." + strofreal(9 - e) + "f", x))
        if (strpos(f, ".")) {
            while (substr(f, strlen(f), 1) == "0") f = substr(f, 1, strlen(f) - 1)
            if (substr(f, strlen(f), 1) == ".") f = substr(f, 1, strlen(f) - 1)
        }
        return(f)
    }
    f = substr(s, 1, pos - 1)
    if (strpos(f, ".")) {
        while (substr(f, strlen(f), 1) == "0") f = substr(f, 1, strlen(f) - 1)
        if (substr(f, strlen(f), 1) == ".") f = substr(f, 1, strlen(f) - 1)
    }
    ex = strofreal(abs(e))
    if (strlen(ex) < 2) ex = "0" + ex
    return(f + "e" + (e < 0 ? "-" : "+") + ex)
}

string scalar rj_join(string rowvector s, string scalar sep)
{
    string scalar o
    real scalar   i
    o = ""
    for (i = 1; i <= cols(s); i++) o = o + (i > 1 ? sep : "") + s[i]
    return(o)
}

// Stata expansion-term name -> Python name
string scalar rj_pyname1(string scalar s)
{
    string scalar lev, v, src, lab
    real scalar   p, val
    p = strpos(s, ".")
    if (p == 0 | !regexm(s, "^[0-9]+\.")) return(s)
    lev = substr(s, 1, p - 1)
    v   = substr(s, p + 1, .)
    val = strtoreal(lev)
    src = v
    lab = lev
    if (v == "dcatn" | v == "dtxtn") {
        src = substr(v, 1, strlen(v) - 1)
        lab = st_vlmap(v, val)
    }
    return("_D_" + src + "_" + lab)
}

string scalar rj_pyname(string scalar s)
{
    string rowvector parts
    string scalar    o
    real scalar      i
    parts = tokens(subinstr(s, "#", " "))
    o = ""
    for (i = 1; i <= cols(parts); i++) o = o + (i > 1 ? " " + uchar(215) + " " : "") + rj_pyname1(parts[i])
    return(o)
}

void rj_ar(string scalar x)
{
    string scalar pre, shp, mth
    real matrix   S
    real scalar   i, ns
    pre = "ar_" + x
    shp = st_global("e(ar_" + x + "_shape)")
    mth = st_global("e(ar_" + x + "_method)")
    rj_put(pre + ".shape", shp)
    rj_put(pre + ".method", mth)
    rj_put(pre + ".contains_zero", rj_bool(en("ar_" + x + "_contains0")))
    rj_put(pre + ".cv", rj_num(en("ar_cv")))
    rj_put(pre + ".df_chi2", rj_num(en("ar_df")))
    rj_put(pre + ".alpha", rj_num((100 - en("ar_level")) / 100))
    rj_put(pre + ".stat_at_zero", rj_num(en("ar_" + x + "_stat0")))
    rj_put(pre + ".pval_at_zero", rj_num(en("ar_" + x + "_p")))
    rj_put(pre + ".lower", rj_num(en("ar_" + x + "_lower")))
    rj_put(pre + ".upper", rj_num(en("ar_" + x + "_upper")))
    rj_put(pre + ".left_upper", rj_num(en("ar_" + x + "_left_upper")))
    rj_put(pre + ".right_lower", rj_num(en("ar_" + x + "_right_lower")))
    if (mth == "grid" & shp != "unavailable") {
        S = J(0, 4, .)
        if (shp == "grid") S = st_matrix("e(ar_" + x + "_seg)")
        ns = rows(S)
        rj_put(pre + ".n_segments", strofreal(ns))
        for (i = 1; i <= ns; i++) {
            rj_put(pre + ".seg" + strofreal(i) + ".lo", rj_num(S[i, 1]))
            rj_put(pre + ".seg" + strofreal(i) + ".hi", rj_num(S[i, 2]))
            rj_put(pre + ".seg" + strofreal(i) + ".open_left", rj_bool(S[i, 3]))
            rj_put(pre + ".seg" + strofreal(i) + ".open_right", rj_bool(S[i, 4]))
        }
    }
}

void rj_eps(string scalar key, string scalar var)
{
    real colvector e
    real scalar    mu
    e = st_data(., var, "_es")
    mu = mean(e)
    rj_put(key + ".mean", rj_num(mu))
    rj_put(key + ".sd", rj_num(sqrt(mean((e :- mu) :^ 2))))
    rj_put(key + ".v1", rj_num(e[1]))
    rj_put(key + ".v2", rj_num(e[2]))
    rj_put(key + ".v3", rj_num(e[3]))
}

void rj_lasso_detail(string scalar pre, real scalar j, string rowvector pyn)
{
    real matrix      S, B, L, I, C
    real colvector   sel
    string rowvector meth, cf, ix
    real scalar      i
    S = st_matrix("e(lasso_S)")
    B = st_matrix("e(lasso_B)")
    L = st_matrix("e(lasso_lambda)")
    I = st_matrix("e(lasso_n_iter)")
    C = st_matrix("e(lasso_converged)")
    meth = tokens(st_global("e(lasso_methods)"))
    sel = selectindex(S[., j])
    rj_put(pre + ".n_selected", strofreal(rows(sel)))
    rj_put(pre + ".lambda_value", rj_num(L[j]))
    rj_put(pre + ".converged", rj_bool(C[j]))
    rj_put(pre + ".n_iter", strofreal(I[j]))
    rj_put(pre + ".method", meth[j])
    cf = J(1, rows(sel), "")
    ix = J(1, rows(sel), "")
    for (i = 1; i <= rows(sel); i++) {
        cf[i] = rj_g10(B[sel[i], j])
        ix[i] = strofreal(sel[i])
    }
    if (rows(sel)) rj_put(pre + ".names", rj_join(pyn[sel'], "|"))
    else           rj_put(pre + ".names", "")
    rj_put(pre + ".coefs", rj_join(cf, "|"))
    rj_put(pre + ".indices", rj_join(ix, "|"))
}

void rj_dumpres()
{
    real scalar      K, z, i, j, lasso
    real matrix      M, A
    string scalar    wt, ab, vce, s
    string rowvector nm, pyn, eqs
    string rowvector est, ses, cis

    K = en("k_inst")
    vce = st_global("e(se_type)")
    rj_put("n", rj_num(en("N")))
    rj_put("n_instruments", rj_num(K))
    rj_put("n_covariates", rj_num(en("k_cov")))
    rj_put("n_dummies", rj_num(en("k_dummies")))
    rj_put("se_type", vce)
    rj_put("alpha", rj_num((100 - en("level")) / 100))
    rj_put("n_fe_groups", rj_num(en("N_fegroups")))
    rj_put("n_dropped_singletons", rj_num(en("N_singleton")))
    rj_put("df_fe_absorbed", rj_num(en("df_fe")))
    wt = st_global("e(wtype)")
    rj_put("weighted", (wt != "" ? "TRUE" : "FALSE"))
    rj_put("wtype", (wt != "" ? wt : "NULL"))
    if (vce == "spatial") {
        rj_put("spatial_cutoff", rj_num(en("spatial_cutoff")))
        rj_put("spatial_kernel", st_global("e(spatial_kernel)"))
        rj_put("spatial_metric", st_global("e(spatial_metric)"))
    }
    else {
        rj_put("spatial_cutoff", "NULL")
        rj_put("spatial_kernel", "NULL")
        rj_put("spatial_metric", "NULL")
    }
    rj_put("theta_hat", rj_num(en("theta")));     rj_put("se_theta", rj_num(en("se_theta")))
    rj_put("tau_hat", rj_num(en("tau")));         rj_put("se_tau", rj_num(en("se_tau")))
    rj_put("gamma_hat", rj_num(en("gamma")));     rj_put("se_gamma", rj_num(en("se_gamma")))
    rj_put("total_hat", rj_num(en("total")));     rj_put("se_total", rj_num(en("se_total")))
    rj_put("indirect", rj_num(en("indirect")));   rj_put("se_indirect", rj_num(en("se_indirect")))
    rj_put("cov_gamma_theta", rj_num(en("cov_gamma_theta")))
    rj_put("decomp_check", rj_num(en("decomp_check")))
    rj_put("decomp_error", rj_num(en("decomp_error")))
    rj_put("F_TZ", rj_num(en("F_TZ")))
    rj_put("F_MZ_T", rj_num(en("F_MZ_T")))
    rj_put("F_TZ_op", rj_num(en("F_TZ_op")))
    rj_put("F_MZ_T_op", rj_num(en("F_MZ_T_op")))
    rj_put("F_TZ_op_critval", rj_num(en("F_TZ_crit")))
    rj_put("F_MZ_T_op_critval", rj_num(en("F_MZ_T_crit")))
    rj_put("op_tau", rj_num(en("op_tau")))
    rj_put("op_alpha", rj_num(en("op_alpha")))
    rj_put("c_TM_hat", rj_num(en("c_TM")))
    rj_put("rho_TM_hat", rj_num(en("c_TM")))
    rj_put("sigma_eT", rj_num(en("sigma_eT")))
    rj_put("sigma_eY", rj_num(en("sigma_eY")))
    rj_put("kappa_star", rj_num(en("kappa_star")))
    rj_put("kappa_star_CI", rj_num(en("kappa_star_CI")))
    rj_put("rho_TY_star", rj_num(en("rho_TY_star")))
    rj_put("f_reduced_form", rj_num(en("f_reduced_form")))
    rj_put("f_first_stage_cond", rj_num(en("f_first_stage_cond")))
    rj_put("RV_theta", rj_num(en("RV_theta")))
    rj_put("XRV_theta", rj_num(en("XRV_theta")))
    rj_put("rv_rf", rj_num(en("rv_rf")))
    rj_put("rv_fs", rj_num(en("rv_fs")))
    rj_put("xrv_rf", rj_num(en("xrv_rf")))
    rj_put("xrv_fs", rj_num(en("xrv_fs")))
    lasso = en("lasso")
    rj_put("lasso", rj_bool(lasso))
    rj_put("lasso_penalty", (lasso ? st_global("e(lasso_penalty)") : "NULL"))
    rj_put("n_clusters", (vce == "cluster" ? rj_num(en("N_clust")) : "NULL"))
    rj_put("K_eff_TZ", rj_num(en("K_eff_TZ")))
    rj_put("K_eff_MZ_T", rj_num(en("K_eff_MZ_T")))
    rj_put("kappa_star_CI_informative", rj_bool(en("kappa_star_CI_informative")))
    rj_put("F_MZ_T_above_crit", rj_bool(en("F_MZ_T_above_crit")))
    M = st_matrix("e(pi1)")
    A = st_matrix("e(se_pi1)")
    if (K == 1) {
        rj_put("pi1", rj_num(M[1]))
        rj_put("se_pi1", rj_num(A[1]))
    }
    else {
        for (i = 1; i <= K; i++) {
            rj_put("pi1." + strofreal(i), rj_num(M[i]))
            rj_put("se_pi1." + strofreal(i), rj_num(A[i]))
        }
    }
    ab = st_global("e(absorb)")
    if (vce == "cluster" & ab != "") {
        M = st_matrix("e(fe_nested)")
        s = ""
        for (i = 1; i <= cols(M); i++) s = s + (i > 1 ? "," : "") + (M[i] ? "TRUE" : "FALSE")
        rj_put("fe_nested_in_cluster", s)
    }
    else rj_put("fe_nested_in_cluster", "NULL")
    rj_put("fe", (ab != "" ? rj_join(tokens(ab), ",") : "NULL"))
    z = en("z_crit")
    est = ("theta", "tau", "gamma", "total", "indirect")
    ses = ("se_theta", "se_tau", "se_gamma", "se_total", "se_indirect")
    cis = ("ci_theta", "ci_tau", "ci_gamma", "ci_total", "ci_indirect")
    for (i = 1; i <= 5; i++) {
        if (en(ses[i]) < .) {
            rj_put(cis[i] + ".1", rj_num(en(est[i]) - z * en(ses[i])))
            rj_put(cis[i] + ".2", rj_num(en(est[i]) + z * en(ses[i])))
        }
        else {
            rj_put(cis[i] + ".1", "NaN")
            rj_put(cis[i] + ".2", "NaN")
        }
    }
    M = st_matrix("e(F_MZ_T_individual)")
    for (i = 1; i <= K; i++) rj_put("F_MZ_T_individual." + strofreal(i), rj_num(M[i]))
    M = st_matrix("e(kappa_adj)")
    for (i = 1; i <= rows(M); i++) {
        s = "kappa_adj." + strofreal(i)
        rj_put(s + ".kappa", rj_num(M[i, 1]))
        rj_put(s + ".kappa_frac", rj_num(M[i, 2]))
        rj_put(s + ".theta_adj", rj_num(M[i, 3]))
        rj_put(s + ".ci_lo", rj_num(M[i, 4]))
        rj_put(s + ".ci_hi", rj_num(M[i, 5]))
    }
    if (en("ar") == 1) {
        rj_ar("theta")
        rj_ar("gamma")
        rj_ar("total")
    }
    else {
        rj_put("ar_theta", "NULL")
        rj_put("ar_gamma", "NULL")
        rj_put("ar_total", "NULL")
    }
    rj_eps("eps_T_hat", "_e_eT")
    rj_eps("eps_M_hat", "_e_eM")
    rj_eps("eps_Y_hat", "_e_eY")
    if (lasso) {
        nm = tokens(st_global("e(lasso_names)"))
        pyn = J(1, cols(nm), "")
        for (i = 1; i <= cols(nm); i++) pyn[i] = rj_pyname(nm[i])
        rj_put("lasso_info.n_selected_union", rj_num(en("lasso_n_selected")))
        rj_put("lasso_info.n_total", rj_num(en("lasso_n_total")))
        M = st_matrix("e(lasso_sel_idx)")
        if (cols(M)) rj_put("lasso_info.union_names", rj_join(pyn[M], "|"))
        else         rj_put("lasso_info.union_names", "")
        eqs = ("Y", "T", "M")
        for (j = 1; j <= K; j++) eqs = (eqs, "Z" + strofreal(j))
        for (j = 1; j <= 3 + K; j++) rj_lasso_detail("lasso_info." + eqs[j], j, pyn)
    }
}

void rj_dumpboot()
{
    real matrix      P, S, C, A, D
    real colvector   v
    string rowvector pn
    real scalar      i, mu, k
    P = st_matrix("e(boot_point)")
    S = st_matrix("e(boot_se)")
    C = st_matrix("e(boot_ci)")
    A = st_matrix("e(boot_analytic_se)")
    D = st_matrix("e(boot_dist)")
    pn = ("total_hat", "tau_hat", "theta_hat", "gamma_hat", "indirect_hat", "kappa_star", "rho_TY_star", "RV_theta", "XRV_theta")
    rj_put("n_failed", rj_num(en("boot_n_failed")))
    rj_put("method", st_global("e(boot_method)"))
    rj_put("B", rj_num(en("boot_B")))
    rj_put("ci_type", st_global("e(boot_citype)"))
    rj_put("alpha", rj_num(en("boot_alpha")))
    rj_put("n_est", rj_num(en("boot_n_est")))
    for (i = 1; i <= rows(P); i++) {
        rj_put("point." + pn[i], rj_num(P[i]))
        rj_put("boot_se." + pn[i], rj_num(S[i]))
        rj_put("boot_ci." + pn[i] + ".1", rj_num(C[i, 1]))
        rj_put("boot_ci." + pn[i] + ".2", rj_num(C[i, 2]))
        if (i <= 4)      rj_put("analytic_se." + pn[i], rj_num(A[i]))
        else if (i == 5) rj_put("analytic_se." + pn[i], rj_numN(A[i]))
        else             rj_put("analytic_se." + pn[i], "NULL")
        v = select(D[., i], D[., i] :< .)
        mu = mean(v)
        rj_put("boot_dist." + pn[i] + ".mean", rj_num(mu))
        rj_put("boot_dist." + pn[i] + ".sd", rj_num(sqrt(sum((v :- mu) :^ 2) / (rows(v) - 1))))
        for (k = 1; k <= 3; k++) rj_put("boot_dist." + pn[i] + ".v" + strofreal(k), rj_num(D[k, i]))
    }
}

// lines of a log file without the command echo and log header
string colvector rj_loglines(string scalar path)
{
    real scalar      fh
    string scalar    line
    string colvector L
    L = J(0, 1, "")
    fh = fopen(path, "r")
    while ((line = fget(fh)) != J(0, 0, "")) {
        if (substr(line, 1, 2) == ". " | substr(line, 1, 2) == "> ") continue
        if (strtrim(line) == "") continue
        if (regexm(line, "^ *(name|log|log type|opened on|closed on): ")) continue
        L = L \ line
    }
    fclose(fh)
    return(L)
}

void rj_stdout(string scalar path)
{
    rj_put("_stdout_lines", strofreal(rows(rj_loglines(path))))
}

// error class from the phrase table (configs.json error_classes, + Stata phrases)
void rj_error(string scalar path, real scalar rc)
{
    string colvector L
    string matrix    T
    string scalar    txt, msg, low
    real scalar      i, j
    L = rj_loglines(path)
    txt = ""
    msg = ""
    for (i = 1; i <= rows(L); i++) {
        low = strlower(strtrim(L[i]))
        if (substr(low, 1, 5) == "note:" | substr(low, 1, 8) == "warning:") continue
        txt = txt + " " + low
        msg = strtrim(L[i])
    }
    T = ("clusters", "at least 2 clusters" \
         "spatial_degenerate", "covers all pairs" \
         "nonfinite", "infinite" \
         "weights", "weights must be non-negative" \
         "weights", "all weights are zero" \
         "too_few_obs", "not enough observations" \
         "too_few_obs", "no observations remaining" \
         "identification", "no variation" \
         "identification", "collinear after" \
         "identification", "collinear with the treatment" \
         "identification", "not identified" \
         "role", "must not contain the outcome" \
         "role", "must differ from the outcome" \
         "role", "three different columns" \
         "role", "duplicate instrument" \
         "option", "used only with se=" \
         "option", "used only with vce(conley)" \
         "option", "must be true or false" \
         "option", "must lie strictly between" \
         "option", "must be one of")
    j = 0
    for (i = 1; i <= rows(T); i++) {
        if (strpos(txt, T[i, 2]) & j == 0) j = i
    }
    rj_put("_status", "error")
    rj_put("_error_class", (j ? T[j, 1] : "other"))
    rj_put("_ERROR", "r(" + strofreal(rc) + "): " + msg)
}

// bootstrap draws: 0-based CSV -> 1-based Stata matrix
void rj_loadidx(string scalar path, string scalar mname)
{
    real scalar      fh
    string scalar    line
    real matrix      M
    real rowvector   r
    M = J(0, 0, .)
    fh = fopen(path, "r")
    while ((line = fget(fh)) != J(0, 0, "")) {
        r = strtoreal(tokens(subinstr(line, ",", " ")))
        if (rows(M) == 0) M = r
        else              M = M \ r
    }
    fclose(fh)
    st_matrix(mname, M :+ 1)
}

// does the Stata estimation sample (marked in _mk) equal the Python rows file?
real scalar rj_checkrows(string scalar path)
{
    real scalar    fh
    string scalar  line
    real colvector py, st
    py = J(0, 1, .)
    fh = fopen(path, "r")
    while ((line = fget(fh)) != J(0, 0, "")) {
        if (strtrim(line) != "") py = py \ strtoreal(line)
    }
    fclose(fh)
    st = st_data(., "_row0", "_mk")
    if (rows(st) != rows(py)) return(0)
    return(all(st :== py))
}

void rj_static(string scalar ver, real matrix T)
{
    string rowvector tl
    real scalar      K, j
    rj_put("version", ver)
    tl = ("0.05", "0.10", "0.20", "0.30")
    for (K = 1; K <= 30; K++) {
        for (j = 1; j <= 4; j++) rj_put("op." + strofreal(K) + "." + tl[j], rj_num(T[K, j]))
    }
}
end

// ---- one config ---------------------------------------------------------------------
// kind: ana (analysis), boot (bootstrap), err (expected to error); for boot,
// bvars lists the columns whose completeness defines the estimation sample
program define runcfg
    args name data cmd kind bvars
    use `data', clear
    if "`kind'" == "boot" {
        mata: rj_loadidx("$REFDIR/bootidx_`name'.csv", "BIDX")
        local cmd `"`cmd' bootidx(BIDX)"'
    }
    else local cmd `"`cmd' residuals(_e)"'
    * printer smoke test: both optional sections on
    local cmd `"`cmd' kappatable r2sensitivity"'
    capture log close _cl
    quietly log using "_cfg.log", text replace name(_cl)
    capture noisily `cmd'
    local rc = _rc
    quietly log close _cl
    mata: rj_begin("`name'")
    if `rc' {
        mata: rj_error("_cfg.log", `rc')
        di as text "[stata] `name': ERROR rc=`rc'"
    }
    else {
        rj_put_ok
        if "`kind'" == "boot" {
            mata: rj_dumpboot()
            mark _mk
            markout _mk `bvars', strok
            foreach v of local bvars {
                capture confirm string variable `v'
                if !_rc qui replace _mk = 0 if strtrim(`v') == ""
            }
            mata: st_local("rowsok", strofreal(rj_checkrows("$REFDIR/bootidx_`name'_rows.csv")))
            mata: rj_put("_rows_check", ("`rowsok'" == "1" ? "ok" : "FAIL"))
            if "`rowsok'" != "1" di as error "[stata] `name': estimation sample differs from the Python rows file"
        }
        else {
            qui gen byte _es = e(sample)
            mata: rj_dumpres()
        }
        mata: rj_stdout("_cfg.log")
        di as text "[stata] `name': ok"
    }
    mata: rj_end()
end

program define rj_put_ok
    mata: rj_put("_status", "ok")
end

program define runna
    args name note
    mata: rj_begin("`name'")
    mata: rj_put("_status", "n/a")
    mata: rj_put("_note", "`note'")
    mata: rj_end()
    di as text "[stata] `name': n/a (`note')"
end

global REFDIR "`REF'"
local REFDIR "`REF'"
mata: rj_open("../out_stata.json")

local M  "mediator(mvar) treatment(tvar)"
local XB "x1 x2 x3"
local IV "ivmediate"

// ============================================================================
//  Analysis configs (configs.json "configs")
// ============================================================================
runcfg base_robust      _main.dta `"`IV' y `XB', `M' instrument(z1) vce(robust)"' ana
local ver "`e(version)'"
runcfg classical        _main.dta `"`IV' y `XB', `M' instrument(z1) vce(unadjusted)"' ana
runcfg cluster          _main.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster clu)"' ana
runcfg spatial_uni_hav  _main.dta `"`IV' y `XB', `M' instrument(z1) vce(conley) latitude(lat) longitude(lon) cutoff(150) kernel(uniform) metric(haversine)"' ana
runcfg spatial_bart_euc _main.dta `"`IV' y `XB', `M' instrument(z1) vce(conley) latitude(lat) longitude(lon) cutoff(1.2) kernel(bartlett) metric(euclidean)"' ana
runcfg wt_aw_robust     _main.dta `"`IV' y `XB' [aw=w], `M' instrument(z1) vce(robust)"' ana
runcfg wt_pw_cluster    _main.dta `"`IV' y `XB' [pw=w], `M' instrument(z1) vce(cluster clu)"' ana
runcfg wt_pw_classical  _main.dta `"`IV' y `XB' [pw=w], `M' instrument(z1) vce(unadjusted)"' ana
runcfg wt_aw_classical  _main.dta `"`IV' y `XB' [aw=w], `M' instrument(z1) vce(unadjusted)"' ana
runcfg fe1_cluster      _main.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster clu) absorb(fe1)"' ana
runcfg fe2way_robust    _main.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) absorb(fe1 fe2)"' ana
runcfg dummies_robust   _main.dta `"`IV' y `XB' i.dcatn i.dnum i.dflo, `M' instrument(z1) vce(robust)"' ana
runcfg glob_string_X    _main.dta `"`IV' y x1 x2 reg_*, `M' instrument(z1) vce(robust)"' ana
runcfg k2_cluster       _main.dta `"`IV' y `XB', `M' instrument(z1 z2) vce(cluster clu) arwidth(50) arpoints(801)"' ana
runcfg weak_ar          _weak.dta `"`IV' y x1, `M' instrument(z1) vce(robust)"' ana
runcfg weak_ar_cluster  _weak.dta `"`IV' y x1, `M' instrument(z1) vce(cluster clu)"' ana
runcfg alpha10_tau20    _main.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) level(90) optau(0.20)"' ana
runcfg op_offgrid       _main.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) opalpha(0.10)"' ana
runcfg noX_robust       _main.dta `"`IV' y, `M' instrument(z1) vce(robust)"' ana
runcfg fe1_weights      _main.dta `"`IV' y `XB' [aw=w], `M' instrument(z1) vce(robust) absorb(fe1)"' ana
runcfg lasso_base       _main.dta `"`IV' y `XB' i.dnum, `M' instrument(z1) vce(robust) lasso showsel"' ana
runcfg lasso_keep_glob  _main.dta `"`IV' y `XB' i.dnum, `M' instrument(z1) vce(robust) lasso lkeep(x1*) expandx(square cube interact_with_X) expanddum(interact_with_X interact_with_dummies)"' ana
runcfg lasso_override   _main.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) lasso llambda(25)"' ana
runcfg lasso_cluster_wt _main.dta `"`IV' y `XB' i.dnum [aw=w], `M' instrument(z1) vce(cluster clu) lasso"' ana
runcfg lasso_spatial    _main.dta `"`IV' y `XB', `M' instrument(z1) vce(conley) latitude(lat) longitude(lon) cutoff(150) kernel(uniform) metric(haversine) lasso"' ana
runna  lasso_expand_override "expand_override (per-variable expansion rules) is not part of the Stata command; expandx()/expanddum() apply to every control"
runcfg lasso_fe         _main.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) absorb(fe1) lasso"' ana
runcfg lasso_k2         _main.dta `"`IV' y `XB', `M' instrument(z1 z2) vce(robust) arwidth(50) arpoints(801) lasso"' ana
runcfg lasso_dcat_probe _main.dta `"`IV' y `XB' i.dcatn, `M' instrument(z1) vce(robust) lasso"' ana
runcfg lasso_dcat_keep  _main.dta `"`IV' y x1 i.dcatn i.dnum, `M' instrument(z1) vce(robust) lasso lkeep(*.dcatn *.dcatn#*) expanddum(interact_with_X interact_with_dummies) showsel"' ana
runcfg glob_class_X     _main.dta `"`IV' y x1 x2 x3, `M' instrument(z1) vce(robust)"' ana
runcfg neg_theta_robust _ext.dta  `"`IV' yneg `XB', `M' instrument(z1) vce(robust)"' ana
runcfg fe1_wt_aw_robust _ext.dta  `"`IV' y `XB' [aw=w], `M' instrument(z1) vce(robust) absorb(fa)"' ana
runcfg fe1_wt_pw_cluster _ext.dta `"`IV' y `XB' [pw=w], `M' instrument(z1) vce(cluster cl) absorb(fa)"' ana
runcfg fe2way_wt_aw_robust _ext.dta `"`IV' y `XB' [aw=w], `M' instrument(z1) vce(robust) absorb(fa fb)"' ana
runcfg fe2way_wt_pw_cluster _ext.dta `"`IV' y `XB' [pw=w], `M' instrument(z1) vce(cluster cl) absorb(fa fb)"' ana
runcfg fe_nested_cluster _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster cl) absorb(fan)"' ana
runcfg fe2way_nested_cluster _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster cl) absorb(fan fb)"' ana
runcfg fe2way_disconnected _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) absorb(fd1 fd2)"' ana
runcfg fe2way_disconnected_cluster _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster cl) absorb(fd1 fd2)"' ana
runcfg k2_het_robust    _ext.dta  `"`IV' y `XB', `M' instrument(z1 z2) vce(robust)"' ana
runcfg k2_het_cluster   _ext.dta  `"`IV' y `XB', `M' instrument(z1 z2) vce(cluster cl)"' ana
runcfg lasso_shift_raw  _ext.dta  `"`IV' y x1-x8, `M' instrument(z1) vce(robust) lasso"' ana
runcfg lasso_shift_ytmz _ext.dta  `"`IV' ys x1-x8, mediator(ms) treatment(ts) instrument(zs) vce(robust) lasso"' ana
runcfg lasso_noexp_raw  _ext.dta  `"`IV' y x1-x8, `M' instrument(z1) vce(robust) lasso expandx(none)"' ana
runcfg lasso_noexp_shifted _ext.dta `"`IV' ys xs1-xs8, mediator(ms) treatment(ts) instrument(zs) vce(robust) lasso expandx(none)"' ana
runcfg lasso_wt_pweight _ext.dta  `"`IV' y x1-x8 [pw=w], `M' instrument(z1) vce(robust) lasso"' ana
runcfg lasso_cluster_only _ext.dta `"`IV' y x1-x8, `M' instrument(z1) vce(cluster cl) lasso"' ana
runcfg lasso_postlasso_loadings _ext.dta `"`IV' y x1-x8 c01-c20, `M' instrument(z1) vce(robust) lasso showsel"' ana
runcfg fe1_zero_weights _ext.dta  `"`IV' y `XB' [aw=wz], `M' instrument(z1) vce(robust) absorb(fz)"' ana
runcfg missing_text_ids _ext.dta  `"`IV' y `XB' i.dtxtn, `M' instrument(z1) vce(cluster cl_txt) absorb(fa_txt)"' ana
runcfg illcond_controls _ext.dta  `"`IV' y x1 year year2, `M' instrument(z1) vce(robust)"' ana
runcfg collinear_control _ext.dta `"`IV' y x1 x2 x12, `M' instrument(z1) vce(robust)"' ana
runcfg fe2way_nested_disconnected _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster cl) absorb(fan fw)"' ana
runcfg spatial_negvar   _negvar.dta `"`IV' y, `M' instrument(z1) vce(conley) latitude(gx) longitude(gy) cutoff(2.5) kernel(uniform) metric(euclidean)"' ana
runcfg fe2_chain        _chain.dta  `"`IV' y, `M' instrument(z1) vce(robust) absorb(wkr firm) noar"' ana

// ============================================================================
//  Bootstrap configs (configs.json "boot_configs"; Python draws replayed)
// ============================================================================
runcfg boot_pairs_basic _main.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) bootstrap bmethod(pairs) bcitype(basic)"' boot "y tvar mvar z1 x1 x2 x3"
runcfg boot_block_pct   _main.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster clu) bootstrap bmethod(block) bcitype(percentile)"' boot "y tvar mvar z1 x1 x2 x3 clu"
runcfg boot_block_stud  _main.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster clu) bootstrap bmethod(block) bcitype(studentized) bdiagnostics"' boot "y tvar mvar z1 x1 x2 x3 clu"
runcfg boot_pairs_missing _ext.dta `"`IV' y x1 x2 xmiss, mediator(mmiss) treatment(tvar) instrument(z1) vce(robust) bootstrap bmethod(pairs) bcitype(percentile)"' boot "y tvar mmiss z1 x1 x2 xmiss"
runcfg boot_block_missing _ext.dta `"`IV' y x1 x2 xmiss, mediator(mmiss) treatment(tvar) instrument(z1) vce(cluster cl_txt) bootstrap bmethod(block) bcitype(studentized)"' boot "y tvar mmiss z1 x1 x2 xmiss cl_txt"
runcfg boot_block_fe_nested _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster cl) absorb(fan) bootstrap bmethod(block) bcitype(studentized)"' boot "y tvar mvar z1 x1 x2 x3 fan cl"

// ============================================================================
//  Degenerate inputs (configs.json "error_configs")
// ============================================================================
runcfg err_one_cluster      _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(cluster cl_one)"' err
runcfg err_z_absorbed_by_fe _ext.dta `"`IV' y `XB', `M' instrument(zg) vce(robust) absorb(fa)"' err
runcfg err_z_collinear      _ext.dta `"`IV' y `XB', `M' instrument(z1 z1dbl) vce(robust)"' err
runcfg err_const_treatment  _ext.dta `"`IV' y `XB', mediator(mvar) treatment(tconst) instrument(z1) vce(robust)"' err
runcfg err_x_contains_y     _ext.dta `"`IV' y x1 y, `M' instrument(z1) vce(robust)"' err
runcfg err_z_equals_t       _ext.dta `"`IV' y `XB', `M' instrument(tvar) vce(robust)"' err
runcfg err_dup_instrument   _ext.dta `"`IV' y `XB', `M' instrument(z1 z1) vce(robust)"' err
runna  err_inf_value "Stata has no infinity: an infinite value cannot be stored in a numeric variable"
runcfg err_negative_weight  _ext.dta `"`IV' y `XB' [aw=wneg], `M' instrument(z1) vce(robust)"' err
runcfg err_n_le_k           _ext.dta `"`IV' y `XB' xsparse, `M' instrument(z1) vce(robust)"' err
runcfg err_spatial_allpairs _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(conley) latitude(lat) longitude(lon) cutoff(500) kernel(uniform) metric(haversine)"' err
runcfg err_spatial_args_robust _ext.dta `"`IV' y `XB', `M' instrument(z1) vce(robust) latitude(lat) longitude(lon)"' err
runcfg err_m_collinear_t    _ext.dta `"`IV' y x1, mediator(mcol) treatment(tvar) instrument(z1) vce(robust)"' err

// ---- static: version, the OP table and computed critical values (from the ado) ---------
ivmediate _optable
matrix OPT = r(table)
mata: rj_begin("_static")
mata: rj_static("`ver'", st_matrix("OPT"))
foreach spec in "1 0.10 0.05 opcv.1.0.1.0.05" "5 0.30 0.05 opcv.5.0.3.0.05" ///
                "30 0.05 0.05 opcv.30.0.05.0.05" "31 0.10 0.05 opcv.31.0.1.0.05" ///
                "1 0.15 0.05 opcv.1.0.15.0.05" "1 0.10 0.01 opcv.1.0.1.0.01" {
    tokenize `spec'
    ivmediate _optable, k(`1') tau(`2') alpha(`3')
    mata: rj_put("`4'", rj_num(st_numscalar("r(cv)")))
}
mata: rj_end()

mata: rj_close()
foreach f in _main _weak _ext _negvar _chain {
    capture erase `f'.dta
}
capture erase _cfg.log
di "HARNESS DONE"
exit, clear
