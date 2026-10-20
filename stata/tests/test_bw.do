* Smoke test: ivmediate reproduces the Becker and Woessmann (2009) estimates of Ferrara, Gold,
* Heblich, and Pinto (2026), Table 6, column 1, to the printed precision.
* Run from stata/tests/:
*     do test_bw.do              uses stata/ivmediate.ado of this repository
*     do test_bw.do installed    uses the installed copy (net install)
* Stops with an error if a check fails.
version 16
clear
discard
set more off
args mode

if "`mode'" != "installed" {
    adopath ++ ".."
    quietly findfile ivmediate.ado
    if `"`r(fn)'"' != "../ivmediate.ado" {
        display as error `"ivmediate resolves to `r(fn)', not to ../ivmediate.ado"'
        adopath - ".."
        exit 601
    }
}
which ivmediate

import delimited using "../../examples/becker_woessmann/becker_woessmann_2009.csv", clear
destring inctaxpc f_prot f_rw kmwittenberg f_young f_jew f_fem f_ortsgeb f_pruss ///
    f_blind f_deaf f_dumb hhsize lnpop gpop, replace force
capture noisily ivmediate inctaxpc f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb ///
    hhsize lnpop gpop, treatment(f_prot) mediator(f_rw) instrument(kmwittenberg)
local rc = _rc
if "`mode'" != "installed" adopath - ".."
if `rc' exit `rc'

* Table 6, column 1: (stored result, printed value, printed decimals)
local nchk 0
local nfail 0
display as text _n "Check against Table 6, column 1"
foreach spec in "total 0.586 3" "se_total 0.242 3" "tau -0.299 3" "se_tau 0.214 3" ///
                "theta 4.709 3" "se_theta 2.164 3" "gamma 0.188 3" "se_gamma 0.027 3" ///
                "indirect 0.886 3" "se_indirect 0.408 3" "ar_theta_lower 0.43 2" ///
                "ar_theta_upper 10.20 2" "F_TZ 75.69 2" "F_MZ_T 16.36 2" "rho_TY_star -0.257 3" {
    tokenize `spec'
    local ok = (abs(round(e(`1'), 10^(-`3')) - `2') < 1e-9)
    local ++nchk
    if !`ok' local ++nfail
    display as text "  " %-16s "`1'" as result %12.6f e(`1') as text "   Table 6: `2'" ///
        cond(`ok', "", "   <- FAIL")
}
local ok = (e(N) == 426 & "`e(cmd)'" == "ivmediate" & "`e(version)'" == "2.0.0" ///
    & "`e(vce)'" == "robust" & e(ar) == 1 & "`e(ar_theta_shape)'" == "interval")
local ++nchk
if !`ok' local ++nfail
display as text "  N = " e(N) ", e(cmd) = `e(cmd)', e(version) = `e(version)', e(vce) = `e(vce)', " ///
    "AR set for theta: `e(ar_theta_shape)'" cond(`ok', "", "   <- FAIL")

if `nfail' {
    display as error "test_bw: `nfail' of `nchk' checks failed"
    exit 9
}
display as result "test_bw: all `nchk' checks passed"
