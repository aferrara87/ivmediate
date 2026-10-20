* Speed comparison, Stata (benchmarks/README.md).
*
* Times ivmediate on data/bench_<N>.csv with the defaults and with the
* post-double-selection LASSO: one untimed warm-up, then 3 timed runs (only the
* estimation command; loading the data and printing are excluded). Writes
* results/stata.csv and results/env_stata.txt.
*
* Run from this folder:  do run_stata.do               (all sizes)
*                        do run_stata.do 1000          (selected sizes)

version 16
clear all
set more off
args sizes_arg
local sizes = cond("`sizes_arg'" == "", "1000 10000 100000", "`sizes_arg'")
capture findfile ivmediate.ado
if _rc {                                   // a checkout without the installed command
    adopath ++ "../stata"
}

tempname res
postfile `res' str8 language long n str8 spec double(t1 t2 t3 median theta gamma tau total) ///
    str12 se_type using "results/_stata_tmp.dta", replace

foreach n of local sizes {
    import delimited using "data/bench_`n'.csv", clear asdouble
    foreach spec in default lasso {
        local opt = cond("`spec'" == "lasso", "lasso", "")
        quietly ivmediate y x1-x50, treatment(t) mediator(m) instrument(z) absorb(g) `opt'
        forvalues i = 1/3 {
            timer clear 1
            timer on 1
            quietly ivmediate y x1-x50, treatment(t) mediator(m) instrument(z) absorb(g) `opt'
            timer off 1
            quietly timer list 1
            local t`i' = r(t1)
        }
        local med = `t1' + `t2' + `t3' - max(`t1', `t2', `t3') - min(`t1', `t2', `t3')
        post `res' ("Stata") (`n') ("`spec'") (`t1') (`t2') (`t3') (`med') ///
            (e(theta)) (e(gamma)) (e(tau)) (e(total)) ("`e(vce)'")
        local ver "`e(version)'"
        display "Stata  n=" %6.0f `n' " " %-7s "`spec'" " median " %8.3f `med' " s"
    }
}
postclose `res'

use "results/_stata_tmp.dta", clear
format t1 t2 t3 median theta gamma tau total %21.17g
local suffix = cond("`sizes_arg'" == "", "", "_" + subinstr("`sizes_arg'", " ", "_", .))
export delimited using "results/stata`suffix'.csv", replace datafmt
erase "results/_stata_tmp.dta"

file open env using "results/env_stata.txt", write replace text
file write env "Stata `c(stata_version)' `c(edition_real)' (processors used: `c(processors)'); ivmediate `ver'" _n
file close env
