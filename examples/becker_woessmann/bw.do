* Becker and Woessmann (2009): IV mediation estimates (paper, Table 6, column 1).
* Run from this folder; the paths below are relative to it.

* Uses stata/ivmediate.ado of this repository. To use an installed copy instead, delete the
* adopath line (install: net install ivmediate, from("https://raw.githubusercontent.com/aferrara87/ivmediate/main/stata")).
adopath ++ "../../stata"
import delimited using "becker_woessmann_2009.csv", clear

* The CSV codes missing values as NA, which import delimited reads as text.
destring inctaxpc f_prot f_rw kmwittenberg f_young f_jew f_fem f_ortsgeb f_pruss ///
    f_blind f_deaf f_dumb hhsize lnpop gpop, replace force

* Default options: robust (HC1) SEs and Anderson-Rubin sets, as in the paper.
ivmediate inctaxpc f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop, ///
    treatment(f_prot) mediator(f_rw) instrument(kmwittenberg)
