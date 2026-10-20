*! ivmediate 2.0.0  17 Oct 2026
*! Andreas Ferrara (maintainer, a.ferrara@pitt.edu), Robert Gold, Stephan Heblich, Rodrigo Pinto
*! Copyright (c) 2026 Andreas Ferrara, Robert Gold, Stephan Heblich, Rodrigo Pinto; MIT license
*! Replaces ivmediate 1.0.4 (SSC, 2020) by Dippel, Ferrara and Heblich; see
*!   help ivmediate, "What's new in 2.0".
*!
*! Causal mediation analysis for linear IV models with a single instrument
*! (Ferrara, Gold, Heblich and Pinto, "Identifying Causal Mediation Effects
*! with a Single Instrument"). Mata port of the Python reference
*! implementation (ivmediate 2.0.0); the cross-language test harness is in
*! validation/ of https://github.com/aferrara87/ivmediate.
*!
*! Implementation: QR least squares on equilibrated designs; collinear
*!   controls dropped; identification and role checks; absorbed-FE df = rank
*!   of the FE dummies, FE nested in the clusters excluded from the CR1
*!   factor; weighted FE demeaning with Irons-Tuck two-way projections;
*!   Olea-Pflueger critical values at the estimated K_eff; closed-form K>=2
*!   AR sets; LASSO with an unpenalized intercept, scale-free weights,
*!   n-based cluster penalty and post-LASSO loadings; bootstrap on the
*!   estimation sample with the user's RNG state restored; e(V) is omitted
*!   when a variance is unavailable (e(V_unavailable) = 1).
*!
*! Conventions: Wald CIs use normal critical values; AR sets use chi2(K);
*!   absorbed FE count as the rank of the FE dummies (G-1 one-way;
*!   G1+G2-c-1 two-way with c connected components); with cluster SEs, FE
*!   dimensions nested in the clusters are not counted in the CR1 factor
*!   (reghdfe convention; with one of two dimensions nested, the other counts
*!   G-1); RV/XRV with K >= 2 use the first instrument.

program define ivmediate, eclass
    version 16.0
    // test hook for the harness: Olea-Pflueger reference table / critical value
    gettoken first rest : 0, parse(" ,")
    if `"`first'"' == "_optable" {
        _ivmediate_optable `rest'
        exit
    }
    _ivmediate_cleanup
    capture noisily _ivmediate_main `0'
    local rc = _rc
    _ivmediate_cleanup
    exit `rc'
end

// remove every scratch object of this command (called on entry and on every exit path)
program define _ivmediate_cleanup
    foreach g in _ivmediate_R _ivmediate_Bdist _ivmediate_Bse _ivmediate_Bok {
        capture mata: rmexternal("`g'")
    }
    foreach s in __ivm_vunavail __ivm_bootB {
        capture scalar drop `s'
    }
end

program define _ivmediate_main, eclass sortpreserve
    version 16.0

    syntax varlist(numeric ts fv min=1) [if] [in] [aweight pweight] , ///
        MEDiator(varname numeric)                                      ///
        TREATment(varname numeric)                                     ///
        INSTrument(varlist numeric min=1)                              ///
        [ Absorb(varlist)                                              ///
          vce(string)                                                  ///
          LATitude(varname numeric)                                    ///
          LONgitude(varname numeric)                                   ///
          CUToff(real -1)                                              ///
          KERNel(string)                                               ///
          METRic(string)                                               ///
          Level(cilevel)                                               ///
          OPTau(real 0.10)                                             ///
          OPAlpha(real 0.05)                                           ///
          noAR                                                         ///
          ARPoints(integer 4001)                                       ///
          ARWidth(real 200)                                            ///
          Bootstrap                                                    ///
          BREps(integer 999)                                           ///
          BSEed(integer 20260317)                                      ///
          BMethod(string)                                              ///
          BCItype(string)                                              ///
          BOOTIDX(name)                                                ///
          BDIAGnostics                                                 ///
          noSTARs                                                      ///
          STARCUToffs(numlist min=3 max=3 >0 <1)                       ///
          LASSO                                                        ///
          LPENALTY(string)                                             ///
          LKEEP(string)                                                ///
          EXPANDX(string)                                              ///
          EXPANDDUM(string)                                            ///
          LCSLACK(real 1.1)                                            ///
          LALPHAP(real -1)                                             ///
          LNITER(integer 15)                                           ///
          LTOL(real 1e-6)                                              ///
          LLAMBDA(real -1)                                             ///
          SHOWSEL                                                      ///
          KAPPAtable                                                   ///
          R2sensitivity                                                ///
          RESIDuals(name) ]

    local cmdline `"ivmediate `0'"'
    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    local doar = ("`ar'" != "noar")
    local lassoflag = ("`lasso'" != "")
    local bootflag  = ("`bootstrap'" != "")

    // ---- roles: compare base names (fv and time-series operators stripped) ----
    _ivmediate_basenames `depvar'
    local bY `s(names)'
    _ivmediate_basenames `treatment'
    local bT `s(names)'
    _ivmediate_basenames `mediator'
    local bM `s(names)'
    if "`bY'" == "`bT'" | "`bY'" == "`bM'" | "`bT'" == "`bM'" {
        di as err "the outcome, treatment and mediator must be three different columns"
        exit 198
    }
    local seen ""
    local dupz ""
    foreach z of local instrument {
        if `: list z in seen' local dupz `dupz' `z'
        local seen `seen' `z'
    }
    if "`dupz'" != "" {
        di as err "duplicate instrument names: `: list uniq dupz'"
        exit 198
    }
    local clash : list instrument & bY
    local clash `clash' `: list instrument & bT'
    local clash `clash' `: list instrument & bM'
    if "`clash'" != "" {
        di as err "instruments must differ from the outcome, treatment and mediator: `clash'"
        exit 198
    }

    // ---- vce(): standard Stata spellings, nothing else ------------------------
    local setype "robust"
    local clustervar ""
    if `"`vce'"' != "" {
        local vce = strtrim(`"`vce'"')
        gettoken vk vrest : vce
        local vk = lower(`"`vk'"')
        local vrest = strtrim(`"`vrest'"')
        local L = strlen("`vk'")
        if (`L' >= 2 & substr("unadjusted", 1, `L') == "`vk'") | "`vk'" == "ols" | "`vk'" == "classical" {
            local setype "classical"
        }
        else if `L' >= 1 & substr("robust", 1, `L') == "`vk'" {
            local setype "robust"
        }
        else if `L' >= 2 & substr("cluster", 1, `L') == "`vk'" {
            local setype "cluster"
            local nclu : word count `vrest'
            if `nclu' != 1 {
                di as err "vce(cluster clustvar) requires exactly one cluster variable"
                exit 198
            }
            unab clustervar : `vrest'
            local vrest ""
        }
        else if "`vk'" == "conley" | "`vk'" == "spatial" {
            local setype "spatial"
        }
        else {
            di as err `"vce() must be one of unadjusted, robust, cluster clustvar or conley; got vce(`vce')"'
            exit 198
        }
        if `"`vrest'"' != "" {
            di as err `"vce(`vk') takes no further arguments: `vrest'"'
            exit 198
        }
    }

    // ---- spatial (Conley) options: an error unless vce(conley) ----------------
    local spgiven ""
    foreach o in latitude longitude kernel metric {
        if "``o''" != "" local spgiven `spgiven' `o'()
    }
    if `cutoff' != -1 local spgiven `spgiven' cutoff()
    local spkern ""
    local spmetric ""
    if "`setype'" == "spatial" {
        if "`latitude'" == "" | "`longitude'" == "" {
            di as err "vce(conley) requires latitude() and longitude()"
            exit 198
        }
        if !(`cutoff' > 0) {
            di as err "vce(conley) requires cutoff() > 0 (km under haversine; coordinate units under euclidean)"
            exit 198
        }
        local spkern = cond("`kernel'" == "", "uniform", lower("`kernel'"))
        local spmetric = cond("`metric'" == "", "haversine", lower("`metric'"))
        if !inlist("`spkern'", "uniform", "bartlett") {
            di as err "kernel() must be one of uniform, bartlett; got `kernel'"
            exit 198
        }
        if !inlist("`spmetric'", "haversine", "euclidean") {
            di as err "metric() must be one of haversine, euclidean; got `metric'"
            exit 198
        }
    }
    else if "`spgiven'" != "" {
        di as err "`spgiven' given, but the spatial options are used only with vce(conley)"
        exit 198
    }

    // ---- controls: patterns never pick up a role column -------------------------
    // Tokens of the varlist as typed that contain * ? ~ or a range (a-b) are
    // patterns; the role columns they match (outcome, treatment, mediator,
    // instruments and the absorb(), cluster, weight and coordinate variables) are
    // excluded with a note. An outcome, treatment, mediator or instrument named
    // explicitly among the controls is an error.
    if "`indepvars'" != "" {
        local roles `bY' `bT' `bM' `instrument'
        local excl `roles' `absorb' `clustervar' `latitude' `longitude'
        if "`weight'" != "" {
            local wexpv = strtrim(subinstr(`"`exp'"', "=", "", 1))
            capture confirm variable `wexpv'
            if !_rc local excl `excl' `wexpv'
        }
        local rest `"`0'"'
        local raw ""
        while 1 {
            gettoken tk rest : rest, parse(" ,[") bind
            if `"`tk'"' == "" | `"`tk'"' == "," | `"`tk'"' == "[" | `"`tk'"' == "if" | `"`tk'"' == "in" continue, break
            local raw `"`raw' `tk'"'
        }
        gettoken rawdep rawx : raw
        local newx ""
        local exclx ""
        foreach tk of local rawx {
            if regexm(`"`tk'"', "[*?~]") | regexm(`"`tk'"', "[A-Za-z0-9_]-[A-Za-z_]") {
                fvunab ex : `tk'
                local keepx ""
                foreach term of local ex {
                    _ivmediate_basenames `term'
                    local bn `s(names)'
                    local hit : list bn & excl
                    if "`hit'" != "" local exclx `exclx' `term'
                    else             local keepx `keepx' `term'
                }
                if "`keepx'" == "" {
                    di as err "controls: the pattern `tk' matched only role columns (`ex'), which cannot be controls"
                    exit 198
                }
                local newx `newx' `keepx'
            }
            else local newx `newx' `tk'
        }
        if "`exclx'" != "" {
            di as text "  Note: control patterns also matched role columns, excluded: `exclx'"
        }
        local indepvars `newx'
        local varlist `depvar' `indepvars'
        _ivmediate_basenames `indepvars'
        local bX `s(names)'
        local bad : list bX & roles
        if "`bad'" != "" {
            di as err "controls must not contain the outcome, treatment, mediator or instruments: `: list uniq bad'"
            exit 198
        }
    }

    // ---- numeric options ----------------------------------------------------------
    if !(`optau' > 0 & `optau' < 1) {
        di as err "optau() must lie strictly between 0 and 1; got `optau'"
        exit 198
    }
    if !(`opalpha' > 0 & `opalpha' < 1) {
        di as err "opalpha() must lie strictly between 0 and 1; got `opalpha'"
        exit 198
    }
    if `arpoints' < 3 {
        di as err "arpoints() must be an integer >= 3"
        exit 198
    }
    if !(`arwidth' > 0) {
        di as err "arwidth() must be positive"
        exit 198
    }

    // ---- enumerated options (case-insensitive; validated whenever given) -------
    if "`bmethod'" != "" {
        local bmethod = lower("`bmethod'")
        if !inlist("`bmethod'", "pairs", "block") {
            di as err "bmethod() must be one of pairs, block; got `bmethod'"
            exit 198
        }
    }
    if "`bcitype'" != "" {
        local bcitype = lower("`bcitype'")
        if !inlist("`bcitype'", "basic", "percentile", "studentized") {
            di as err "bcitype() must be one of basic, percentile, studentized; got `bcitype'"
            exit 198
        }
    }
    if "`lpenalty'" != "" {
        local lpenalty = lower("`lpenalty'")
        if "`lpenalty'" == "cv" {
            di as err "lpenalty(cv) is not supported: cross-validated penalty selection is RNG-dependent and cannot be reproduced across platforms; use lpenalty(plugin) (the default)"
            exit 198
        }
        if "`lpenalty'" != "plugin" {
            di as err "lpenalty() must be one of plugin; got `lpenalty'"
            exit 198
        }
    }
    else local lpenalty "plugin"
    _ivmediate_expops x `"`expandx'"'
    local expx `"`s(ops)'"'
    _ivmediate_expops dum `"`expanddum'"'
    local expd `"`s(ops)'"'

    // ---- bootstrap options --------------------------------------------------------
    if `bootflag' {
        if `lassoflag' {
            di as err "bootstrap with lasso is not supported: post-double-selection uses analytic SEs (bootstrapping the selection under-covers)"
            exit 198
        }
        if "`bmethod'" == "" local bmethod = cond("`setype'" == "cluster", "block", "pairs")
        if "`bmethod'" == "block" & "`setype'" != "cluster" {
            di as err "bmethod(block) requires vce(cluster clustvar)"
            exit 198
        }
        if "`bcitype'" == "" local bcitype "basic"
        if `breps' < 50 {
            di as err "breps() must be at least 50 (999 recommended)"
            exit 198
        }
        if `bseed' < 0 {
            di as err "bseed() must be a non-negative integer"
            exit 198
        }
        if "`bootidx'" != "" {
            capture confirm matrix `bootidx'
            if _rc {
                di as err "bootidx(): matrix `bootidx' not found"
                exit 111
            }
        }
    }
    local bdiag = ("`bdiagnostics'" != "")

    // ---- LASSO tuning ---------------------------------------------------------------
    if `lassoflag' {
        if !(`lcslack' > 0) {
            di as err "lcslack() must be positive"
            exit 198
        }
        if `lalphap' != -1 & !(`lalphap' > 0 & `lalphap' < 1) {
            di as err "lalphap() must lie strictly between 0 and 1"
            exit 198
        }
        if `lniter' < 1 {
            di as err "lniter() must be an integer >= 1"
            exit 198
        }
        if !(`ltol' > 0) {
            di as err "ltol() must be positive"
            exit 198
        }
        if `llambda' != -1 & !(`llambda' >= 0) {
            di as err "llambda() must be >= 0"
            exit 198
        }
    }

    if `: word count `absorb'' > 2 {
        di as err "absorb(): at most two-way fixed effects are supported"
        exit 198
    }

    if "`residuals'" != "" {
        foreach s in eT eM eY {
            confirm new variable `residuals'_`s'
        }
    }

    // ---- weights ------------------------------------------------------------------
    local wtype ""
    tempvar wv obsno
    qui gen long `obsno' = _n
    if "`weight'" != "" {
        local wtype "`weight'"
        qui gen double `wv' `exp'
    }
    if "`wtype'" == "pweight" & "`setype'" == "classical" {
        di as text "  Note: pweight + unadjusted SE is inconsistent. Upgrading to vce(robust) (HC1 sandwich)."
        local setype "robust"
    }
    else if "`wtype'" == "aweight" & "`setype'" == "classical" {
        di as text "  Note: aweight + unadjusted SE assumes the aweight variance model is correct;"
        di as text "        vce(robust) or vce(cluster) is a defensive alternative."
    }

    // ---- estimation sample (before singletons: the sample the bootstrap resamples) ----
    // (mark/markout rather than marksample, which stops on negative weights with
    //  its own message before the weights are checked here)
    tempvar touse
    mark `touse' `if' `in'
    markout `touse' `varlist'
    markout `touse' `mediator' `treatment' `instrument'
    if "`setype'" == "spatial" markout `touse' `latitude' `longitude'
    foreach v in `clustervar' `absorb' {
        markout `touse' `v', strok
        capture confirm string variable `v'
        if !_rc qui replace `touse' = 0 if strtrim(`v') == ""
    }
    if "`wtype'" != "" {
        markout `touse' `wv'
        qui count if `wv' < 0 & `touse'
        if r(N) > 0 {
            di as err "weights must be non-negative"
            exit 402
        }
        qui count if `wv' > 0 & `touse'
        if r(N) == 0 {
            qui count if `touse'
            if r(N) > 0 {
                di as err "all weights are zero -- no observations have positive weight"
                exit 2000
            }
        }
        qui count if `wv' == 0 & `touse'
        local nzero = r(N)
        if `nzero' > 0 {
            di as text "  Note: dropping `nzero' observation(s) with zero weight."
            qui replace `touse' = 0 if `wv' == 0
        }
    }
    qui count if `touse'
    if r(N) == 0 {
        di as err "no observations remaining after dropping missing values"
        exit 2000
    }
    if "`setype'" == "spatial" & "`spmetric'" == "haversine" {
        qui count if (abs(`latitude') > 90 | abs(`longitude') > 180) & `touse'
        if r(N) > 0 {
            di as err "haversine metric requires latitude in [-90,90] and longitude in [-180,180];"
            di as err "`r(N)' observation(s) are out of range. Use metric(euclidean) for projected coordinates."
            exit 198
        }
    }

    // ---- controls -------------------------------------------------------------------
    // Non-LASSO: fvrevar on the full list, aligned with fvexpand; base and
    // omitted factor levels (all-zero columns by construction) are removed.
    // LASSO: continuous terms and main-effect i. factors (one 0/1 indicator per
    // non-base level); the expansion is built in Mata.
    local xvars ""
    local xnames ""
    local ndum 0
    local xc_vars ""
    local xc_names ""
    local xd_vars ""
    local xd_names ""
    if "`indepvars'" != "" & !`lassoflag' {
        fvexpand `indepvars' if `touse'
        local terms `r(varlist)'
        fvrevar `indepvars' if `touse'
        local rvars `r(varlist)'
        local nt : word count `terms'
        forvalues j = 1/`nt' {
            local tm : word `j' of `terms'
            local rv : word `j' of `rvars'
            local isbase 0
            local isfac 0
            local parts = subinstr("`tm'", "#", " ", .)
            local npart : word count `parts'
            foreach pt of local parts {
                local dp = strrpos("`pt'", ".")
                if `dp' > 0 {
                    local pre = substr("`pt'", 1, `dp' - 1)
                    if strpos("`pre'", "o") local isbase 1
                    if strpos("`pre'", "b") & !strpos("`pre'", "bn") local isbase 1
                    if regexm("`pre'", "^[0-9]") local isfac 1
                }
            }
            if !`isbase' {
                local xvars `xvars' `rv'
                local xnames `xnames' `tm'
                if `isfac' & `npart' == 1 local ++ndum
            }
        }
    }
    else if "`indepvars'" != "" {
        foreach tm of local indepvars {
            if strpos("`tm'", "#") > 0 {
                di as err "lasso: control term {bf:`tm'} is an interaction term; supply only continuous variables and i.factor terms (the LASSO builds interactions itself)"
                exit 198
            }
            fvexpand `tm' if `touse'
            local ex `r(varlist)'
            if `: word count `ex'' > 1 | regexm("`ex'", "^[0-9]+[a-z]*\.") {
                foreach term of local ex {
                    if regexm("`term'", "^([0-9]+)([a-z]*)\.(.+)$") {
                        local lev = regexs(1)
                        local mk  = regexs(2)
                        local fv  = regexs(3)
                        if strpos("`mk'", "b") == 0 & strpos("`mk'", "o") == 0 {
                            tempvar di_
                            qui gen byte `di_' = (`fv' == `lev') if `touse'
                            local xd_vars `xd_vars' `di_'
                            local xd_names `xd_names' `lev'.`fv'
                            local ++ndum
                        }
                    }
                }
            }
            else {
                fvrevar `tm' if `touse'
                local xc_vars `xc_vars' `r(varlist)'
                local xc_names `xc_names' `tm'
            }
        }
    }

    // ---- FE and cluster ids (dense codes on the sample) -----------------------------
    local feids ""
    foreach fv of local absorb {
        tempvar gid
        qui egen long `gid' = group(`fv') if `touse'
        local feids `feids' `gid'
    }
    local clid ""
    if "`setype'" == "cluster" {
        tempvar clid
        qui egen long `clid' = group(`clustervar') if `touse'
    }

    // ---- run the engine -----------------------------------------------------------
    sort `obsno'
    tempvar keep
    qui gen byte `keep' = 0
    tempname bmat Vmat
    local yvar   "`depvar'"
    if strpos("`depvar'", ".") {
        tsrevar `depvar'
        local yvar `r(varlist)'
    }
    local resstub "`residuals'"
    mata: _ivmediate_main()

    qui replace `touse' = 0 if !`keep'
    qui count if `touse'
    local N = r(N)

    // ---- bootstrap (optional): replicates replay the analysis on the estimation
    //      sample; a failed replicate is counted, not fatal. The user's RNG state
    //      is saved before set seed and restored afterward, also after errors.
    if `bootflag' {
        local rngstate = c(rngstate)
        capture noisily {
            if "`bootidx'" == "" set seed `bseed'
            mata: _ivmediate_boot_init()
            local B = scalar(__ivm_bootB)
            di as text "  Bootstrap (`bmethod', B = `B')..."
            forvalues b = 1/`B' {
                capture mata: _ivmediate_boot_rep(`b')
            }
            mata: _ivmediate_boot_finish()
        }
        local brc = _rc
        set rngstate `rngstate'
        if `brc' exit `brc'
    }

    // ---- post e(): without e(V) when a variance is unavailable (no r(504)) -------
    mata: _ivmediate_bV("`bmat'", "`Vmat'")
    local vunavail = scalar(__ivm_vunavail)
    if `vunavail' {
        ereturn post `bmat', esample(`touse') obs(`N') depname("`depvar'")
    }
    else {
        ereturn post `bmat' `Vmat', esample(`touse') obs(`N') depname("`depvar'")
    }
    ereturn scalar V_unavailable = `vunavail'
    mata: _ivmediate_post()
    ereturn local title    "Linear IV Mediation Analysis"
    ereturn local cmd      "ivmediate"
    ereturn local cmdline  `"`cmdline'"'
    ereturn local version  "2.0.0"
    ereturn local depvar   "`depvar'"
    ereturn local treat    "`treatment'"
    ereturn local med      "`mediator'"
    ereturn local inst     "`instrument'"
    // e(vce): Stata spelling; e(se_type): the name used by the other implementations
    local vce_st = cond("`setype'" == "classical", "unadjusted", cond("`setype'" == "spatial", "conley", "`setype'"))
    ereturn local vce      "`vce_st'"
    ereturn local vcetype  = cond("`setype'" == "classical", "", cond("`setype'" == "spatial", "Conley", "Robust"))
    ereturn local se_type  "`setype'"
    if "`setype'" == "cluster" ereturn local clustvar "`clustervar'"
    if "`setype'" == "spatial" {
        ereturn scalar spatial_cutoff = `cutoff'
        ereturn local  spatial_kernel "`spkern'"
        ereturn local  spatial_metric "`spmetric'"
        ereturn local  latitude  "`latitude'"
        ereturn local  longitude "`longitude'"
    }
    if "`wtype'" != "" {
        ereturn local wtype "`wtype'"
        ereturn local wexp  `"`exp'"'
    }
    if "`absorb'" != "" ereturn local absorb "`absorb'"
    ereturn local properties = cond(`vunavail', "b", "b V")
    ereturn scalar op_tau   = `optau'
    ereturn scalar op_alpha = `opalpha'
    ereturn scalar level    = `level'
    if `lassoflag' {
        ereturn local lasso_penalty   "`lpenalty'"
        ereturn local lasso_expandx   "`expx'"
        ereturn local lasso_expanddum "`expd'"
    }
    if `bootflag' {
        ereturn local boot_method "`bmethod'"
        ereturn local boot_citype "`bcitype'"
        if "`bootidx'" == "" ereturn scalar boot_seed = `bseed'
        else                 ereturn scalar boot_seed = .
    }
    if "`residuals'" != "" {
        tempvar es
        qui gen byte `es' = e(sample)
        mata: _ivmediate_resid("`residuals'", "`es'")
    }

    // ---- display --------------------------------------------------------------------
    local dostars = ("`stars'" == "")
    if "`starcutoffs'" != "" {
        numlist "`starcutoffs'", sort
        local sc `r(numlist)'
        local cL : word 1 of `sc'
        local cM : word 2 of `sc'
        local cH : word 3 of `sc'
    }
    else {
        local cL 0.01
        local cM 0.05
        local cH 0.10
    }
    mata: _ivmediate_show(`dostars', `cL', `cM', `cH', `=("`showsel'" != "")', `=("`kappatable'" != "")', `=("`r2sensitivity'" != "")')
end

// reference (not used by the estimation): Olea-Pflueger (2013) Table 1 in
// r(table), and with k() the simplified-test critical value at K_eff = k in r(cv)
program define _ivmediate_optable, rclass
    version 16.0
    syntax [, K(integer 0) TAU(real 0.1) Alpha(real 0.05)]
    tempname T
    mata: st_matrix("`T'", _ivm_optable())
    if `k' > 0 {
        mata: st_local("cv", strofreal(_ivm_opcrit(`k', `tau', `alpha'), "%21.17g"))
        return scalar cv = `cv'
    }
    return matrix table = `T'
end

// base variable names of a list of fv/ts terms (each component of an interaction)
program define _ivmediate_basenames, sclass
    local out ""
    foreach tm of local 0 {
        local parts = subinstr("`tm'", "#", " ", .)
        foreach pt of local parts {
            local nm = regexr("`pt'", "^.*\.", "")
            local out `out' `nm'
        }
    }
    sreturn local names `: list uniq out'
end

// validate and canonicalize the LASSO expansion operators (case-insensitive)
program define _ivmediate_expops, sclass
    args kind ops
    if `"`ops'"' == "" {
        if "`kind'" == "x" sreturn local ops "square cube"
        else               sreturn local ops "interact_with_X"
        exit
    }
    local out ""
    foreach op of local ops {
        local lo = lower("`op'")
        if "`lo'" == "none" continue
        if "`kind'" == "x" & inlist("`lo'", "square", "cube") local out `out' `lo'
        else if "`lo'" == "interact_with_x" local out `out' interact_with_X
        else if "`kind'" == "dum" & "`lo'" == "interact_with_dummies" local out `out' interact_with_dummies
        else {
            if "`kind'" == "x" di as err "expandx(): unknown operator {bf:`op'}; must be one of square, cube, interact_with_X, none"
            else               di as err "expanddum(): unknown operator {bf:`op'}; must be one of interact_with_X, interact_with_dummies, none"
            exit 198
        }
    }
    if "`out'" == "" local out "none"
    sreturn local ops "`out'"
end

// =============================================================================
//  Mata engine
// =============================================================================
mata:
mata set matastrict off

// ---- variance regime shared by every regression of one analysis -------------
struct _ivm_se {
    string scalar  type            // classical robust cluster spatial none
    real colvector clo             // cluster: row order (sorted by cluster, stable)
    real matrix    clinfo          // cluster: panelsetup() of the sorted ids
    real scalar    G               // cluster: number of clusters
    real colvector sri, scj, swt   // spatial: kernel triplets (row, col, weight), self pairs included
}

// ---- one least-squares / 2SLS fit -----------------------------------------------
struct _ivm_reg {
    real colvector coef, se, resid
    real matrix    vcov, Q, T      // Xhat*T = Q, (Xhat'Xhat)^-1 = T*T'
    real scalar    n, k
}

// ---- group structure of one id column on the current sample --------------------
struct _ivm_grp {
    real colvector o, gid          // stable order by id; dense group index per row
    real matrix    info            // panelsetup() of id[o]
    real scalar    G
}

// ---- error exit: message in red, return code to the caller ---------------------
void _ivm_err(real scalar rc, string scalar msg)
{
    errprintf("%s\n", msg)
    exit(rc)
}

string scalar _ivm_g(real scalar x)
{
    return(strtrim(sprintf("%21.17g", x)))
}

// ---- QR of the equilibrated design ------------------------------------------------
// A*diag(scale) = Q*R with unit-norm columns (zero columns keep scale 1). A design
// that is numerically singular (|R_jj| <= 1e-12*max|R_ii|) is an error, as in
// the reference implementation (numpy LinAlgError).
void _ivm_eqr(real matrix A, real matrix Q, real matrix R, real rowvector scale)
{
    real rowvector nrm, tau
    real matrix    H, A2
    real colvector d
    real scalar    j

    nrm   = sqrt(colsum(A :* A))
    scale = J(1, cols(A), 1)
    for (j = 1; j <= cols(A); j++) {
        if (nrm[j] > 0 & nrm[j] < .) scale[j] = 1 / nrm[j]
    }
    A2 = A :* scale
    H = .
    tau = .
    R = .
    hqrd(A2, H, tau, R)
    Q = hqrdq1(H, tau)
    d = abs(diagonal(R))
    if (rows(d) > 0) {
        if (missing(d) | min(d) <= 1e-12 * max((max(d), 1e-300))) {
            _ivm_err(459, "singular matrix in a least-squares solve (the design is not of full column rank)")
        }
    }
}

// residuals of the columns of A on B (QR of the equilibrated B)
real matrix _ivm_resid(real matrix A, real matrix B)
{
    real matrix    Q, R
    real rowvector sc
    if (cols(B) == 0) return(A)
    _ivm_eqr(B, Q, R, sc)
    return(A - Q * (Q' * A))
}

// square roots of variances; negative or missing entries give missing
real colvector _ivm_sqrtv(real colvector d)
{
    real colvector s
    real scalar    i
    s = J(rows(d), 1, .)
    for (i = 1; i <= rows(d); i++) {
        if (d[i] >= 0 & d[i] < .) s[i] = sqrt(d[i])
    }
    return(s)
}

// finite and positive definite (a usable variance block)
real scalar _ivm_is_pd(real matrix V)
{
    real colvector ev
    if (rows(V) == 0) return(0)
    if (missing(V)) return(0)
    ev = symeigenvalues((V + V') / 2)
    if (missing(ev)) return(0)
    return(min(ev) > 0)
}

// ---- sandwich variance in the rotated coordinates of the QR solve --------------
real matrix _ivm_vcov(real colvector e, real matrix Q, real matrix T, real scalar n,
                      real scalar k, struct _ivm_se scalar S, real scalar dfa)
{
    real matrix Sc, cs, meat
    real scalar s2, G

    if (S.type == "none") return(J(k, k, .))
    if (S.type == "classical") {
        s2 = (e' * e) / (n - k - dfa)
        return(s2 * (T * T'))
    }
    Sc = Q :* e
    if (S.type == "robust") {
        meat = (Sc' * Sc) * (n / (n - k - dfa))
    }
    else if (S.type == "cluster") {
        G    = S.G
        cs   = panelsum(Sc[S.clo, .], S.clinfo)
        meat = (cs' * cs) * ((G / (G - 1)) * ((n - 1) / (n - k - dfa)))
    }
    else {
        meat = _ivm_spmeat(Sc, Sc, S)
        meat = meat * ((n - 1) / max((n - k - dfa, 1)))
    }
    return(T * meat * T')
}

// spatial meat sum over the kernel pairs of w_ij a_i' b_j, in chunks of 1e6 pairs
real matrix _ivm_spmeat(real matrix A, real matrix B, struct _ivm_se scalar S)
{
    real scalar nnz, a0, a1, cs
    real matrix meat
    nnz  = rows(S.sri)
    cs   = 1000000
    meat = J(cols(A), cols(B), 0)
    for (a0 = 1; a0 <= nnz; a0 = a0 + cs) {
        a1 = min((nnz, a0 + cs - 1))
        meat = meat + cross(A[S.sri[|a0 \ a1|], .], S.swt[|a0 \ a1|], B[S.scj[|a0 \ a1|], .])
    }
    return(meat)
}

// ---- 2SLS (or OLS when is_ols) with the analysis-wide variance regime -----------
struct _ivm_reg scalar _ivm_2sls(real colvector y, real matrix X, real matrix Z,
                                 real scalar is_ols, struct _ivm_se scalar S,
                                 real scalar dfa)
{
    struct _ivm_reg scalar r
    real matrix    Qz, Rz, Xh, Q, R
    real rowvector sz, sc

    r.n = rows(y)
    r.k = cols(X)
    if (is_ols) Xh = X
    else {
        _ivm_eqr(Z, Qz, Rz, sz)
        Xh = Qz * (Qz' * X)
    }
    _ivm_eqr(Xh, Q, R, sc)
    r.coef  = sc' :* solveupper(R, Q' * y)
    r.resid = y - X * r.coef
    r.T     = sc' :* solveupper(R, I(r.k))
    r.Q     = Q
    r.vcov  = _ivm_vcov(r.resid, Q, r.T, r.n, r.k, S, dfa)
    r.se    = _ivm_sqrtv(diagonal(r.vcov))
    return(r)
}

// ---- cross-equation covariance Cov(b1[i1], b2[i2]) (stacked GMM) ----------------
real scalar _ivm_xcov(struct _ivm_reg scalar r1, struct _ivm_reg scalar r2,
                      real scalar i1, real scalar i2, struct _ivm_se scalar S,
                      real scalar dfa)
{
    real colvector psi1, psi2, s1, s2
    real matrix    cb
    real scalar    n, ke, G

    if (S.type == "none") return(.)
    n  = r1.n
    ke = max((r1.k, r2.k))
    psi1 = (r1.Q * r1.T[i1, .]') :* r1.resid
    psi2 = (r2.Q * r2.T[i2, .]') :* r2.resid
    if (S.type == "classical") {
        cb = r1.T * (r1.Q' * r2.Q) * r2.T'
        return(((r1.resid' * r2.resid) / max((n - ke - dfa, 1))) * cb[i1, i2])
    }
    if (S.type == "robust") {
        return(sum(psi1 :* psi2) * (n / (n - ke - dfa)))
    }
    if (S.type == "cluster") {
        G  = S.G
        s1 = panelsum(psi1[S.clo], S.clinfo)
        s2 = panelsum(psi2[S.clo], S.clinfo)
        return(sum(s1 :* s2) * ((G / (G - 1)) * ((n - 1) / (n - ke - dfa))))
    }
    return(_ivm_spmeat(psi1, psi2, S) * ((n - 1) / max((n - ke - dfa, 1))))
}

// ---- spatial kernel triplets (Conley): all pairs within the cutoff ---------------
// Haversine distances in km (degrees in, R = 6371; the great-circle formula of
// a BallTree radius query), euclidean in coordinate units. Built in row blocks.
void _ivm_spatial(real colvector lat, real colvector lon, real scalar cutoff,
                  string scalar kernel, string scalar metric,
                  real colvector ri, real colvector cj, real colvector wt)
{
    real scalar    n, bs, a0, a1, R, nb
    real matrix    D, A, Lb, Lr
    real colvector la, lo, cl, idx, dv, rr
    real rowvector one_n

    n  = rows(lat)
    R  = 6371.0
    bs = max((1, floor(2000000 / n)))
    ri = J(0, 1, .)
    cj = J(0, 1, .)
    wt = J(0, 1, .)
    one_n = J(1, n, 1)
    if (metric == "haversine") {
        la = lat :* (pi() / 180)
        lo = lon :* (pi() / 180)
        cl = cos(la)
    }
    for (a0 = 1; a0 <= n; a0 = a0 + bs) {
        a1 = min((n, a0 + bs - 1))
        nb = a1 - a0 + 1
        if (metric == "euclidean") {
            Lb = lat[|a0 \ a1|] * one_n
            Lr = J(nb, 1, 1) * lat'
            D  = (Lb - Lr) :^ 2
            Lb = lon[|a0 \ a1|] * one_n
            Lr = J(nb, 1, 1) * lon'
            D  = sqrt(D + (Lb - Lr) :^ 2)
        }
        else {
            Lb = la[|a0 \ a1|] * one_n
            Lr = J(nb, 1, 1) * la'
            A  = sin((Lb - Lr) :/ 2) :^ 2
            Lb = lo[|a0 \ a1|] * one_n
            Lr = J(nb, 1, 1) * lo'
            A  = A + (cl[|a0 \ a1|] * cl') :* (sin((Lb - Lr) :/ 2) :^ 2)
            A  = A :* (A :<= 1) + (A :> 1)
            A  = A :* (A :>= 0)
            D  = (2 * R) :* asin(sqrt(A))
        }
        dv  = vec(D')
        idx = selectindex(dv :<= cutoff)
        if (length(idx) == 0) continue
        rr  = ceil(idx :/ n)
        ri  = ri \ (rr :+ (a0 - 1))
        cj  = cj \ (idx - (rr :- 1) :* n)
        if (kernel == "uniform") wt = wt \ J(rows(idx), 1, 1)
        else {
            dv = 1 :- dv[idx] :/ cutoff
            wt = wt \ (dv :* (dv :> 0))
        }
        if (rows(ri) > 1e8) {
            _ivm_err(910, "vce(conley): more than 100 million pairs lie within the cutoff (memory budget); choose a smaller cutoff() or cluster on spatial cells")
        }
    }
}

// ---- group structure ------------------------------------------------------------
struct _ivm_grp scalar _ivm_grp_setup(real colvector id)
{
    struct _ivm_grp scalar g
    real colvector ids, fl, gs
    real scalar    n

    n     = rows(id)
    g.o   = order((id, (1::n)), (1, 2))
    ids   = id[g.o]
    g.info = panelsetup(ids, 1)
    g.G   = rows(g.info)
    fl    = J(n, 1, 0)
    fl[g.info[., 1]] = J(g.G, 1, 1)
    gs    = runningsum(fl)
    g.gid = J(n, 1, .)
    g.gid[g.o] = gs
    return(g)
}

// weighted group means of the columns of V, expanded back to the rows
real matrix _ivm_gmeans(real matrix V, real colvector w, struct _ivm_grp scalar g,
                        real colvector sw)
{
    real matrix m
    m = panelsum((V :* w)[g.o, .], g.info) :/ sw
    return(m[g.gid, .])
}

// population (ddof 0) column standard deviations
real rowvector _ivm_colsd(real matrix V)
{
    real rowvector mu
    mu = colsum(V) / rows(V)
    return(sqrt(colsum((V :- mu) :^ 2) / rows(V)))
}

// root label of every node of a graph given by its edges (union-find)
real colvector _ivm_roots(real colvector a, real colvector b, real scalar N)
{
    real colvector par
    real scalar    e, x, y, i

    par = (1::N)
    for (e = 1; e <= rows(a); e++) {
        x = a[e]
        while (par[x] != x) {
            par[x] = par[par[x]]
            x = par[x]
        }
        y = b[e]
        while (par[y] != y) {
            par[y] = par[par[y]]
            y = par[y]
        }
        if (x != y) par[x] = y
    }
    for (i = 1; i <= N; i++) {
        x = i
        while (par[x] != x) x = par[x]
        par[i] = x
    }
    return(par)
}

// exact two-way FE residualization (weighted): demean within the larger dimension A,
// then solve the Laplacian system of the smaller dimension B (one level fixed per
// connected component) and remove its projection. ok = 0 when not applicable.
real matrix _ivm_demean2_direct(real matrix V, real colvector id1, real colvector id2,
                                real colvector w, real scalar maxgb, real scalar ok)
{
    struct _ivm_grp scalar gA, gB
    real colvector sa, sb, key, o, nab, pa, pb, bs, nv, roots, keepb, kb, seen
    real matrix    infoP, infoA, L, V1, rhs, beta, sol, Bv
    real scalar    n, k, i0, i1, Ga, Gb, r

    ok = 0
    n  = rows(V)
    gA = _ivm_grp_setup(id1)
    gB = _ivm_grp_setup(id2)
    if (gB.G > gA.G) {
        gA = _ivm_grp_setup(id2)
        gB = _ivm_grp_setup(id1)
    }
    Ga = gA.G
    Gb = gB.G
    if (Gb > maxgb) return(V)
    sa = panelsum(w[gA.o], gA.info)
    sb = panelsum(w[gB.o], gB.info)
    // weight of every (A level, B level) pair, sorted by A level
    key = (gA.gid :- 1) :* Gb :+ gB.gid
    o = order((key, (1::n)), (1, 2))
    infoP = panelsetup(key[o], 1)
    nab = panelsum(w[o], infoP)
    pa  = gA.gid[o[infoP[., 1]]]
    pb  = gB.gid[o[infoP[., 1]]]
    L = diag(sb)
    infoA = panelsetup(pa, 1)
    for (k = 1; k <= rows(infoA); k++) {
        i0 = infoA[k, 1]
        i1 = infoA[k, 2]
        bs = pb[|i0 \ i1|]
        nv = nab[|i0 \ i1|]
        L[bs, bs] = L[bs, bs] - (nv * nv') :/ sa[pa[i0]]
    }
    L = (L + L') / 2
    // one B level per connected component is fixed at 0
    roots = _ivm_roots(pa, Ga :+ pb, Ga + Gb)
    keepb = J(Gb, 1, 1)
    seen  = J(Ga + Gb, 1, 0)
    for (k = 1; k <= Gb; k++) {
        r = roots[Ga + k]
        if (!seen[r]) {
            seen[r] = 1
            keepb[k] = 0
        }
    }
    V1  = V - _ivm_gmeans(V, w, gA, sa)
    rhs = panelsum((V1 :* w)[gB.o, .], gB.info)
    beta = J(Gb, cols(V), 0)
    kb = selectindex(keepb)
    if (length(kb)) {
        sol = cholsolve(L[kb, kb], rhs[kb, .])
        if (missing(sol)) sol = lusolve(L[kb, kb], rhs[kb, .])
        if (missing(sol)) return(V)
        beta[kb, .] = sol
    }
    Bv = beta[gB.gid, .]
    ok = 1
    return(V1 - (Bv - _ivm_gmeans(Bv, w, gA, sa)))
}

// one- or two-way FE demeaning with (weighted) group means. Two-way: the exact
// projection (_ivm_demean2_direct) when the smaller dimension has at most 3,000
// levels; otherwise alternating projections with Irons-Tuck acceleration, with
// convergence when the largest change in one sweep, divided by the SD of the
// original column, is below 1e-12 for every column that varies (at most 10,000
// iterations).
real matrix _ivm_demean(real matrix V0, real matrix FE, real colvector w, real scalar quiet)
{
    struct _ivm_grp scalar g1, g2
    real colvector sw1, sw2
    real matrix    V, F, Gs, D1, D2, DD
    real rowvector sd, act, scl, den, cf, mx, aidx
    real scalar    conv, it, chg, tol, maxit, j, ok

    if (cols(FE) == 0 | cols(V0) == 0) return(V0)
    g1  = _ivm_grp_setup(FE[., 1])
    sw1 = panelsum(w[g1.o], g1.info)
    if (cols(FE) == 1) return(V0 - _ivm_gmeans(V0, w, g1, sw1))

    V = _ivm_demean2_direct(V0, FE[., 1], FE[., 2], w, 3000, ok)
    if (ok) return(V)
    g2  = _ivm_grp_setup(FE[., 2])
    sw2 = panelsum(w[g2.o], g2.info)
    tol   = 1e-12
    maxit = 10000
    V   = V0
    sd  = _ivm_colsd(V)
    mx  = colmax(abs(V))
    act = J(1, cols(V), 0)
    scl = J(1, cols(V), 1e-300)
    for (j = 1; j <= cols(V); j++) {
        act[j] = (sd[j] > 1e-14 * max((mx[j], 1e-300)))
        if (sd[j] > 1e-300) scl[j] = sd[j]
    }
    aidx = selectindex(act)
    conv = (length(aidx) == 0)
    chg  = 0
    it   = 0
    while (!conv & it < maxit) {
        F = V - _ivm_gmeans(V, w, g1, sw1)
        F = F - _ivm_gmeans(F, w, g2, sw2)
        it++
        chg = max(colmax(abs(F[., aidx] - V[., aidx])) :/ scl[aidx])
        if (chg < tol) {
            V = F
            conv = 1
            break
        }
        Gs = F - _ivm_gmeans(F, w, g1, sw1)
        Gs = Gs - _ivm_gmeans(Gs, w, g2, sw2)
        it++
        D1  = F - V
        D2  = Gs - F
        DD  = D2 - D1
        den = colsum(DD :* DD)
        cf  = J(1, cols(V), 0)
        for (j = 1; j <= cols(V); j++) {
            if (den[j] > 0) cf[j] = sum(D2[., j] :* DD[., j]) / den[j]
        }
        V = Gs - cf :* D2
    }
    if (!conv & !quiet) {
        printf("{txt}  Warning: two-way FE demeaning did not converge after %g iterations (last relative\n", maxit)
        printf("{txt}  change = %9.2e > tol = 1e-12). Results may be imprecise; check for nearly-collinear FE dimensions.\n", chg)
    }
    return(V)
}

// iterative singleton removal over every FE dimension (reghdfe convention);
// returns the keep mask
real colvector _ivm_singletons(real matrix FE)
{
    real colvector keep, idx, ids, o, sz, sing
    real matrix    info
    real scalar    n, d, dpass, nd

    n = rows(FE)
    keep = J(n, 1, 1)
    if (cols(FE) == 0) return(keep)
    do {
        dpass = 0
        for (d = 1; d <= cols(FE); d++) {
            do {
                idx = selectindex(keep)
                if (length(idx) == 0) return(keep)
                ids  = FE[idx, d]
                o    = order((ids, (1::rows(ids))), (1, 2))
                info = panelsetup(ids[o], 1)
                sz   = info[., 2] - info[., 1] :+ 1
                sing = selectindex(sz :== 1)
                nd   = length(sing)
                if (nd > 0) keep[idx[o[info[sing, 1]]]] = J(nd, 1, 0)
                dpass = dpass + nd
            } while (nd > 0)
        }
    } while (dpass > 0)
    return(keep)
}

// number of connected components of the bipartite FE graph (union-find)
real scalar _ivm_ncomp(real colvector a, real colvector b, real scalar N)
{
    real colvector par
    real scalar    e, x, y, c, i

    par = (1::N)
    for (e = 1; e <= rows(a); e++) {
        x = a[e]
        while (par[x] != x) {
            par[x] = par[par[x]]
            x = par[x]
        }
        y = b[e]
        while (par[y] != y) {
            par[y] = par[par[y]]
            y = par[y]
        }
        if (x != y) par[x] = y
    }
    c = 0
    for (i = 1; i <= N; i++) c = c + (par[i] == i)
    return(c)
}

// absorbed FE df (rank of the FE dummies beyond the intercept), the df that count
// toward CR1 (FE nested in the clusters excluded), total groups, nested flags
void _ivm_fedf(real matrix FE, real colvector cl, real scalar docl, real scalar d,
               real scalar dcr1, real scalar nG, real rowvector nested)
{
    struct _ivm_grp scalar g1, g2
    real matrix    E
    real rowvector G
    real scalar    k, c

    nested = J(1, cols(FE), 0)
    G = J(1, cols(FE), 0)
    for (k = 1; k <= cols(FE); k++) {
        G[k] = rows(uniqrows(FE[., k]))
        if (docl) nested[k] = (rows(uniqrows((FE[., k], cl))) == G[k])
    }
    nG = sum(G)
    if (cols(FE) == 1) {
        d = G[1] - 1
        dcr1 = (nested[1] ? 0 : d)
        return
    }
    g1 = _ivm_grp_setup(FE[., 1])
    g2 = _ivm_grp_setup(FE[., 2])
    E  = uniqrows((g1.gid, g2.gid))
    c  = _ivm_ncomp(E[., 1], E[., 2] :+ g1.G, g1.G + g2.G)
    d  = G[1] + G[2] - c - 1
    if (nested[1] & nested[2])  dcr1 = 0
    else if (nested[1])         dcr1 = G[2] - 1
    else if (nested[2])         dcr1 = G[1] - 1
    else                        dcr1 = d
}

// ---- collinearity: columns kept in order (greedy Gram-Schmidt) -------------------
// Each column, divided by its reference scale, is kept when its residual on
// the previously kept columns (and keep_first) has norm > tol.
real rowvector _ivm_greedy(real matrix A, real matrix K0, real rowvector refs, real scalar tol)
{
    real matrix    B
    real colvector v
    real rowvector kept
    real scalar    j, b, pass, nv, ref, nb

    B  = J(rows(A), 0, .)
    nb = 0
    for (j = 1; j <= cols(K0); j++) {
        v = K0[., j] / max((norm(K0[., j]), 1e-300))
        for (b = 1; b <= nb; b++) v = v - (B[., b]' * v) * B[., b]
        nv = norm(v)
        if (nv > tol & nv < .) {
            B = (B, v / nv)
            nb++
        }
    }
    kept = J(1, cols(A), 0)
    for (j = 1; j <= cols(A); j++) {
        if (cols(refs)) ref = refs[j]
        else            ref = norm(A[., j])
        if (!(ref > 0 & ref < .)) continue
        v = A[., j] / ref
        for (pass = 1; pass <= 2; pass++) {
            for (b = 1; b <= nb; b++) v = v - (B[., b]' * v) * B[., b]
        }
        nv = norm(v)
        if (nv > tol & nv < .) {
            B = (B, v / nv)
            nb++
            kept[j] = 1
        }
    }
    return(kept)
}

// ---- Hajek-form weighted moments ---------------------------------------------------
real scalar _ivm_wcov(real colvector x, real colvector y, real colvector w)
{
    real scalar ws, mx, my
    ws = sum(w)
    mx = sum(w :* x) / ws
    my = sum(w :* y) / ws
    return(sum(w :* (x :- mx) :* (y :- my)) / ws)
}

// weighted least squares on raw data: residuals y - X*b (QR of sqrt(w)-scaled X)
real colvector _ivm_wls_resid(real matrix X, real colvector y, real colvector w)
{
    real colvector sw, b
    real matrix    Q, R
    real rowvector sc
    sw = sqrt(w)
    _ivm_eqr(X :* sw, Q, R, sc)
    b = sc' :* solveupper(R, Q' * (y :* sw))
    return(y - X * b)
}

// ---- Olea-Pflueger (2013) simplified test at the estimated K_eff -----------------
// W2 = Q^(1/2) V Q^(1/2); K_eff = [tr W2]^2 (1 + 2x) / (tr(W2'W2) + 2x tr(W2) lmax(W2)),
// x = 1/tau; critical value = invnchi2(K_eff, x*K_eff, 1 - alpha) / K_eff.
void _ivm_keff(real matrix V, real matrix Qzz, real scalar tau, real scalar alpha,
               real scalar keff, real scalar crit)
{
    real matrix    U, Qh, W2
    real rowvector ev
    real scalar    x, tr, lmax, den

    x = 1 / tau
    keff = .
    crit = .
    if (rows(V) == 1) keff = 1
    else {
        if (!_ivm_is_pd(V) | missing(Qzz)) return
        U = .
        ev = .
        symeigensystem((Qzz + Qzz') / 2, U, ev)
        Qh = (U :* sqrt(ev :* (ev :> 0))) * U'
        W2 = Qh * V * Qh
        W2 = (W2 + W2') / 2
        tr   = trace(W2)
        lmax = max(symeigenvalues(W2))
        den  = trace(W2 * W2) + 2 * x * tr * lmax
        if (!(den > 0 & den < . & tr > 0 & tr < .)) return
        keff = tr * tr * (1 + 2 * x) / den
    }
    crit = invnchi2(keff, x * keff, 1 - alpha) / keff
}

// ---- Olea-Pflueger (2013) Table 1 (alpha = 0.05; tau = 0.05 0.10 0.20 0.30) ----------
// Reference only: the analysis uses the simplified test at the estimated K_eff
// (_ivm_keff), which reproduces this table at K_eff = K to its two decimals.
real matrix _ivm_optable()
{
    real matrix T
    T = (37.42, 23.11, 15.06, 12.05 \ 32.32, 19.29, 12.17,  9.57 \
         30.13, 17.67, 10.95,  8.53 \ 28.85, 16.72, 10.23,  7.92 \
         27.98, 16.08,  9.75,  7.51 \ 27.35, 15.62,  9.40,  7.21 \
         26.86, 15.26,  9.14,  6.98 \ 26.47, 14.97,  8.92,  6.80 \
         26.15, 14.73,  8.74,  6.65 \ 25.87, 14.53,  8.59,  6.52 \
         25.64, 14.36,  8.47,  6.41 \ 25.44, 14.21,  8.36,  6.32 \
         25.26, 14.08,  8.26,  6.24 \ 25.10, 13.96,  8.17,  6.16 \
         24.96, 13.86,  8.10,  6.10 \ 24.83, 13.77,  8.03,  6.04 \
         24.71, 13.68,  7.96,  5.99 \ 24.60, 13.60,  7.91,  5.94 \
         24.50, 13.53,  7.85,  5.89 \ 24.41, 13.46,  7.80,  5.85 \
         24.33, 13.40,  7.76,  5.81 \ 24.25, 13.35,  7.72,  5.78 \
         24.18, 13.29,  7.68,  5.74 \ 24.11, 13.24,  7.64,  5.71 \
         24.05, 13.20,  7.61,  5.68 \ 23.98, 13.15,  7.57,  5.66 \
         23.93, 13.11,  7.54,  5.63 \ 23.87, 13.07,  7.51,  5.61 \
         23.82, 13.04,  7.49,  5.58 \ 23.77, 13.00,  7.46,  5.56)
    return(T)
}

// simplified-test critical value at K_eff = K (the case the table reports)
real scalar _ivm_opcrit(real scalar K, real scalar tau, real scalar alpha)
{
    return(invnchi2(K, K / tau, 1 - alpha) / K)
}

// Wald statistic over K; missing unless the variance block is positive definite
real scalar _ivm_waldK(real colvector b, real matrix V)
{
    if (!_ivm_is_pd(V)) return(.)
    return((b' * lusolve(V, b)) / rows(b))
}

// ---- quantiles (numpy default: linear interpolation, type 7) ----------------------
real scalar _ivm_quantile(real colvector v, real scalar p)
{
    real colvector s
    real scalar    n, h, lo, hi, g, a, b
    s = sort(v, 1)
    n = rows(s)
    if (n == 0) return(.)
    h  = (n * p + (1 + p * (-1))) - 1
    lo = floor(h)
    hi = lo + 1
    if (lo < 0) lo = 0
    if (hi > n - 1) hi = n - 1
    if (lo > n - 1) lo = n - 1
    g = h - floor(h)
    a = s[lo + 1]
    b = s[hi + 1]
    if (g >= 0.5) return(b - (b - a) * (1 - g))
    return(a + (b - a) * g)
}

// sample SD (ddof 1) of the non-missing entries
real scalar _ivm_sd1(real colvector v)
{
    real colvector u
    real scalar    mu
    u = select(v, v :< .)
    if (rows(u) < 2) return(.)
    mu = sum(u) / rows(u)
    return(sqrt(sum((u :- mu) :^ 2) / (rows(u) - 1)))
}

end

// =============================================================================
//  Anderson-Rubin confidence sets
// =============================================================================
mata:

// shape: interval rays rays_left_only rays_right_only line grid empty unavailable
struct _ivm_ar {
    string scalar shape, method
    real scalar   lower, upper, left_upper, right_lower
    real scalar   cv, df, alpha, stat0, pval0
    real scalar   contains0          // 1, 0, or missing (unavailable)
    real matrix   seg                // nseg x 4: lo, hi, open_left, open_right
}

// everything the AR statistic at b0 needs
struct _ivm_arctx {
    real colvector ydep, xend
    real matrix    W, Q, Rq, T
    real rowvector scl, zc
    real scalar    n, k, dfa, cv, closed, tc
    struct _ivm_se scalar S
    real colvector b_c, Bv
    real matrix    V_c, V1, V2
}

// AR regression of (ydep - b0*xend) on W at b0: instrument coefficients and
// their variance block (the QR of W is computed once)
void _ivm_ar_eval(struct _ivm_arctx scalar C, real scalar b0, real colvector bz, real matrix Vz)
{
    real colvector u, cf, e
    real matrix    V
    u  = C.ydep - b0 * C.xend
    cf = C.scl' :* solveupper(C.Rq, C.Q' * u)
    e  = u - C.W * cf
    V  = _ivm_vcov(e, C.Q, C.T, C.n, C.k, C.S, C.dfa)
    bz = cf[C.zc]
    Vz = V[C.zc, C.zc]
}

real scalar _ivm_ar_brute(struct _ivm_arctx scalar C, real scalar b0)
{
    real colvector bz
    real matrix    Vz
    _ivm_ar_eval(C, b0, bz, Vz)
    return(bz' * lusolve(Vz, bz))
}

// closed form: b(b0) = b_c - s*Bv, V(b0) = V_c + s*V1 + s^2*V2, s = b0 - tc
void _ivm_ar_closedbV(struct _ivm_arctx scalar C, real scalar b0, real colvector bb, real matrix VV)
{
    real scalar s
    s  = b0 - C.tc
    bb = C.b_c - s * C.Bv
    VV = C.V_c + s * C.V1 + (s * s) * C.V2
}

real scalar _ivm_ar_closed(struct _ivm_arctx scalar C, real scalar b0)
{
    real colvector bb, sol
    real matrix    VV
    _ivm_ar_closedbV(C, b0, bb, VV)
    sol = lusolve(VV, bb)
    if (missing(sol)) return(_ivm_ar_brute(C, b0))
    return(bb' * sol)
}

// AR(b0) - cv
real scalar _ivm_ar_f(struct _ivm_arctx scalar C, real scalar b0)
{
    if (C.closed) return(_ivm_ar_closed(C, b0) - C.cv)
    return(_ivm_ar_brute(C, b0) - C.cv)
}

// Brent's root finder (port of scipy.optimize.brentq); returns missing when the
// bracket does not change sign
real scalar _ivm_brentq(struct _ivm_arctx scalar C, real scalar xa, real scalar xb,
                        real scalar xtol, real scalar rtol, real scalar maxit)
{
    real scalar xpre, xcur, xblk, fpre, fcur, fblk, spre, scur, sbis, delta
    real scalar stry, dpre, dblk, i

    xpre = xa; xcur = xb
    xblk = 0; fblk = 0; spre = 0; scur = 0
    fpre = _ivm_ar_f(C, xpre)
    fcur = _ivm_ar_f(C, xcur)
    if (fpre >= . | fcur >= .) return(.)
    if (fpre == 0) return(xpre)
    if (fcur == 0) return(xcur)
    if ((fpre < 0) == (fcur < 0)) return(.)
    for (i = 1; i <= maxit; i++) {
        if (fpre != 0 & fcur != 0 & ((fpre < 0) != (fcur < 0))) {
            xblk = xpre
            fblk = fpre
            spre = xcur - xpre
            scur = spre
        }
        if (abs(fblk) < abs(fcur)) {
            xpre = xcur; xcur = xblk; xblk = xpre
            fpre = fcur; fcur = fblk; fblk = fpre
        }
        delta = (xtol + rtol * abs(xcur)) / 2
        sbis  = (xblk - xcur) / 2
        if (fcur == 0 | abs(sbis) < delta) return(xcur)
        if (abs(spre) > delta & abs(fcur) < abs(fpre)) {
            if (xpre == xblk) {
                stry = -fcur * (xcur - xpre) / (fcur - fpre)
            }
            else {
                dpre = (fpre - fcur) / (xpre - xcur)
                dblk = (fblk - fcur) / (xblk - xcur)
                stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
            }
            if (2 * abs(stry) < min((abs(spre), 3 * abs(sbis) - delta))) {
                spre = scur
                scur = stry
            }
            else {
                spre = sbis
                scur = sbis
            }
        }
        else {
            spre = sbis
            scur = sbis
        }
        xpre = xcur
        fpre = fcur
        if (abs(scur) > delta) xcur = xcur + scur
        else                   xcur = xcur + (sbis > 0 ? delta : -delta)
        fcur = _ivm_ar_f(C, xcur)
        if (fcur >= .) return(.)
    }
    return(xcur)
}

// bounded scalar minimization (port of scipy's fminbound); returns x and f(x)
void _ivm_fminbound(struct _ivm_arctx scalar C, real scalar x1, real scalar x2,
                    real scalar xatol, real scalar maxfun, real scalar xmin, real scalar fmin)
{
    real scalar sqrt_eps, gm, a, b, fulc, nfc, xf, rat, e, x, fx, num, fu
    real scalar ffulc, fnfc, xm, tol1, tol2, golden, r, q, p, si

    sqrt_eps = sqrt(2.2e-16)
    gm = 0.5 * (3 - sqrt(5))
    a = x1; b = x2
    fulc = a + gm * (b - a)
    nfc = fulc; xf = fulc
    rat = 0; e = 0
    x  = xf
    fx = _ivm_ar_f(C, x)
    num = 1
    fu = .
    ffulc = fx; fnfc = fx
    xm = 0.5 * (a + b)
    tol1 = sqrt_eps * abs(xf) + xatol / 3
    tol2 = 2 * tol1
    while (abs(xf - xm) > (tol2 - 0.5 * (b - a))) {
        golden = 1
        if (abs(e) > tol1) {
            golden = 0
            r = (xf - nfc) * (fx - ffulc)
            q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r
            q = 2 * (q - r)
            if (q > 0) p = -p
            q = abs(q)
            r = e
            e = rat
            if ((abs(p) < abs(0.5 * q * r)) & (p > q * (a - xf)) & (p < q * (b - xf))) {
                rat = p / q
                x = xf + rat
                if (((x - a) < tol2) | ((b - x) < tol2)) {
                    si = sign(xm - xf) + ((xm - xf) == 0)
                    rat = tol1 * si
                }
            }
            else golden = 1
        }
        if (golden) {
            if (xf >= xm) e = a - xf
            else          e = b - xf
            rat = gm * e
        }
        si = sign(rat) + (rat == 0)
        x  = xf + si * max((abs(rat), tol1))
        fu = _ivm_ar_f(C, x)
        num++
        if (fu >= .) break
        if (fu <= fx) {
            if (x >= xf) a = xf
            else         b = xf
            fulc = nfc; ffulc = fnfc
            nfc = xf;   fnfc = fx
            xf = x;     fx = fu
        }
        else {
            if (x < xf) a = x
            else        b = x
            if ((fu <= fnfc) | (nfc == xf)) {
                fulc = nfc; ffulc = fnfc
                nfc = x;    fnfc = fu
            }
            else if ((fu <= ffulc) | (fulc == xf) | (fulc == nfc)) {
                fulc = x; ffulc = fu
            }
        }
        xm = 0.5 * (a + b)
        tol1 = sqrt_eps * abs(xf) + xatol / 3
        tol2 = 2 * tol1
        if (num >= maxfun) break
    }
    xmin = xf
    fmin = fx
}

// AR confidence set for the coefficient on xend; zc are the instrument columns
// of W (the AR regression's design, instruments included)
struct _ivm_ar scalar _ivm_arset(real colvector ydep, real colvector xend, real matrix W,
        real rowvector zc, struct _ivm_se scalar S, real scalar dfa, real scalar alpha,
        real scalar center, real scalar scale, real scalar hw, real scalar gp)
{
    struct _ivm_ar    scalar out
    struct _ivm_arctx scalar C
    real colvector bz0, b_m, b_c, b_p, b_k, grid, st, lam_lo, lam_sc, ev, bb, sol
    real matrix    Vz0, V_m, V_c, V_p, V_k, VV, segs
    real rowvector cand
    real scalar    K, tc, sc, s_hw, check_neg, use_grid, d, bm, b0c, bp, Vm, V0, Vp
    real scalar    B, bv, cvv, qa, qb, qc, q_pred, q_act, q_scale, v_neg, tiny, disc
    real scalar    sq, q_, r1, r2, s_lo, s_hi, i, j, ng, chk_b, chk_act, chk_pred, exact
    real scalar    unavailable, lo, hi, rt, a_, b_, xatol, xm, fm, lo_, hi_, step, test_in
    real rowvector ss
    real colvector stv, ord

    K  = cols(zc)
    out.cv    = invchi2(K, 1 - alpha)
    out.df    = K
    out.alpha = alpha
    out.lower = .; out.upper = .; out.left_upper = .; out.right_lower = .
    out.seg   = J(0, 4, .)
    out.contains0 = .
    out.method = ""
    out.shape  = ""

    C.ydep = ydep
    C.xend = xend
    C.W    = W
    C.n    = rows(W)
    C.k    = cols(W)
    C.zc   = zc
    C.S    = S
    C.dfa  = dfa
    C.cv   = out.cv
    C.closed = 0
    _ivm_eqr(W, C.Q, C.Rq, C.scl)
    C.T    = C.scl' :* solveupper(C.Rq, I(C.k))

    // AR test of H0: beta = 0
    _ivm_ar_eval(C, 0, bz0, Vz0)
    if (_ivm_is_pd(Vz0)) {
        out.stat0 = bz0' * lusolve(Vz0, bz0)
        out.pval0 = chi2tail(K, out.stat0)
    }
    else {
        out.stat0 = .
        out.pval0 = .
    }

    tc = (center < . ? center : 0)
    sc = ((scale < . & scale > 0) ? scale : max((abs(tc), 1)))
    s_hw = hw * sc
    check_neg = (S.type == "spatial")
    C.tc = tc
    use_grid = (K > 1)

    if (!use_grid) {
        // ---- exact Fieller inversion (K = 1) ----
        d = sc
        _ivm_ar_eval(C, tc - d, b_m, V_m)
        _ivm_ar_eval(C, tc,     b_c, V_c)
        _ivm_ar_eval(C, tc + d, b_p, V_p)
        bm = b_m[1]; b0c = b_c[1]; bp = b_p[1]
        Vm = V_m[1, 1]; V0 = V_c[1, 1]; Vp = V_p[1, 1]
        B   = (bm - bp) / (2 * d)
        bv  = (Vp - Vm) / (2 * d)
        cvv = (Vp - 2 * V0 + Vm) / (2 * d * d)
        qa = B * B - out.cv * cvv
        qb = -2 * b0c * B - out.cv * bv
        qc = b0c * b0c - out.cv * V0
        _ivm_ar_eval(C, tc + 2 * d, b_k, V_k)
        q_pred = qa * (2 * d)^2 + qb * (2 * d) + qc
        q_act  = b_k[1]^2 - out.cv * V_k[1, 1]
        q_scale = max((abs(qa) * (2 * d)^2, abs(qb) * (2 * d), abs(qc), abs(q_act), 1e-300))
        if (q_pred >= . | abs(q_pred - q_act) > 1e-6 * q_scale) {
            use_grid = 1
        }
        else {
            out.method = "fieller"
            v_neg = 0
            if (check_neg) {
                ss = (-s_hw, 0, s_hw)
                if (cvv > 0 & abs(bv / (2 * cvv)) <= s_hw) ss = (ss, -bv / (2 * cvv))
                for (i = 1; i <= cols(ss); i++) {
                    if (V0 + bv * ss[i] + cvv * ss[i] * ss[i] <
                        -1e-12 * (abs(V0) + abs(bv * ss[i]) + abs(cvv) * ss[i] * ss[i])) v_neg = 1
                }
            }
            tiny = 1e-14 * max((abs(qb) * d, abs(qc), 1e-300))
            if (v_neg) out.shape = "unavailable"
            else if (abs(qa) * d * d <= tiny) {
                if (abs(qb) * d <= tiny) out.shape = (qc <= 0 ? "line" : "empty")
                else if (qb > 0) {
                    out.shape = "rays_left_only"
                    out.left_upper = tc + (-qc / qb)
                }
                else {
                    out.shape = "rays_right_only"
                    out.right_lower = tc + (-qc / qb)
                }
            }
            else {
                disc = qb * qb - 4 * qa * qc
                if (disc < 0) out.shape = (qa > 0 ? "empty" : "line")
                else {
                    sq = sqrt(disc)
                    if (qb != 0) q_ = -0.5 * (qb + (qb >= 0 ? sq : -sq))
                    else         q_ = 0.5 * sq
                    r1 = q_ / qa
                    r2 = (q_ != 0 ? qc / q_ : -qb / qa - r1)
                    s_lo = min((r1, r2))
                    s_hi = max((r1, r2))
                    if (qa > 0) {
                        out.shape = "interval"
                        out.lower = tc + s_lo
                        out.upper = tc + s_hi
                    }
                    else {
                        out.shape = "rays"
                        out.left_upper  = tc + s_lo
                        out.right_lower = tc + s_hi
                    }
                }
            }
        }
    }

    if (use_grid) {
        // ---- closed-form statistic on a grid, boundaries refined (K > 1) ----
        out.method = "grid"
        ng   = gp
        step = (2 * hw) / (ng - 1)
        grid = (-hw :+ (0::ng - 1) :* step)
        grid[ng] = hw
        grid = tc :+ grid :* sc
        d = sc
        _ivm_ar_eval(C, tc - d, b_m, V_m)
        _ivm_ar_eval(C, tc,     b_c, V_c)
        _ivm_ar_eval(C, tc + d, b_p, V_p)
        C.b_c = b_c
        C.V_c = V_c
        C.Bv  = (b_m - b_p) / (2 * d)
        C.V1  = (V_p - V_m) / (2 * d)
        C.V2  = (V_p - 2 * V_c + V_m) / (2 * d * d)
        chk_b   = tc + 2 * d
        chk_act = _ivm_ar_brute(C, chk_b)
        chk_pred = _ivm_ar_closed(C, chk_b)
        exact = (chk_act < . & chk_pred < . & abs(chk_pred - chk_act) <= 1e-6 * max((abs(chk_act), 1)))
        C.closed = exact
        st = J(ng, 1, .)
        lam_lo = J(ng, 1, .)
        lam_sc = J(ng, 1, .)
        for (i = 1; i <= ng; i++) {
            if (exact) {
                _ivm_ar_closedbV(C, grid[i], bb, VV)
                sol = lusolve(VV, bb)
                st[i] = (missing(sol) ? _ivm_ar_brute(C, grid[i]) : bb' * sol)
            }
            else {
                _ivm_ar_eval(C, grid[i], bb, VV)
                st[i] = (_ivm_is_pd(VV) ? bb' * lusolve(VV, bb) : .)
            }
            if (check_neg) {
                ev = symeigenvalues((VV + VV') / 2)'
                lam_lo[i] = min(ev)
                lam_sc[i] = max(abs(ev))
                if (missing(ev)) lam_lo[i] = .
            }
        }
        unavailable = 0
        if (check_neg) {
            for (i = 1; i <= ng; i++) {
                if (lam_lo[i] >= . | lam_lo[i] < -1e-12 * lam_sc[i]) unavailable = 1
            }
        }
        if (unavailable) out.shape = "unavailable"
        else {
            segs = J(0, 4, .)
            i = 1
            while (i <= ng) {
                if (st[i] <= out.cv) {
                    j = i
                    while (j + 1 <= ng) {
                        if (st[j + 1] <= out.cv) j++
                        else break
                    }
                    lo = grid[i]
                    hi = grid[j]
                    if (i > 1) {
                        rt = _ivm_brentq(C, grid[i - 1], grid[i],
                                         1e-14 * max((abs(grid[i - 1]), abs(grid[i]), sc)), 1e-10, 200)
                        if (rt < .) lo = rt
                    }
                    if (j < ng) {
                        rt = _ivm_brentq(C, grid[j], grid[j + 1],
                                         1e-14 * max((abs(grid[j]), abs(grid[j + 1]), sc)), 1e-10, 200)
                        if (rt < .) hi = rt
                    }
                    segs = segs \ (lo, hi, (i == 1), (j == ng))
                    i = j + 1
                }
                else i++
            }
            // accepted intervals narrower than the grid step: minimize the exact
            // statistic around every interior local minimum (lowest first, <= 50)
            if (ng >= 3) {
                cand = J(1, 0, .)
                for (i = 2; i <= ng - 1; i++) {
                    if (st[i] < . & st[i - 1] < . & st[i + 1] < . & st[i] <= st[i - 1] & st[i] <= st[i + 1] & !(st[i] <= out.cv)) {
                        cand = (cand, i)
                    }
                }
                if (cols(cand) > 0) {
                    stv = st[cand']
                    ord = order((stv, (1::rows(stv))), (1, 2))
                    cand = cand[ord']
                    if (cols(cand) > 50) cand = cand[|1 \ 50|]
                    for (j = 1; j <= cols(cand); j++) {
                        i  = cand[j]
                        a_ = grid[i - 1]
                        b_ = grid[i + 1]
                        xatol = 1e-12 * max((abs(a_), abs(b_), sc))
                        _ivm_fminbound(C, a_, b_, xatol, 500, xm, fm)
                        if (fm < 0) {
                            lo_ = _ivm_brentq(C, a_, xm, 1e-14 * max((abs(a_), abs(xm), sc)), 1e-10, 200)
                            hi_ = _ivm_brentq(C, xm, b_, 1e-14 * max((abs(xm), abs(b_), sc)), 1e-10, 200)
                            if (lo_ < . & hi_ < .) segs = segs \ (lo_, hi_, 0, 0)
                        }
                    }
                }
                if (rows(segs) > 1) {
                    ord  = order((segs[., 1], (1::rows(segs))), (1, 2))
                    segs = segs[ord, .]
                }
            }
            out.seg   = segs
            out.shape = (rows(segs) > 0 ? "grid" : "empty")
        }
    }

    // ---- does the set contain zero? ----
    if (out.shape == "unavailable") {
        out.contains0 = .
        return(out)
    }
    if (out.shape == "interval")             out.contains0 = (out.lower <= 0 & 0 <= out.upper)
    else if (out.shape == "rays")            out.contains0 = (0 <= out.left_upper | 0 >= out.right_lower)
    else if (out.shape == "rays_left_only")  out.contains0 = (0 <= out.left_upper)
    else if (out.shape == "rays_right_only") out.contains0 = (0 >= out.right_lower)
    else if (out.shape == "line")            out.contains0 = 1
    else if (out.shape == "grid") {
        out.contains0 = 0
        for (i = 1; i <= rows(out.seg); i++) {
            if (out.seg[i, 1] <= 0 & 0 <= out.seg[i, 2]) out.contains0 = 1
        }
    }
    else out.contains0 = 0
    if (out.stat0 < .) {
        test_in = (out.stat0 <= out.cv)
        if (test_in != out.contains0) out.contains0 = test_in
    }
    return(out)
}

// unbounded (rays, line, a grid segment reaching the grid boundary) or disconnected
real scalar _ivm_ar_unbdisc(struct _ivm_ar scalar ar)
{
    if (anyof(("rays", "rays_left_only", "rays_right_only", "line"), ar.shape)) return(1)
    if (ar.shape == "grid") {
        if (rows(ar.seg) > 1) return(1)
        if (rows(ar.seg) == 1) return(ar.seg[1, 3] | ar.seg[1, 4])
    }
    return(0)
}

end

// =============================================================================
//  Post-double-selection LASSO (BCH plug-in penalty; paper App. H)
// =============================================================================
mata:

struct _ivm_lres {
    real colvector beta, sel          // p x 1 (unscaled coefficients; selection 0/1)
    real scalar    lambda, n_iter, converged
    string scalar  method
}

// coordinate descent for (1/(2n))||y - Xb||^2 + alpha*||b||_1: a line-by-line
// port of scikit-learn's enet_coordinate_descent (cyclic, l2 = 0) with its
// duality-gap stopping rule
real colvector _ivm_cd(real matrix X, real colvector y, real scalar alpha_sk,
                       real scalar maxit, real scalar tol0)
{
    real scalar    n, p, it, ii, w_ii, tmp, w_max, d_w_max, d_w_ii, alpha, tol
    real scalar    dual_norm, R_norm2, cnst, A_norm2, gap
    real colvector w, normc, R, XtA

    n = rows(X)
    p = cols(X)
    alpha = alpha_sk * n
    w = J(p, 1, 0)
    normc = colsum(X :^ 2)'
    R = y
    tol = tol0 * (y' * y)
    for (it = 1; it <= maxit; it++) {
        w_max = 0
        d_w_max = 0
        for (ii = 1; ii <= p; ii++) {
            if (normc[ii] == 0) continue
            w_ii = w[ii]
            if (w_ii != 0) R = R + w_ii * X[., ii]
            tmp = X[., ii]' * R
            w[ii] = sign(tmp) * max((abs(tmp) - alpha, 0)) / normc[ii]
            if (w[ii] != 0) R = R - w[ii] * X[., ii]
            d_w_ii = abs(w[ii] - w_ii)
            if (d_w_ii > d_w_max) d_w_max = d_w_ii
            if (abs(w[ii]) > w_max) w_max = abs(w[ii])
        }
        if (w_max == 0 | d_w_max / w_max < tol0 | it == maxit) {
            XtA = X' * R
            dual_norm = max(abs(XtA))
            R_norm2 = R' * R
            if (dual_norm > alpha) {
                cnst = alpha / dual_norm
                A_norm2 = R_norm2 * (cnst^2)
                gap = 0.5 * (R_norm2 + A_norm2)
            }
            else {
                cnst = 1
                gap = R_norm2
            }
            gap = gap + (alpha * sum(abs(w)) - cnst * (R' * y))
            if (gap < tol) break
        }
    }
    return(w)
}

// per-feature loadings from the current residuals (score w_i W_ij u_i):
// heteroskedastic, cluster or spatial; floored at 1e-10*max(1, max|W|)
real rowvector _ivm_psi(real matrix W, real colvector u, real colvector w,
                        struct _ivm_se scalar S)
{
    real colvector wu
    real matrix    wuW, cs
    real rowvector psi
    real scalar    nw, j, fl

    nw = sum(w)
    wu = w :* u
    if (S.type == "spatial") {
        wuW = W :* wu
        psi = J(1, cols(W), 0)
        for (j = 1; j <= cols(W); j++) {
            psi[j] = _ivm_spmeat(wuW[., j], wuW[., j], S)
        }
        psi = sqrt((psi :* (psi :> 0)) / nw)
    }
    else if (S.type == "cluster") {
        wuW = W :* wu
        cs  = panelsum(wuW[S.clo, .], S.clinfo)
        psi = sqrt(colsum(cs :^ 2) / nw)
    }
    else {
        psi = sqrt(colsum((W :^ 2) :* (wu :^ 2)) / nw)
    }
    fl = 1e-10 * max((1, max(abs(W))))
    for (j = 1; j <= cols(psi); j++) {
        if (psi[j] < fl) psi[j] = fl
    }
    return(psi)
}

// residuals of a weighted least-squares fit of y on the columns of X (SVD)
real colvector _ivm_lsq_resid(real colvector y, real matrix X, real colvector sw)
{
    real colvector b
    if (cols(X) == 0) return(y)
    b = svsolve(X :* sw, y :* sw)
    return(y - X * b)
}

// BCH plug-in LASSO of y on W with weights, an unpenalized intercept (weighted
// centering), keep columns partialed out, post-LASSO residual loadings
struct _ivm_lres scalar _ivm_bch(real colvector y_raw, real matrix W_raw, real colvector w0,
        struct _ivm_se scalar S, real rowvector keep, real scalar niter_max,
        real scalar tol, real scalar c_slack, real scalar alpha_p0, real scalar lam_over)
{
    struct _ivm_lres scalar r
    real colvector w, sw, y, yl, ylp, u, unew, bl, bnew, bk
    real matrix    Wc, Wl, Wlp, Kp, Wk, Bk
    real rowvector nonkeep, psi, csc, kidx, sel
    real scalar    n, p, pl, alpha_p, neff, lam, it, mc, j, w1

    n = rows(W_raw)
    p = cols(W_raw)
    alpha_p = (alpha_p0 < . ? alpha_p0 : 0.1 / max((ln(max((n, 2))), 1)))
    r.beta = J(p, 1, 0)
    r.sel  = J(p, 1, 0)
    r.n_iter = 0
    r.converged = 0
    if (p == 0) {
        r.lambda = 0
        r.converged = 1
        r.method = "trivial"
        return(r)
    }
    w  = w0 / mean(w0)
    y  = y_raw :- sum(w :* y_raw) / sum(w)
    Wc = W_raw :- colsum(W_raw :* w) / sum(w)
    sw = sqrt(w)
    kidx = J(1, 0, .)
    if (cols(keep)) {
        if (any(keep)) kidx = selectindex(keep)
    }
    if (length(kidx) > 0) {
        nonkeep = selectindex(!keep)
        Kp = Wc[., kidx] :* sw
        bk = svsolve(Kp, y :* sw)
        yl = y - Wc[., kidx] * bk
        if (length(nonkeep) > 0) {
            Bk = svsolve(Kp, Wc[., nonkeep] :* sw)
            Wl = Wc[., nonkeep] - Wc[., kidx] * Bk
        }
        else Wl = J(n, 0, .)
    }
    else {
        nonkeep = (1..p)
        yl = y
        Wl = Wc
    }
    pl  = cols(Wl)
    ylp = yl :* sw
    Wlp = Wl :* sw
    // n_eff: n, or the Kish size with unequal weights (also under clustering)
    w1 = w[1]
    if (all(abs(w :- w1) :<= 1e-8 + 1e-5 * abs(w1))) neff = n
    else neff = sum(w)^2 / sum(w :^ 2)

    if (pl == 0) {
        if (length(kidx) > 0) r.sel[kidx'] = J(length(kidx), 1, 1)
        r.lambda = 0
        r.converged = 1
        r.method = "keep-only"
        return(r)
    }
    if (lam_over < .) {
        psi = _ivm_psi(Wl, yl, w, S)
        csc = 1 :/ psi
        lam = lam_over
        bl  = _ivm_cd(Wlp :* csc, ylp, lam / (2 * n), 20000, 1e-6) :* csc'
        r.beta[nonkeep'] = bl
        r.sel[nonkeep']  = (bl :!= 0)
        if (length(kidx) > 0) r.sel[kidx'] = J(length(kidx), 1, 1)
        r.lambda = lam
        r.n_iter = 1
        r.converged = 1
        r.method = "override"
        return(r)
    }
    lam = 2 * c_slack * sqrt(neff) * invnormal(1 - alpha_p / (2 * pl))
    u  = yl
    bl = J(pl, 1, 0)
    for (it = 1; it <= niter_max; it++) {
        r.n_iter = it
        psi  = _ivm_psi(Wl, u, w, S)
        csc  = 1 :/ psi
        bnew = _ivm_cd(Wlp :* csc, ylp, lam / (2 * n), 20000, 1e-6) :* csc'
        sel  = selectindex(bnew' :!= 0)
        if (length(sel)) unew = _ivm_lsq_resid(yl, Wl[., sel], sw)
        else             unew = yl
        mc   = max(abs(bnew - bl))
        bl   = bnew
        u    = unew
        if (mc < tol) {
            r.converged = 1
            break
        }
    }
    r.beta[nonkeep'] = bl
    r.sel[nonkeep']  = (bl :!= 0)
    if (length(kidx) > 0) r.sel[kidx'] = J(length(kidx), 1, 1)
    r.lambda = lam
    r.method = "plugin"
    return(r)
}

// rich expansion of the controls (paper App. H): originals (x, then dummies),
// squares and cubes of continuous x, x*x pairs, dummy*x, dummy*dummy (not
// within one factor); near-constant columns dropped
real matrix _ivm_rich(real matrix Xc, real matrix Xd, real rowvector dsrc,
        string rowvector xcn, string rowvector xdn, real scalar sq, real scalar cb,
        real scalar ixx, real scalar dix, real scalar ddd, string rowvector names)
{
    real scalar    n, kx, kd, j, u, i, nc
    real matrix    W
    real rowvector iscont, mu, cv, am, keep
    string rowvector nm

    n  = max((rows(Xc), rows(Xd)))
    kx = cols(Xc)
    kd = cols(Xd)
    iscont = J(1, kx, 0)
    for (j = 1; j <= kx; j++) iscont[j] = (rows(uniqrows(Xc[., j])) > 2)
    nc = kx + kd
    nc = nc + (sq + cb) * sum(iscont)
    if (ixx) nc = nc + kx * (kx - 1) / 2
    if (dix) nc = nc + kd * kx
    if (ddd) nc = nc + kd * (kd - 1) / 2
    W  = J(n, nc, .)
    nm = J(1, nc, "")
    nc = 0
    for (j = 1; j <= kx; j++) {
        nc++; W[., nc] = Xc[., j]; nm[nc] = xcn[j]
    }
    for (j = 1; j <= kd; j++) {
        nc++; W[., nc] = Xd[., j]; nm[nc] = xdn[j]
    }
    for (j = 1; j <= kx; j++) {
        if (!iscont[j]) continue
        if (sq) {
            nc++; W[., nc] = Xc[., j] :^ 2; nm[nc] = xcn[j] + "^2"
        }
        if (cb) {
            nc++; W[., nc] = Xc[., j] :^ 3; nm[nc] = xcn[j] + "^3"
        }
    }
    if (ixx) {
        for (j = 1; j <= kx; j++) {
            for (u = j + 1; u <= kx; u++) {
                nc++; W[., nc] = Xc[., j] :* Xc[., u]; nm[nc] = xcn[j] + "#" + xcn[u]
            }
        }
    }
    if (dix) {
        for (i = 1; i <= kd; i++) {
            for (j = 1; j <= kx; j++) {
                nc++; W[., nc] = Xd[., i] :* Xc[., j]; nm[nc] = xdn[i] + "#" + xcn[j]
            }
        }
    }
    if (ddd) {
        for (i = 1; i <= kd; i++) {
            for (u = i + 1; u <= kd; u++) {
                if (dsrc[i] == dsrc[u]) continue
                nc++; W[., nc] = Xd[., i] :* Xd[., u]; nm[nc] = xdn[i] + "#" + xdn[u]
            }
        }
    }
    if (nc < cols(W)) {
        if (nc > 0) {
            W  = W[., (1..nc)]
            nm = nm[(1..nc)]
        }
        else {
            W  = J(n, 0, .)
            nm = J(1, 0, "")
        }
    }
    if (cols(W) > 0) {
        mu   = colsum(W) / n
        cv   = colsum((W :- mu) :^ 2) / n
        am   = abs(mu)
        am   = am :* (am :> 1) + (am :<= 1)
        keep = cv :> (1e-12 :* am :^ 2)
        if (!all(keep)) {
            W  = select(W, keep)
            nm = select(nm, keep)
        }
    }
    names = nm
    return(W)
}

end

// =============================================================================
//  The analysis (one call per estimation; also each bootstrap replicate)
// =============================================================================
mata:

struct _ivm_opt {
    string scalar    se              // classical robust cluster spatial
    real scalar      point_only      // bootstrap replicates: point estimates only
    real scalar      quiet           // no notes
    real scalar      weighted
    real scalar      alpha, optau, opalpha
    real scalar      doar, arpoints, arwidth
    real scalar      cutoff
    string scalar    kernel, metric
    real scalar      lasso, sq, cb, ixx, dix, ddd
    string rowvector keeppat
    real scalar      lniter, ltol, lcslack, lalphap, llambda
}

struct _ivm_dat {
    real colvector   y, t, m, w
    real matrix      Z, X            // LASSO: X = [continuous terms, dummy indicators]
    real scalar      kxc             // LASSO: number of continuous columns of X
    real rowvector   dsrc            // LASSO: factor id of each dummy column
    real scalar      ndum            // number of dummy indicator columns
    string rowvector xnames, znames
    string scalar    tname, mname
    real matrix      FE
    real colvector   cl, lat, lon
}

struct _ivm_res {
    real scalar      n, K, kcov, ndum, nsing, nfeg, dffe, dfvcov, G, zc
    real rowvector   nested
    string rowvector dropped, xsel
    real colvector   keep
    real scalar      theta, se_theta, tau, se_tau, gamma, se_gamma, total, se_total
    real scalar      indirect, se_indirect, cov_gt
    real rowvector   pi1, se_pi1
    real scalar      decomp_check, decomp_error, decomp_rel
    real scalar      F_TZ, F_TZ_op, K_eff_TZ, F_TZ_crit, F_MZ_T, F_MZ_T_op, K_eff_MZ_T, F_MZ_T_crit
    real colvector   F_ind
    real scalar      c_TM, sigma_eT, sigma_eY, kappa_star, kappa_star_CI, ks_inf, rho_TY_star, kt_rel
    real matrix      kadj
    real scalar      F_rf, F_fs, RV, XRV, rv_rf, rv_fs, xrv_rf, xrv_fs
    real scalar      ar_done, spneg, nb_share
    struct _ivm_ar scalar ar_theta, ar_gamma, ar_total
    real colvector   eT, eM, eY
    // LASSO
    real scalar      lasso_p, lasso_nkeep, lasso_nsel
    string rowvector rich_names
    real matrix      lasso_S, lasso_B
    real rowvector   lasso_lambda, lasso_niter, lasso_conv, lasso_union
    string rowvector lasso_method
}

// weighted variation around the weighted mean (reference scale of a column)
real rowvector _ivm_refnorm(real matrix V, real colvector w)
{
    real rowvector mu
    if (cols(V) == 0) return(J(1, 0, .))
    mu = colsum(V :* w) / sum(w)
    return(sqrt(colsum(w :* (V :- mu) :^ 2)))
}

string scalar _ivm_join(string rowvector s, string scalar sep)
{
    string scalar out
    real scalar   i
    out = ""
    for (i = 1; i <= cols(s); i++) out = out + (i > 1 ? sep : "") + s[i]
    return(out)
}

struct _ivm_res scalar _ivm_analysis(struct _ivm_dat scalar D0, struct _ivm_opt scalar O)
{
    struct _ivm_res scalar R
    struct _ivm_se  scalar S, SL
    struct _ivm_reg scalar rfs, rg, rth, rtot, rc, rrf, rfc
    struct _ivm_lres scalar LY, LT, LM, LZ
    real colvector y, t, m, w, sw, cl, lat, lon, keep, tw, mw, yw, r_, eT, tor, mor, sc_
    real matrix    Z, X, FE, V, Zw, Xw, exog, Zb, Xg, Xth, Zth, Xc, Wf, Xz, Qzz, Ztl, A
    real matrix    Wr, Xcon, Xdum, Xp, Zor, Sm, Bm
    real rowvector refc, refz, kept, zi, keepm, un, zi2, sel
    real scalar    n0, n, K, j, nx, refT, refM, kmax, dcr1, nG, d, nnz, share, iT, iM
    real scalar    vind, ftheta, x_, kst, dfr, df_rv, tstar, fstar, fb, fqa, i, ip, nm
    real scalar    sgn, nfrac, kf
    string rowvector xn, rn, pats, xcn_, xdn_
    string scalar  msg
    real rowvector fracs

    R.ar_done = 0
    R.spneg = 0
    R.nb_share = .
    R.lasso_p = 0
    R.lasso_nkeep = 0
    R.lasso_nsel = 0
    R.dropped = J(1, 0, "")
    R.rich_names = J(1, 0, "")
    n0 = rows(D0.y)
    K  = cols(D0.Z)
    R.K = K
    R.ndum = D0.ndum

    // ---- singleton FE groups (iteratively, every dimension) ----
    keep = _ivm_singletons(D0.FE)
    R.keep  = keep
    R.nsing = n0 - sum(keep)
    if (sum(keep) == 0) {
        _ivm_err(2000, "no observations remaining after dropping singleton groups: every absorb() group has a single observation")
    }
    if (R.nsing > 0 & !O.quiet) {
        printf("{txt}  Note: dropped %g singleton observations across %g FE dimension(s).\n", R.nsing, cols(D0.FE))
    }
    if (R.nsing > 0) {
        sel = selectindex(keep)
        y = D0.y[sel]; t = D0.t[sel]; m = D0.m[sel]; w = D0.w[sel]
        Z = D0.Z[sel, .]
        X = J(rows(sel), 0, .)
        FE = J(rows(sel), 0, .)
        cl = J(0, 1, .)
        lat = J(0, 1, .)
        lon = J(0, 1, .)
        if (cols(D0.X))  X  = D0.X[sel, .]
        if (cols(D0.FE)) FE = D0.FE[sel, .]
        if (rows(D0.cl))  cl  = D0.cl[sel]
        if (rows(D0.lat)) lat = D0.lat[sel]
        if (rows(D0.lon)) lon = D0.lon[sel]
    }
    else {
        y = D0.y; t = D0.t; m = D0.m; w = D0.w; Z = D0.Z; X = D0.X; FE = D0.FE
        cl = D0.cl; lat = D0.lat; lon = D0.lon
        if (cols(X) == 0) X = J(n0, 0, .)
        if (cols(FE) == 0) FE = J(n0, 0, .)
    }
    n = rows(y)
    R.n = n

    // ---- clusters ----
    R.G = .
    if (O.se == "cluster") {
        R.G = rows(uniqrows(cl))
        if (R.G < 2) _ivm_err(459, "cluster SEs need at least 2 clusters; the estimation sample has " + strofreal(R.G))
    }

    // ---- absorbed FE degrees of freedom ----
    R.dffe = 0
    R.dfvcov = 0
    R.nfeg = 0
    R.nested = J(1, 0, .)
    if (cols(FE) > 0) {
        _ivm_fedf(FE, cl, (O.se == "cluster"), d, dcr1, nG, R.nested)
        R.dffe = d
        R.nfeg = nG
        R.dfvcov = (O.se == "cluster" ? dcr1 : d)
        if (O.se != "cluster") R.nested = J(1, 0, .)
    }

    // ---- LASSO: rich expansion on level data, before FE demeaning ----
    xn = D0.xnames
    if (O.lasso) {
        Xcon = J(n, 0, .)
        Xdum = J(n, 0, .)
        xcn_ = J(1, 0, "")
        xdn_ = J(1, 0, "")
        if (D0.kxc > 0) {
            Xcon = X[., (1..D0.kxc)]
            xcn_ = xn[(1..D0.kxc)]
        }
        if (cols(X) > D0.kxc) {
            Xdum = X[., ((D0.kxc + 1)..cols(X))]
            xdn_ = xn[((D0.kxc + 1)..cols(xn))]
        }
        X = _ivm_rich(Xcon, Xdum, D0.dsrc, xcn_, xdn_, O.sq, O.cb, O.ixx, O.dix, O.ddd, rn)
        xn = rn
        R.rich_names = rn
        R.lasso_p = cols(X)
        // lkeep(): patterns resolved against the expansion names
        keepm = J(1, cols(X), 0)
        pats = O.keeppat
        for (i = 1; i <= cols(pats); i++) {
            nm = 0
            for (j = 1; j <= cols(rn); j++) {
                if (strmatch(rn[j], pats[i])) {
                    keepm[j] = 1
                    nm++
                }
            }
            if (nm == 0) {
                msg = "lkeep(): " + pats[i] + " matches no term of the LASSO expansion"
                if (cols(rn) > 0) msg = msg + "; available terms include: " + _ivm_join(rn[(1..min((8, cols(rn))))], " ")
                _ivm_err(198, msg)
            }
        }
        R.lasso_nkeep = sum(keepm)
        if (!O.quiet & cols(X) > 0) {
            printf("{txt}  Rich expansion built: p = %g\n", cols(X))
        }
    }

    // ---- reference scales (weighted variation before the FE are absorbed) ----
    refc = _ivm_refnorm(X, w)
    refT = _ivm_refnorm(t, w)
    refM = _ivm_refnorm(m, w)
    refz = _ivm_refnorm(Z, w)

    // ---- FE demeaning (weighted group means), all analysis columns together ----
    if (cols(FE) > 0) {
        V = _ivm_demean((y, t, m, Z, X), FE, w, O.quiet)
        y = V[., 1]; t = V[., 2]; m = V[., 3]
        Z = V[., (4..(3 + K))]
        if (cols(X)) X = V[., ((4 + K)..cols(V))]
    }

    // ---- variance regime ----
    S.type = (O.point_only ? "none" : O.se)
    S.G = R.G
    if (O.se == "cluster") {
        S.clo = order((cl, (1::n)), (1, 2))
        S.clinfo = panelsetup(cl[S.clo], 1)
    }
    if (O.se == "spatial" & (!O.point_only | O.lasso)) {
        _ivm_spatial(lat, lon, O.cutoff, O.kernel, O.metric, S.sri, S.scj, S.swt)
        nnz = rows(S.sri)
        share = nnz / n / n
        R.nb_share = share
        if (O.kernel == "uniform" & nnz >= n * n) {
            _ivm_err(459, "the spatial cutoff covers all pairs; Conley SEs are degenerate (0 by construction). Choose a smaller cutoff().")
        }
        if (share > 0.5 & !O.quiet) {
            printf("{txt}  Warning: the spatial cutoff puts %2.0f%% of all pairs within the cutoff (average neighbor share > 50%%);\n", 100 * share)
            printf("{txt}  Conley SEs are unreliable in this regime. Consider a smaller cutoff.\n")
        }
    }

    // ---- LASSO selection (four LASSOs on the FE-demeaned raw data) ----
    if (O.lasso) {
        SL = S
        if (O.se != "cluster" & O.se != "spatial") SL.type = "robust"
        else SL.type = O.se
        R.lasso_lambda = J(1, 3 + K, .)
        R.lasso_niter  = J(1, 3 + K, .)
        R.lasso_conv   = J(1, 3 + K, .)
        R.lasso_method = J(1, 3 + K, "")
        R.lasso_S = J(cols(X), 3 + K, 0)
        R.lasso_B = J(cols(X), 3 + K, 0)
        R.lasso_union = J(1, cols(X), 0)
        if (cols(X) > 0) {
            for (j = 1; j <= 3 + K; j++) {
                if (j == 1)      LY = _ivm_bch(y, X, w, SL, keepm, O.lniter, O.ltol, O.lcslack, O.lalphap, O.llambda)
                else if (j == 2) LY = _ivm_bch(t, X, w, SL, keepm, O.lniter, O.ltol, O.lcslack, O.lalphap, O.llambda)
                else if (j == 3) LY = _ivm_bch(m, X, w, SL, keepm, O.lniter, O.ltol, O.lcslack, O.lalphap, O.llambda)
                else             LY = _ivm_bch(Z[., j - 3], X, w, SL, keepm, O.lniter, O.ltol, O.lcslack, O.lalphap, O.llambda)
                R.lasso_S[., j] = LY.sel
                R.lasso_B[., j] = LY.beta
                R.lasso_lambda[j] = LY.lambda
                R.lasso_niter[j]  = LY.n_iter
                R.lasso_conv[j]   = LY.converged
                R.lasso_method[j] = LY.method
            }
            R.lasso_union = (rowsum(R.lasso_S) :> 0)'
            R.lasso_union = R.lasso_union :| keepm
            sel = selectindex(R.lasso_union)
            R.lasso_nsel = length(sel)
            if (length(sel) == 0) {
                X = J(n, 0, .)
                xn = J(1, 0, "")
                refc = J(1, 0, .)
            }
            else {
                X = X[., sel]
                xn = xn[sel]
                refc = refc[sel]
            }
            if (!O.quiet) {
                printf("{txt}  LASSO selected %g of %g candidate variables (plugin penalty).\n", R.lasso_nsel, R.lasso_p)
            }
        }
    }

    // ---- collinear controls: drop exactly redundant columns of [1, X] ----
    sw = sqrt(w)
    if (cols(X) > 0) {
        kept = _ivm_greedy(X :* sw, sw, refc, 1e-9)
        if (!all(kept)) {
            R.dropped = select(xn, !kept)
            if (!O.quiet) printf("{txt}  Note: dropped (collinear): %s\n", _ivm_join(R.dropped, ", "))
            sel = selectindex(kept)
            if (length(sel)) {
                X  = X[., sel]
                xn = xn[sel]
            }
            else {
                X  = J(n, 0, .)
                xn = J(1, 0, "")
            }
        }
    }
    R.xsel = xn
    R.kcov = cols(X)

    // ---- sqrt(w) premultiplication ----
    yw = y :* sw
    tw = t :* sw
    mw = m :* sw
    Zw = Z :* sw
    Xw = X :* sw
    exog = (sw, Xw)
    Zb   = (sw, Zw, Xw)
    R.zc = invnormal(1 - O.alpha / 2)
    nx = cols(exog)

    // ---- identification (paper Theorem 1, eq. 26) ----
    for (j = 1; j <= K; j++) {
        r_ = _ivm_resid(Zw[., j], exog)
        if (!(refz[j] > 0) | (r_' * r_) / refz[j]^2 < 1e-12) {
            _ivm_err(459, "instrument " + D0.znames[j] + " has no variation after the controls/fixed effects (sigma_ZZ = 0, paper eq. 26)")
        }
    }
    if (K > 1) {
        Ztl = _ivm_resid(Zw, exog)
        kept = _ivm_greedy(Ztl, J(n, 0, .), refz, 1e-9)
        if (!all(kept)) {
            _ivm_err(459, "instruments are collinear after the controls/fixed effects (" + _ivm_join(select(D0.znames, !kept), ", ") + ")")
        }
    }
    r_ = _ivm_resid(tw, exog)
    if (!(refT > 0) | (r_' * r_) / refT^2 < 1e-12) {
        _ivm_err(459, "treatment " + D0.tname + " has no variation after the controls/fixed effects")
    }
    r_ = _ivm_resid(mw, exog)
    if (!(refM > 0) | (r_' * r_) / refM^2 < 1e-12) {
        _ivm_err(459, "mediator " + D0.mname + " has no variation after the controls/fixed effects")
    }
    r_ = _ivm_resid(mw, (exog, tw))
    if ((r_' * r_) / refM^2 < 1e-12) {
        _ivm_err(459, "the mediator is collinear with the treatment after the controls/fixed effects")
    }
    kmax = nx + K + 1
    if (n - kmax - R.dffe < 1) {
        _ivm_err(2001, "not enough observations: n = " + strofreal(n) + ", but the largest regression has " +
                 strofreal(kmax) + " coefficients plus " + strofreal(R.dffe) +
                 " absorbed fixed-effect degrees of freedom (n must exceed their sum)")
    }

    // ---- 1. first stage T ~ Z | X: pi1 and F_{T,Z|X} ----
    rfs = _ivm_2sls(tw, Zb, Zb, 1, S, R.dfvcov)
    zi  = (2..(1 + K))
    R.pi1    = rfs.coef[zi]'
    R.se_pi1 = rfs.se[zi]'
    R.F_TZ = _ivm_waldK(rfs.coef[zi], rfs.vcov[zi, zi])
    Ztl = _ivm_resid(Zw, exog)
    Qzz = (Ztl' * Ztl) / n
    if (K == 1) R.F_TZ_op = R.F_TZ
    else {
        R.F_TZ_op = (_ivm_is_pd(rfs.vcov[zi, zi]) ?
            (rfs.coef[zi]' * Qzz * rfs.coef[zi]) / trace(Qzz * rfs.vcov[zi, zi]) : .)
    }
    _ivm_keff(rfs.vcov[zi, zi], Qzz, O.optau, O.opalpha, R.K_eff_TZ, R.F_TZ_crit)

    // ---- 2. gamma: M ~ T | Z, X ----
    Xg = (exog, tw)
    rg = _ivm_2sls(mw, Xg, Zb, 0, S, R.dfvcov)
    R.gamma = rg.coef[nx + 1]
    R.se_gamma = rg.se[nx + 1]

    // ---- 3. conditional 2SLS: Y ~ M + T | Z, T, X ----
    Xth = (exog, mw, tw)
    Zth = (Zb, tw)
    rth = _ivm_2sls(yw, Xth, Zth, 0, S, R.dfvcov)
    iM = nx + 1
    iT = nx + 2
    R.theta = rth.coef[iM]
    R.tau   = rth.coef[iT]
    R.se_theta = rth.se[iM]
    R.se_tau   = rth.se[iT]

    // ---- 4. total effect: Y ~ T | Z, X ----
    rtot = _ivm_2sls(yw, Xg, Zb, 0, S, R.dfvcov)
    R.total = rtot.coef[nx + 1]
    R.se_total = rtot.se[nx + 1]

    // ---- 5. conditional first stage M ~ Z | T, X ----
    Xc = (exog, Zw, tw)
    rc = _ivm_2sls(mw, Xc, Xc, 1, S, R.dfvcov)
    zi2 = ((nx + 1)..(nx + K))
    R.F_MZ_T = _ivm_waldK(rc.coef[zi2], rc.vcov[zi2, zi2])
    R.F_ind = J(K, 1, .)
    for (j = 1; j <= K; j++) {
        if (rc.vcov[nx + j, nx + j] > 0 & rc.vcov[nx + j, nx + j] < .) {
            R.F_ind[j] = rc.coef[nx + j]^2 / rc.vcov[nx + j, nx + j]
        }
    }
    Ztl = _ivm_resid(Zw, (exog, tw))
    Qzz = (Ztl' * Ztl) / n
    if (K == 1) R.F_MZ_T_op = R.F_MZ_T
    else {
        R.F_MZ_T_op = (_ivm_is_pd(rc.vcov[zi2, zi2]) ?
            (rc.coef[zi2]' * Qzz * rc.coef[zi2]) / trace(Qzz * rc.vcov[zi2, zi2]) : .)
    }
    _ivm_keff(rc.vcov[zi2, zi2], Qzz, O.optau, O.opalpha, R.K_eff_MZ_T, R.F_MZ_T_crit)

    // ---- 6. decomposition check and the delta-method indirect effect ----
    R.indirect = R.theta * R.gamma
    R.decomp_check = R.tau + R.indirect
    R.decomp_error = abs(R.decomp_check - R.total)
    R.decomp_rel = R.decomp_error / max((abs(R.total), abs(R.tau) + abs(R.indirect), 1e-300))
    R.cov_gt = _ivm_xcov(rg, rth, nx + 1, iM, S, R.dfvcov)
    vind = R.gamma^2 * R.se_theta^2 + R.theta^2 * R.se_gamma^2 + 2 * R.theta * R.gamma * R.cov_gt
    if (vind < 0) {
        if (vind > -1e-8 * max(((R.gamma * R.se_theta)^2, 1e-300))) vind = 0
        else {
            if (!O.quiet) {
                printf("{txt}  Warning: delta-method variance for the indirect effect is negative (%9.2e). This typically\n", vind)
                printf("{txt}  reflects a poorly identified first stage. Reporting missing; consider the bootstrap.\n")
            }
            vind = .
        }
    }
    R.se_indirect = (vind < . ? sqrt(vind) : .)

    // ---- 7. kappa-based sensitivity (Hajek-weighted moments, raw scale) ----
    Xp = (J(n, 1, 1), X)
    tor = _ivm_wls_resid(Xp, t, w)
    mor = _ivm_wls_resid(Xp, m, w)
    Zor = J(n, K, .)
    for (j = 1; j <= K; j++) Zor[., j] = _ivm_wls_resid(Xp, Z[., j], w)
    eT = _ivm_wls_resid((J(n, 1, 1), Zor), tor, w)
    R.eT = eT
    R.eM = mor - R.gamma * tor
    R.eY = rth.resid :/ sw
    R.c_TM = _ivm_wcov(R.eT, R.eM, w)
    R.sigma_eT = sqrt(_ivm_wcov(R.eT, R.eT, w))
    R.sigma_eY = sqrt(_ivm_wcov(R.eY, R.eY, w))
    R.kappa_star = R.theta * R.c_TM
    sgn = sign(R.theta)
    R.kappa_star_CI = (R.theta - sgn * R.zc * R.se_theta) * R.c_TM
    R.ks_inf = (R.se_theta < . ? (abs(R.theta) > R.zc * R.se_theta) : .)
    if (R.sigma_eT > 0 & R.sigma_eY > 0 & R.sigma_eT < . & R.sigma_eY < .) {
        R.rho_TY_star = R.kappa_star / (R.sigma_eT * R.sigma_eY)
    }
    else R.rho_TY_star = .
    fracs = (0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00)
    if (R.kappa_star >= . | R.kappa_star == 0) fracs = (0)
    R.kadj = J(cols(fracs), 5, .)
    for (i = 1; i <= cols(fracs); i++) {
        kst = (cols(fracs) == 1 ? 0 : fracs[i] * R.kappa_star)
        ftheta = R.theta - kst / R.c_TM
        R.kadj[i, .] = (kst, fracs[i], ftheta, ftheta - R.zc * R.se_theta, ftheta + R.zc * R.se_theta)
    }
    R.kt_rel = ((R.F_MZ_T_op < . & R.F_MZ_T_crit < .) ? (R.F_MZ_T_op >= R.F_MZ_T_crit) : .)

    // ---- 8. R^2-based sensitivity (first instrument) ----
    Wf = (exog, tw)
    Sm = _ivm_resid((yw, mw, Zw[., 1]), Wf)
    Xz = (sw, Sm[., 3])
    dfr = (cols(Wf) - 1) + R.dfvcov
    if (!O.point_only) {
        rrf = _ivm_2sls(Sm[., 1], Xz, Xz, 1, S, dfr)
        rfc = _ivm_2sls(Sm[., 2], Xz, Xz, 1, S, dfr)
        R.F_rf = (rrf.coef[2] / rrf.se[2])^2
        R.F_fs = (rfc.coef[2] / rfc.se[2])^2
    }
    else {
        R.F_rf = .
        R.F_fs = .
    }
    df_rv = n - cols(Wf) - 1 - R.dffe
    tstar = abs(invttail(df_rv - 1, O.alpha / 2))
    fstar = tstar / sqrt(df_rv - 1)
    if (R.F_rf < . & R.F_fs < .) {
        R.xrv_rf = max(((abs(R.F_rf) / df_rv - fstar^2) / (1 + abs(R.F_rf) / df_rv), 0))
        R.xrv_fs = max(((abs(R.F_fs) / df_rv - fstar^2) / (1 + abs(R.F_fs) / df_rv), 0))
        fb = sqrt(abs(R.F_rf) / df_rv)
        fqa = fb - fstar
        if (fqa < 0) R.rv_rf = 0
        else if (fb < 1 / fstar) R.rv_rf = 0.5 * (sqrt(fqa^4 + 4 * fqa^2) - fqa^2)
        else R.rv_rf = R.xrv_rf
        fb = sqrt(abs(R.F_fs) / df_rv)
        fqa = fb - fstar
        if (fqa < 0) R.rv_fs = 0
        else if (fb < 1 / fstar) R.rv_fs = 0.5 * (sqrt(fqa^4 + 4 * fqa^2) - fqa^2)
        else R.rv_fs = R.xrv_fs
        R.XRV = min((R.xrv_rf, R.xrv_fs))
        R.RV  = min((R.rv_rf, R.rv_fs))
    }
    else {
        R.xrv_rf = .; R.xrv_fs = .; R.rv_rf = .; R.rv_fs = .; R.XRV = .; R.RV = .
    }

    // ---- 9. Anderson-Rubin sets (same SE regime and alpha) ----
    if (O.doar & !O.point_only) {
        R.ar_done = 1
        R.ar_theta = _ivm_arset(yw, mw, (exog, tw, Zw), ((nx + 2)..(nx + 1 + K)), S, R.dfvcov,
                                O.alpha, R.theta, R.se_theta, O.arwidth, O.arpoints)
        R.ar_gamma = _ivm_arset(mw, tw, Zb, (2..(1 + K)), S, R.dfvcov,
                                O.alpha, R.gamma, R.se_gamma, O.arwidth, O.arpoints)
        R.ar_total = _ivm_arset(yw, tw, Zb, (2..(1 + K)), S, R.dfvcov,
                                O.alpha, R.total, R.se_total, O.arwidth, O.arpoints)
    }

    // ---- negative spatial variances: one warning, reported quantities only ----
    if (O.se == "spatial" & !O.point_only) {
        if (rg.vcov[nx + 1, nx + 1] < 0 | rth.vcov[iM, iM] < 0 | rth.vcov[iT, iT] < 0 |
            rtot.vcov[nx + 1, nx + 1] < 0 | rrf.vcov[2, 2] < 0 | rfc.vcov[2, 2] < 0 |
            !_ivm_is_pd(rfs.vcov[zi, zi]) | !_ivm_is_pd(rc.vcov[zi2, zi2])) {
            R.spneg = 1
            if (!O.quiet) {
                printf("{txt}  Warning: the spatial (Conley) variance was negative for some reported coefficients (the kernel\n")
                printf("{txt}  matrix is not positive semi-definite); their SEs and the statistics that use them are missing.\n")
            }
        }
    }
    return(R)
}

end

// =============================================================================
//  Stata interface: data in, bootstrap, e() out
// =============================================================================
mata:

struct _ivm_state {
    struct _ivm_dat scalar D          // estimation sample before singleton removal
    struct _ivm_opt scalar O
    struct _ivm_res scalar R
    // bootstrap
    string scalar    bmethod, bcitype
    real scalar      B, need_se, diag, has_idx, n_est, Gb, nfailed
    real matrix      idx
    real colvector   corder, cstart, csize      // block: rows grouped by cluster (first appearance)
    real rowvector   fenest, fenlev
    real matrix      bpoint, bse, bci, bana
    string rowvector pnames
}

void _ivm_put(string scalar nm, transmorphic x)
{
    pointer() scalar p
    p = findexternal(nm)
    if (p == NULL) p = crexternal(nm)
    *p = x
}

transmorphic _ivm_get(string scalar nm)
{
    pointer() scalar p
    p = findexternal(nm)
    if (p == NULL) _ivm_err(3499, "ivmediate: internal state " + nm + " not found")
    return(*p)
}

void _ivmediate_main()
{
    struct _ivm_dat   scalar D
    struct _ivm_opt   scalar O
    struct _ivm_res   scalar R
    struct _ivm_state scalar G
    string scalar    touse, s, expx, expd
    string rowvector xv, xc, xd, xcn, xdn, fe, srcs
    real colvector   obs
    real scalar      n, j, k, lv, ll

    touse = st_local("touse")
    D.y = st_data(., st_local("yvar"), touse)
    n = rows(D.y)
    D.t = st_data(., st_local("treatment"), touse)
    D.m = st_data(., st_local("mediator"), touse)
    D.znames = tokens(st_local("instrument"))
    D.Z = st_data(., D.znames, touse)
    D.tname = st_local("treatment")
    D.mname = st_local("mediator")
    D.ndum = strtoreal(st_local("ndum"))
    O.lasso = strtoreal(st_local("lassoflag"))
    if (!O.lasso) {
        xv = tokens(st_local("xvars"))
        D.xnames = tokens(st_local("xnames"))
        D.X = J(n, 0, .)
        if (cols(xv)) D.X = st_data(., xv, touse)
        D.kxc = cols(xv)
        D.dsrc = J(1, 0, .)
    }
    else {
        xc  = tokens(st_local("xc_vars"))
        xd  = tokens(st_local("xd_vars"))
        xcn = tokens(st_local("xc_names"))
        xdn = tokens(st_local("xd_names"))
        for (j = 1; j <= cols(xcn); j++) {
            if (substr(xcn[j], 1, 2) == "c.") xcn[j] = substr(xcn[j], 3, .)
        }
        D.X = J(n, 0, .)
        if (cols(xc)) D.X = st_data(., xc, touse)
        if (cols(xd)) D.X = (D.X, st_data(., xd, touse))
        D.kxc = cols(xc)
        D.xnames = (xcn, xdn)
        srcs = J(1, 0, "")
        D.dsrc = J(1, cols(xdn), .)
        for (j = 1; j <= cols(xdn); j++) {
            s = substr(xdn[j], strpos(xdn[j], ".") + 1, .)
            k = 0
            for (lv = 1; lv <= cols(srcs); lv++) if (srcs[lv] == s) k = lv
            if (k == 0) {
                srcs = (srcs, s)
                k = cols(srcs)
            }
            D.dsrc[j] = k
        }
    }
    fe = tokens(st_local("feids"))
    D.FE = J(n, 0, .)
    if (cols(fe)) D.FE = st_data(., fe, touse)
    D.cl = J(0, 1, .)
    if (st_local("clid") != "") D.cl = st_data(., st_local("clid"), touse)
    D.w = J(n, 1, 1)
    O.weighted = (st_local("wtype") != "")
    if (O.weighted) D.w = st_data(., st_local("wv"), touse)
    O.se = st_local("setype")
    D.lat = J(0, 1, .)
    D.lon = J(0, 1, .)
    if (O.se == "spatial") {
        D.lat = st_data(., st_local("latitude"), touse)
        D.lon = st_data(., st_local("longitude"), touse)
        O.cutoff = strtoreal(st_local("cutoff"))
        O.kernel = st_local("spkern")
        O.metric = st_local("spmetric")
    }
    O.point_only = 0
    O.quiet = 0
    O.alpha = (100 - strtoreal(st_local("level"))) / 100
    O.optau = strtoreal(st_local("optau"))
    O.opalpha = strtoreal(st_local("opalpha"))
    O.doar = strtoreal(st_local("doar"))
    O.arpoints = strtoreal(st_local("arpoints"))
    O.arwidth = strtoreal(st_local("arwidth"))
    if (O.lasso) {
        expx = st_local("expx")
        expd = st_local("expd")
        O.sq  = (strpos(" " + expx + " ", " square ") > 0)
        O.cb  = (strpos(" " + expx + " ", " cube ") > 0)
        O.ixx = (strpos(" " + expx + " ", " interact_with_X ") > 0)
        O.dix = (strpos(" " + expd + " ", " interact_with_X ") > 0)
        O.ddd = (strpos(" " + expd + " ", " interact_with_dummies ") > 0)
        // lkeep(): term names carry no c. prefix (c.turn^2 is the term turn^2)
        O.keeppat = tokens(st_local("lkeep"))
        for (j = 1; j <= cols(O.keeppat); j++) {
            O.keeppat[j] = subinstr(O.keeppat[j], "#c.", "#")
            if (substr(O.keeppat[j], 1, 2) == "c.") O.keeppat[j] = substr(O.keeppat[j], 3, .)
        }
        O.lniter  = strtoreal(st_local("lniter"))
        O.ltol    = strtoreal(st_local("ltol"))
        O.lcslack = strtoreal(st_local("lcslack"))
        ll = strtoreal(st_local("lalphap"))
        O.lalphap = (ll == -1 ? . : ll)
        ll = strtoreal(st_local("llambda"))
        O.llambda = (ll == -1 ? . : ll)
    }
    else O.keeppat = J(1, 0, "")

    R = _ivm_analysis(D, O)

    obs = selectindex(st_data(., touse))
    st_store(obs[selectindex(R.keep)], st_local("keep"), J(sum(R.keep), 1, 1))
    G.D = D
    G.O = O
    G.R = R
    G.B = 0
    _ivm_put("_ivmediate_R", G)
}

// ---- e(b), e(V): total, direct and indirect effect ------------------------------
void _ivmediate_bV(string scalar bnm, string scalar Vnm)
{
    struct _ivm_state scalar G
    real rowvector b, s
    string matrix  cs
    real scalar    vu

    G = _ivm_get("_ivmediate_R")
    b = (G.R.total, G.R.tau, G.R.indirect)
    s = (G.R.se_total, G.R.se_tau, G.R.se_indirect)
    vu = missing(s)
    cs = (J(3, 1, ""), ("total_effect" \ "direct_effect" \ "indirect_effect"))
    st_matrix(bnm, b)
    st_matrixcolstripe(bnm, cs)
    if (!vu) {
        st_matrix(Vnm, diag(s :^ 2))
        st_matrixcolstripe(Vnm, cs)
        st_matrixrowstripe(Vnm, cs)
    }
    st_numscalar("__ivm_vunavail", vu)
}

void _ivm_es(string scalar nm, real scalar v)
{
    st_numscalar("e(" + nm + ")", v)
}

void _ivm_post_ar(string scalar pf, struct _ivm_ar scalar a)
{
    string matrix cs
    st_global("e(ar_" + pf + "_shape)", a.shape)
    st_global("e(ar_" + pf + "_method)", a.method)
    _ivm_es("ar_" + pf + "_stat0", a.stat0)
    _ivm_es("ar_" + pf + "_p", a.pval0)
    _ivm_es("ar_" + pf + "_contains0", a.contains0)
    _ivm_es("ar_" + pf + "_lower", a.lower)
    _ivm_es("ar_" + pf + "_upper", a.upper)
    _ivm_es("ar_" + pf + "_left_upper", a.left_upper)
    _ivm_es("ar_" + pf + "_right_lower", a.right_lower)
    if (a.shape == "grid") {
        st_matrix("e(ar_" + pf + "_seg)", a.seg)
        cs = (J(4, 1, ""), ("lower" \ "upper" \ "open_left" \ "open_right"))
        st_matrixcolstripe("e(ar_" + pf + "_seg)", cs)
    }
}

void _ivmediate_post()
{
    struct _ivm_state scalar G
    struct _ivm_res   scalar R
    string matrix    cs
    string rowvector lab, nm
    real scalar      K, zc, j

    G = _ivm_get("_ivmediate_R")
    R = G.R
    K = R.K
    zc = R.zc
    _ivm_es("k_inst", K)
    _ivm_es("k_cov", R.kcov)
    _ivm_es("k_dummies", R.ndum)
    _ivm_es("df_fe", R.dffe)
    _ivm_es("df_fe_vce", R.dfvcov)
    _ivm_es("N_fegroups", R.nfeg)
    _ivm_es("N_singleton", R.nsing)
    if (G.O.se == "cluster") _ivm_es("N_clust", R.G)
    if (cols(R.nested)) {
        st_matrix("e(fe_nested)", R.nested)
        st_matrixcolstripe("e(fe_nested)", (J(cols(R.nested), 1, ""), tokens(st_local("absorb"))'))
    }
    if (cols(R.dropped)) st_global("e(dropped_collinear)", _ivm_join(R.dropped, " "))
    _ivm_es("theta", R.theta);       _ivm_es("se_theta", R.se_theta)
    _ivm_es("tau", R.tau);           _ivm_es("se_tau", R.se_tau)
    _ivm_es("gamma", R.gamma);       _ivm_es("se_gamma", R.se_gamma)
    _ivm_es("total", R.total);       _ivm_es("se_total", R.se_total)
    _ivm_es("indirect", R.indirect); _ivm_es("se_indirect", R.se_indirect)
    _ivm_es("cov_gamma_theta", R.cov_gt)
    _ivm_es("z_crit", zc)
    _ivm_es("decomp_check", R.decomp_check)
    _ivm_es("decomp_error", R.decomp_error)
    _ivm_es("decomp_error_rel", R.decomp_rel)
    if (R.total != 0 & R.total < .) _ivm_es("mepct", 100 * R.indirect / R.total)
    st_matrix("e(pi1)", R.pi1)
    st_matrixcolstripe("e(pi1)", (J(K, 1, ""), G.D.znames'))
    st_matrix("e(se_pi1)", R.se_pi1)
    st_matrixcolstripe("e(se_pi1)", (J(K, 1, ""), G.D.znames'))
    _ivm_es("F_TZ", R.F_TZ)
    _ivm_es("F_TZ_op", R.F_TZ_op)
    _ivm_es("K_eff_TZ", R.K_eff_TZ)
    _ivm_es("F_TZ_crit", R.F_TZ_crit)
    _ivm_es("F_MZ_T", R.F_MZ_T)
    _ivm_es("F_MZ_T_op", R.F_MZ_T_op)
    _ivm_es("K_eff_MZ_T", R.K_eff_MZ_T)
    _ivm_es("F_MZ_T_crit", R.F_MZ_T_crit)
    st_matrix("e(F_MZ_T_individual)", R.F_ind)
    st_matrixrowstripe("e(F_MZ_T_individual)", (J(K, 1, ""), G.D.znames'))
    _ivm_es("c_TM", R.c_TM)
    _ivm_es("sigma_eT", R.sigma_eT)
    _ivm_es("sigma_eY", R.sigma_eY)
    _ivm_es("kappa_star", R.kappa_star)
    _ivm_es("kappa_star_CI", R.kappa_star_CI)
    _ivm_es("kappa_star_CI_informative", R.ks_inf)
    _ivm_es("rho_TY_star", R.rho_TY_star)
    _ivm_es("F_MZ_T_above_crit", R.kt_rel)
    st_matrix("e(kappa_adj)", R.kadj)
    st_matrixcolstripe("e(kappa_adj)", (J(5, 1, ""), ("kappa0" \ "fraction" \ "theta_adj" \ "ci_lo" \ "ci_hi")))
    _ivm_es("F_reduced_form", R.F_rf)
    _ivm_es("F_first_stage_cond", R.F_fs)
    _ivm_es("f_reduced_form", R.F_rf)
    _ivm_es("f_first_stage_cond", R.F_fs)
    _ivm_es("RV_theta", R.RV)
    _ivm_es("XRV_theta", R.XRV)
    _ivm_es("rv_rf", R.rv_rf)
    _ivm_es("rv_fs", R.rv_fs)
    _ivm_es("xrv_rf", R.xrv_rf)
    _ivm_es("xrv_fs", R.xrv_fs)
    _ivm_es("spatial_negvar", R.spneg)
    if (G.O.se == "spatial") _ivm_es("spatial_share", R.nb_share)
    _ivm_es("ar", R.ar_done)
    if (R.ar_done) {
        _ivm_es("ar_level", 100 * (1 - G.O.alpha))
        _ivm_es("ar_df", R.ar_theta.df)
        _ivm_es("ar_cv", R.ar_theta.cv)
        _ivm_es("ar_grid_points", G.O.arpoints)
        _ivm_es("ar_grid_halfwidth", G.O.arwidth)
        _ivm_post_ar("theta", R.ar_theta)
        _ivm_post_ar("gamma", R.ar_gamma)
        _ivm_post_ar("total", R.ar_total)
    }
    _ivm_es("lasso", G.O.lasso)
    if (G.O.lasso) {
        _ivm_es("lasso_n_total", R.lasso_p)
        _ivm_es("lasso_n_selected", R.lasso_nsel)
        _ivm_es("lasso_n_keep", R.lasso_nkeep)
        st_global("e(lasso_names)", _ivm_join(R.rich_names, " "))
        if (R.lasso_p > 0) {
            nm = select(R.rich_names, R.lasso_union)
            st_global("e(lasso_selected)", _ivm_join(nm, " "))
            if (cols(nm)) st_matrix("e(lasso_sel_idx)", selectindex(R.lasso_union))
            lab = ("Y", "T", "M")
            for (j = 1; j <= K; j++) lab = (lab, "Z" + strofreal(j))
            st_matrix("e(lasso_S)", R.lasso_S)
            st_matrixcolstripe("e(lasso_S)", (J(3 + K, 1, ""), lab'))
            st_matrix("e(lasso_B)", R.lasso_B)
            st_matrixcolstripe("e(lasso_B)", (J(3 + K, 1, ""), lab'))
            st_matrix("e(lasso_lambda)", R.lasso_lambda)
            st_matrixcolstripe("e(lasso_lambda)", (J(3 + K, 1, ""), lab'))
            st_matrix("e(lasso_n_iter)", R.lasso_niter)
            st_matrixcolstripe("e(lasso_n_iter)", (J(3 + K, 1, ""), lab'))
            st_matrix("e(lasso_converged)", R.lasso_conv)
            st_matrixcolstripe("e(lasso_converged)", (J(3 + K, 1, ""), lab'))
            st_global("e(lasso_methods)", _ivm_join(R.lasso_method, " "))
        }
    }
    st_global("e(controls)", _ivm_join(R.xsel, " "))
    _ivm_es("bootstrap", G.B > 0)
    if (G.B > 0) {
        _ivm_es("boot_B", G.B)
        _ivm_es("boot_n_failed", G.nfailed)
        _ivm_es("boot_n_est", G.n_est)
        _ivm_es("boot_alpha", G.O.alpha)
        cs = (J(cols(G.pnames), 1, ""), G.pnames')
        st_matrix("e(boot_point)", G.bpoint)
        st_matrixrowstripe("e(boot_point)", cs)
        st_matrix("e(boot_se)", G.bse)
        st_matrixrowstripe("e(boot_se)", cs)
        st_matrix("e(boot_ci)", G.bci)
        st_matrixrowstripe("e(boot_ci)", cs)
        st_matrixcolstripe("e(boot_ci)", (J(2, 1, ""), ("lower" \ "upper")))
        st_matrix("e(boot_analytic_se)", G.bana)
        st_matrixrowstripe("e(boot_analytic_se)", cs)
        st_matrix("e(boot_dist)", _ivm_get("_ivmediate_Bdist"))
        st_matrixcolstripe("e(boot_dist)", cs)
        _ivm_es("boot_se_total", G.bse[1])
        _ivm_es("boot_se_tau", G.bse[2])
        _ivm_es("boot_se_theta", G.bse[3])
        _ivm_es("boot_se_gamma", G.bse[4])
        _ivm_es("boot_se_indirect", G.bse[5])
    }
}

// residual variables stub_eT, stub_eM, stub_eY on the estimation sample
void _ivmediate_resid(string scalar stub, string scalar es)
{
    struct _ivm_state scalar G
    real colvector obs
    real scalar    v

    G = _ivm_get("_ivmediate_R")
    obs = selectindex(st_data(., es))
    v = st_addvar("double", stub + "_eT")
    st_store(obs, v, G.R.eT)
    v = st_addvar("double", stub + "_eM")
    st_store(obs, v, G.R.eM)
    v = st_addvar("double", stub + "_eY")
    st_store(obs, v, G.R.eY)
}

// ---- bootstrap -------------------------------------------------------------------
void _ivmediate_boot_init()
{
    struct _ivm_state scalar G
    real colvector cl, o, ids, fl, gfirst, cord
    real matrix    info
    real scalar    n, k, G0, g, mx, B, j
    string scalar  bi

    G = _ivm_get("_ivmediate_R")
    G.bmethod = st_local("bmethod")
    G.bcitype = st_local("bcitype")
    G.diag    = strtoreal(st_local("bdiag"))
    G.need_se = (G.bcitype == "studentized" | G.diag)
    n = rows(G.D.y)
    G.n_est = n
    G.pnames = ("total", "tau", "theta", "gamma", "indirect", "kappa_star", "rho_TY_star")
    if (G.diag) G.pnames = (G.pnames, "RV_theta", "XRV_theta")
    if (G.bmethod == "block") {
        // clusters in order of first appearance in the estimation sample
        cl  = G.D.cl
        o   = order((cl, (1::n)), (1, 2))
        info = panelsetup(cl[o], 1)
        G0  = rows(info)
        gfirst = o[info[., 1]]
        cord = order(gfirst, 1)
        G.cstart = J(G0, 1, .)
        G.csize  = J(G0, 1, .)
        G.corder = J(n, 1, .)
        k = 0
        for (g = 1; g <= G0; g++) {
            j = cord[g]
            G.cstart[g] = k + 1
            G.csize[g]  = info[j, 2] - info[j, 1] + 1
            G.corder[|k + 1 \ k + G.csize[g]|] = sort(o[|info[j, 1] \ info[j, 2]|], 1)
            k = k + G.csize[g]
        }
        G.Gb = G0
        // FE nested in the clusters get one level per drawn copy
        G.fenest = J(1, cols(G.D.FE), 0)
        G.fenlev = J(1, cols(G.D.FE), 0)
        for (j = 1; j <= cols(G.D.FE); j++) {
            mx = rows(uniqrows(G.D.FE[., j]))
            G.fenest[j] = (rows(uniqrows((G.D.FE[., j], cl))) == mx)
            G.fenlev[j] = max(G.D.FE[., j])
        }
    }
    bi = st_local("bootidx")
    G.has_idx = (bi != "")
    if (G.has_idx) {
        G.idx = st_matrix(bi)
        B = rows(G.idx)
        k = (G.bmethod == "block" ? G.Gb : n)
        if (cols(G.idx) != k) {
            _ivm_err(503, "bootidx(): matrix " + bi + " has " + strofreal(cols(G.idx)) + " columns; it needs " +
                     strofreal(k) + (G.bmethod == "block" ? " (G, the clusters of the estimation sample)" :
                     " (n, the size of the estimation sample)"))
        }
        if (B < 1) _ivm_err(503, "bootidx(): matrix " + bi + " has no rows")
        if (missing(G.idx) | any(G.idx :< 1) | any(G.idx :> k) | any(G.idx :!= floor(G.idx))) {
            _ivm_err(125, "bootidx(): entries must be integers between 1 and " + strofreal(k))
        }
    }
    else B = strtoreal(st_local("breps"))
    G.B = B
    st_numscalar("__ivm_bootB", B)
    _ivm_put("_ivmediate_Bdist", J(B, cols(G.pnames), .))
    _ivm_put("_ivmediate_Bse", J(B, 5, .))
    _ivm_put("_ivmediate_Bok", J(B, 1, 0))
    _ivm_put("_ivmediate_R", G)
}

void _ivmediate_boot_rep(real scalar b)
{
    struct _ivm_state scalar G
    struct _ivm_dat   scalar Db
    struct _ivm_opt   scalar Ob
    struct _ivm_res   scalar Rb
    real colvector    idx, g, rows_, copy, fev
    real scalar       n, G0, i, k, j, nb
    real rowvector    vals
    pointer() scalar  p

    G = _ivm_get("_ivmediate_R")
    n = rows(G.D.y)
    if (G.bmethod == "pairs") {
        if (G.has_idx) idx = G.idx[b, .]'
        else           idx = runiformint(n, 1, 1, n)
        rows_ = idx
        copy  = J(0, 1, .)
    }
    else {
        G0 = G.Gb
        if (G.has_idx) g = G.idx[b, .]'
        else           g = runiformint(G0, 1, 1, G0)
        nb = sum(G.csize[g])
        rows_ = J(nb, 1, .)
        copy  = J(nb, 1, .)
        k = 0
        for (i = 1; i <= G0; i++) {
            j = g[i]
            rows_[|k + 1 \ k + G.csize[j]|] = G.corder[|G.cstart[j] \ G.cstart[j] + G.csize[j] - 1|]
            copy[|k + 1 \ k + G.csize[j]|]  = J(G.csize[j], 1, i)
            k = k + G.csize[j]
        }
    }
    Db = G.D
    Db.y = G.D.y[rows_]
    Db.t = G.D.t[rows_]
    Db.m = G.D.m[rows_]
    Db.w = G.D.w[rows_]
    Db.Z = G.D.Z[rows_, .]
    Db.X = J(rows(rows_), 0, .)
    if (cols(G.D.X)) Db.X = G.D.X[rows_, .]
    Db.FE = J(rows(rows_), 0, .)
    if (cols(G.D.FE)) {
        Db.FE = G.D.FE[rows_, .]
        if (G.bmethod == "block") {
            for (j = 1; j <= cols(Db.FE); j++) {
                if (G.fenest[j]) Db.FE[., j] = copy :* (G.fenlev[j] + 1) + Db.FE[., j]
            }
        }
    }
    if (rows(G.D.cl)) {
        if (G.bmethod == "block") Db.cl = copy
        else                      Db.cl = G.D.cl[rows_]
    }
    if (rows(G.D.lat)) {
        Db.lat = G.D.lat[rows_]
        Db.lon = G.D.lon[rows_]
    }
    Ob = G.O
    Ob.point_only = !G.need_se
    Ob.quiet = 1
    Ob.doar = 0
    Rb = _ivm_analysis(Db, Ob)
    vals = (Rb.total, Rb.tau, Rb.theta, Rb.gamma, Rb.indirect, Rb.kappa_star, Rb.rho_TY_star)
    if (G.diag) vals = (vals, Rb.RV, Rb.XRV)
    p = findexternal("_ivmediate_Bdist")
    (*p)[b, .] = vals
    p = findexternal("_ivmediate_Bse")
    (*p)[b, .] = (Rb.se_total, Rb.se_tau, Rb.se_theta, Rb.se_gamma, Rb.se_indirect)
    p = findexternal("_ivmediate_Bok")
    (*p)[b] = 1
}

void _ivmediate_boot_finish()
{
    struct _ivm_state scalar G
    real matrix    dist, ser
    real colvector ok, v, tb, se_b, msk
    real rowvector pt, ana
    real scalar    P, B, p, q_lo, q_hi, a, nf, kse

    G = _ivm_get("_ivmediate_R")
    dist = _ivm_get("_ivmediate_Bdist")
    ser  = _ivm_get("_ivmediate_Bse")
    ok   = _ivm_get("_ivmediate_Bok")
    B = G.B
    P = cols(G.pnames)
    a = G.O.alpha
    nf = B - sum(ok)
    G.nfailed = nf
    if (nf >= B) {
        _ivm_err(459, "all " + strofreal(B) + " bootstrap replicates failed; check that the model specification is consistent with the data")
    }
    if (nf > 0.10 * B) {
        printf("{txt}  Warning: %g/%g bootstrap replicates failed (%4.1f%%). Bootstrap CIs may be unreliable.\n", nf, B, 100 * nf / B)
    }
    pt  = (G.R.total, G.R.tau, G.R.theta, G.R.gamma, G.R.indirect, G.R.kappa_star, G.R.rho_TY_star)
    ana = (G.R.se_total, G.R.se_tau, G.R.se_theta, G.R.se_gamma, G.R.se_indirect, ., .)
    if (G.diag) {
        pt  = (pt, G.R.RV, G.R.XRV)
        ana = (ana, ., .)
    }
    G.bpoint = pt'
    G.bana   = ana'
    G.bse = J(P, 1, .)
    G.bci = J(P, 2, .)
    for (p = 1; p <= P; p++) {
        v = select(dist[., p], dist[., p] :< .)
        G.bse[p] = _ivm_sd1(dist[., p])
        if (rows(v) == 0) continue
        if (G.bcitype == "percentile") {
            G.bci[p, .] = (_ivm_quantile(v, a / 2), _ivm_quantile(v, 1 - a / 2))
        }
        else if (G.bcitype == "studentized" & p <= 5) {
            se_b = ser[., p]
            msk  = (dist[., p] :< .) :& (se_b :< .)
            if (sum(msk) < 10) continue
            tb = (select(dist[., p], msk) :- pt[p]) :/ select(se_b, msk)
            tb = select(tb, tb :< .)
            if (rows(tb) < 10) continue
            q_lo = _ivm_quantile(tb, a / 2)
            q_hi = _ivm_quantile(tb, 1 - a / 2)
            if (ana[p] < .) G.bci[p, .] = (pt[p] - q_hi * ana[p], pt[p] - q_lo * ana[p])
        }
        else {
            q_lo = _ivm_quantile(v, a / 2)
            q_hi = _ivm_quantile(v, 1 - a / 2)
            G.bci[p, .] = (2 * pt[p] - q_hi, 2 * pt[p] - q_lo)
        }
    }
    _ivm_put("_ivmediate_R", G)
}

end

// =============================================================================
//  Display
// =============================================================================
mata:

// SE regime as typed in Stata's vce()
string scalar _ivm_selab(string scalar se)
{
    if (se == "classical") return("unadjusted")
    if (se == "spatial")   return("conley")
    return(se)
}

string scalar _ivm_version()
{
    return("2.0.0")
}

// fixed-width number: 4 decimals, or scientific notation when it would not fit
string scalar _ivm_num(real scalar x, real scalar w)
{
    if (x >= .) return(sprintf("%" + strofreal(w) + "s", "."))
    if (x == 0) return(sprintf("%" + strofreal(w) + ".4f", 0))
    if (abs(x) >= 1e-4 & abs(x) < 10^(w - 6)) return(sprintf("%" + strofreal(w) + ".4f", x))
    return(sprintf("%" + strofreal(w) + ".3e", x))
}

string scalar _ivm_numc(real scalar x)
{
    return(strtrim(_ivm_num(x, 12)))
}

string scalar _ivm_f2(real scalar x)
{
    if (x >= .) return(".")
    if (abs(x) < 1e7) return(strtrim(sprintf("%12.2f", x)))
    return(strtrim(sprintf("%12.3e", x)))
}

// a probability or cutoff for text: 0.05, not .05
string scalar _ivm_p(real scalar x)
{
    string scalar s
    s = strtrim(strofreal(x, "%10.0g"))
    if (substr(s, 1, 1) == ".") s = "0" + s
    if (substr(s, 1, 2) == "-.") s = "-0" + substr(s, 2, .)
    return(s)
}

string scalar _ivm_pv(real scalar p)
{
    if (p >= .) return("n/a")
    return(strtrim(sprintf("%10.4f", p)))
}

string scalar _ivm_stars(real scalar est, real scalar se, real rowvector cut)
{
    real scalar p
    if (est >= . | se >= . | se <= 0) return("")
    p = 2 * normal(-abs(est / se))
    if (p < cut[1]) return("***")
    if (p < cut[2]) return("**")
    if (p < cut[3]) return("*")
    return("")
}

void _ivm_row(string scalar lab, real scalar est, real scalar se, real scalar zc,
              real scalar dostars, real rowvector cut)
{
    string scalar st
    if (se < . & se > 0) {
        st = (dostars ? _ivm_stars(est, se, cut) : "")
        printf("{txt}%-24s{res}%s %s  %s %s {txt}%s\n", lab, _ivm_num(est, 10), _ivm_num(se, 10),
               _ivm_num(est - zc * se, 10), _ivm_num(est + zc * se, 10), st)
    }
    else {
        printf("{txt}%-24s{res}%s {txt}  SE unavailable\n", lab, _ivm_num(est, 10))
    }
}

string scalar _ivm_arstr(struct _ivm_ar scalar a)
{
    string scalar s, lo, hi
    real scalar   i, anyopen
    if (a.shape == "interval") return("[" + _ivm_numc(a.lower) + ", " + _ivm_numc(a.upper) + "]")
    if (a.shape == "rays") return("(-inf, " + _ivm_numc(a.left_upper) + "] U [" + _ivm_numc(a.right_lower) +
                                  ", +inf)   <- unbounded: weak-IV signature")
    if (a.shape == "rays_left_only")  return("(-inf, " + _ivm_numc(a.left_upper) + "]   <- unbounded: weak-IV signature")
    if (a.shape == "rays_right_only") return("[" + _ivm_numc(a.right_lower) + ", +inf)   <- unbounded: weak-IV signature")
    if (a.shape == "line") return("(-inf, +inf)   <- AR cannot reject any value at this level (weak-IV signature)")
    if (a.shape == "grid") {
        s = ""
        anyopen = 0
        for (i = 1; i <= rows(a.seg); i++) {
            lo = (a.seg[i, 3] ? "< " : "") + _ivm_numc(a.seg[i, 1])
            hi = (a.seg[i, 4] ? "> " : "") + _ivm_numc(a.seg[i, 2])
            if (a.seg[i, 3] | a.seg[i, 4]) anyopen = 1
            s = s + (i > 1 ? " U " : "") + "[" + lo + ", " + hi + "]"
        }
        if (anyopen) s = s + "   ('<' / '>': the segment reaches the search grid's boundary)"
        return(s)
    }
    if (a.shape == "unavailable") return("unavailable: negative variance estimate")
    return("empty set (AR rejects every value: model misspecification signal)")
}

void _ivm_selline(string scalar lab, real colvector s, real scalar lam, real scalar it,
                  real scalar conv, string rowvector nm)
{
    string scalar cv
    string rowvector sn
    real scalar   k
    k = sum(s)
    cv = (conv ? "converged in " : "did not converge after ") + strofreal(it) + " iteration" + (it != 1 ? "s" : "")
    printf("{txt}  S_%s (%g selected, lambda = %s, %s)\n", lab, k, _ivm_f2(lam), cv)
    if (k > 0) {
        sn = select(nm, s')
        if (cols(sn) > 20) printf("{res}    %s, ... (+%g more)\n", _ivm_join(sn[(1..20)], ", "), cols(sn) - 20)
        else               printf("{res}    %s\n", _ivm_join(sn, ", "))
    }
}

void _ivmediate_show(real scalar dostars, real scalar c1, real scalar c2, real scalar c3, real scalar showsel,
                  real scalar kappatable, real scalar r2sens)
{
    struct _ivm_state scalar G
    struct _ivm_res   scalar R
    string scalar    se, lvl, ab, s, ocrit
    string rowvector zn, absn, lab
    real rowvector   cut
    real scalar      K, zc, i, any_ub, P, j, ks, inf

    G = _ivm_get("_ivmediate_R")
    R = G.R
    K = R.K
    zc = R.zc
    zn = G.D.znames
    cut = (c1, c2, c3)
    se = _ivm_selab(G.O.se)
    lvl = strtrim(strofreal(100 * (1 - G.O.alpha), "%9.0g")) + "%"
    absn = tokens(st_local("absorb"))

    printf("\n{txt}{hline 78}\n")
    printf("{txt}Linear IV Mediation Analysis{col 58}ivmediate %s\n", _ivm_version())
    printf("{txt}{hline 78}\n")
    printf("{txt}Outcome:   {res}%-35s{txt}{col 48}Obs           = {res}%9.0fc\n", st_local("depvar"), R.n)
    printf("{txt}Treatment: {res}%-35s{txt}{col 48}Instruments   = {res}%9.0f\n", G.D.tname, K)
    printf("{txt}Mediator:  {res}%-35s{txt}{col 48}Covariates    = {res}%9.0f\n", G.D.mname, R.kcov)
    if (se == "cluster") {
        printf("{txt}SE:        {res}cluster (%s){txt}{col 48}Clusters      = {res}%9.0fc\n", st_local("clustervar"), R.G)
    }
    else if (se == "conley") {
        printf("{txt}SE:        {res}Conley (cutoff = %s, kernel = %s, metric = %s)\n",
               strtrim(strofreal(G.O.cutoff, "%12.0g")), G.O.kernel, G.O.metric)
    }
    else if (se == "robust") printf("{txt}SE:        {res}robust (HC1)\n")
    else                     printf("{txt}SE:        {res}unadjusted\n")
    if (cols(absn)) {
        printf("{txt}Absorb:    {res}%-35s{txt}{col 48}FE groups     = {res}%9.0fc\n", _ivm_join(absn, " "), R.nfeg)
        s = "df absorbed = " + strofreal(R.dffe) + "; singletons dropped = " + strofreal(R.nsing)
        if (cols(R.nested)) {
            if (any(R.nested)) {
                s = s + "; nested in the clusters: " + _ivm_join(select(absn, R.nested), ", ") + " (not counted in the CR1 df)"
            }
        }
        printf("{txt}           %s\n", s)
    }
    if (G.O.weighted) printf("{txt}Weights:   {res}%s\n", st_local("wtype"))
    if (R.ndum > 0) printf("{txt}Dummy indicators: {res}%g\n", R.ndum)
    if (cols(R.dropped)) printf("{txt}Dropped (collinear): {res}%s\n", _ivm_join(R.dropped, ", "))
    if (G.O.lasso) {
        printf("{txt}LASSO:     {res}post-double selection, plugin penalty, selected %g of %g candidates\n",
               R.lasso_nsel, R.lasso_p)
    }

    // ---- mediation estimates ----
    printf("\n{txt}Mediation estimates\n")
    printf("{txt}{hline 78}\n")
    printf("{txt}%-24s%10s %10s  %21s\n", "Parameter", "Estimate", "Std. err.", "[" + lvl + " conf. interval]")
    printf("{txt}{hline 78}\n")
    _ivm_row("tau_total (Wald)",       R.total,    R.se_total,    zc, dostars, cut)
    _ivm_row("tau (direct T->Y)",      R.tau,      R.se_tau,      zc, dostars, cut)
    _ivm_row("theta (M->Y | T)",       R.theta,    R.se_theta,    zc, dostars, cut)
    _ivm_row("gamma (T->M)",           R.gamma,    R.se_gamma,    zc, dostars, cut)
    _ivm_row("theta*gamma (indirect)", R.indirect, R.se_indirect, zc, dostars, cut)
    printf("{txt}{hline 78}\n")
    if (dostars) {
        printf("{txt}*** p<%s, ** p<%s, * p<%s (two-sided Wald, normal critical values)\n",
               _ivm_p(c1), _ivm_p(c2), _ivm_p(c3))
    }
    if (K == 1) {
        if (R.decomp_rel > 1e-8 & R.decomp_rel < .) {
            printf("{err}  WARNING: decomposition identity violated: tau + theta*gamma = %s vs tau_total = %s\n",
                   _ivm_numc(R.decomp_check), _ivm_numc(R.total))
            printf("{err}  (relative error %9.2e > 1e-08). For one instrument the identity is exact; please report.\n", R.decomp_rel)
        }
    }
    else if (R.decomp_rel > 1e-2 & R.decomp_rel < .) {
        printf("{txt}  Note: with K>1 instruments, tau + theta*gamma = %s and the Wald tau_total = %s differ\n",
               _ivm_numc(R.decomp_check), _ivm_numc(R.total))
        printf("{txt}  by %8.1e (relative). Expected: the equations are estimated one by one, not by joint GMM.\n", R.decomp_rel)
    }
    if (R.total != 0 & R.total < .) {
        printf("{txt}  Mediator explains {res}%s{txt}%% of the total effect.\n", strtrim(sprintf("%12.2f", 100 * R.indirect / R.total)))
    }

    // ---- first stage ----
    printf("\n{txt}First-stage diagnostics  (SE: %s)\n", se)
    printf("{txt}{hline 78}\n")
    s = ""
    for (i = 1; i <= K; i++) {
        s = s + (i > 1 ? ", " : "") + _ivm_numc(R.pi1[i]) + " (" + _ivm_numc(R.se_pi1[i]) + ")"
    }
    printf("{txt}  pi1 (T on Z | X)         = {res}%s\n", s)
    ocrit = (R.F_TZ_crit < . ? "   OP critical value = " + _ivm_f2(R.F_TZ_crit) + " (K_eff = " +
             strtrim(sprintf("%9.2f", R.K_eff_TZ)) + ", tau = " + _ivm_p(G.O.optau) + ", alpha = " +
             _ivm_p(G.O.opalpha) + ")" : "   OP critical value: unavailable")
    if (K > 1) {
        printf("{txt}  F(T, Z|X)                = {res}%10s{txt}  (joint Wald)\n", _ivm_f2(R.F_TZ))
        printf("{txt}  F(T, Z|X), OP effective  = {res}%10s{txt}%s\n", _ivm_f2(R.F_TZ_op), ocrit)
    }
    else printf("{txt}  F(T, Z|X)                = {res}%10s{txt}%s\n", _ivm_f2(R.F_TZ), ocrit)
    ocrit = (R.F_MZ_T_crit < . ? "   OP critical value = " + _ivm_f2(R.F_MZ_T_crit) + " (K_eff = " +
             strtrim(sprintf("%9.2f", R.K_eff_MZ_T)) + ", tau = " + _ivm_p(G.O.optau) + ", alpha = " +
             _ivm_p(G.O.opalpha) + ")" : "   OP critical value: unavailable")
    if (K > 1) {
        printf("{txt}  F(M, Z|T,X)              = {res}%10s{txt}  (joint Wald)\n", _ivm_f2(R.F_MZ_T))
        printf("{txt}  F(M, Z|T,X), OP effective= {res}%10s{txt}%s\n", _ivm_f2(R.F_MZ_T_op), ocrit)
        s = ""
        for (i = 1; i <= K; i++) s = s + "  " + strtrim(_ivm_f2(R.F_ind[i]))
        printf("{txt}    per-instrument t^2:   {res}%s\n", s)
    }
    else printf("{txt}  F(M, Z|T,X)              = {res}%10s{txt}%s\n", _ivm_f2(R.F_MZ_T), ocrit)

    // ---- Anderson-Rubin ----
    if (R.ar_done) {
        printf("\n{txt}Anderson-Rubin %s confidence sets  (weak-IV-robust, SE: %s)\n", lvl, se)
        printf("{txt}{hline 78}\n")
        any_ub = 0
        for (i = 1; i <= 3; i++) {
            if (i == 1) {
                printf("{txt}  %-18s{res}%s\n", "theta (M->Y | T)", _ivm_arstr(R.ar_theta))
                printf("{txt}  %-18s  AR p-value at 0 = %s   (%s)\n", "", _ivm_pv(R.ar_theta.pval0),
                       (R.ar_theta.method == "fieller" ? "exact Fieller inversion" : "exact statistic, grid search"))
                any_ub = any_ub | _ivm_ar_unbdisc(R.ar_theta)
            }
            else if (i == 2) {
                printf("{txt}  %-18s{res}%s\n", "gamma (T->M)", _ivm_arstr(R.ar_gamma))
                printf("{txt}  %-18s  AR p-value at 0 = %s   (%s)\n", "", _ivm_pv(R.ar_gamma.pval0),
                       (R.ar_gamma.method == "fieller" ? "exact Fieller inversion" : "exact statistic, grid search"))
                any_ub = any_ub | _ivm_ar_unbdisc(R.ar_gamma)
            }
            else {
                printf("{txt}  %-18s{res}%s\n", "tau_total", _ivm_arstr(R.ar_total))
                printf("{txt}  %-18s  AR p-value at 0 = %s   (%s)\n", "", _ivm_pv(R.ar_total.pval0),
                       (R.ar_total.method == "fieller" ? "exact Fieller inversion" : "exact statistic, grid search"))
                any_ub = any_ub | _ivm_ar_unbdisc(R.ar_total)
            }
        }
        if (any_ub) {
            printf("{txt}  Note: an unbounded or disconnected AR set is the expected outcome under weak\n")
            printf("{txt}  identification: the data cannot bound the parameter at this confidence level.\n")
        }
        printf("{txt}  AR sets (chi2(%g) critical value %s) remain valid whatever the first-stage strength;\n",
               R.ar_theta.df, strtrim(sprintf("%9.3f", R.ar_theta.cv)))
        printf("{txt}  the Wald CIs above require the effective F to clear the OP critical value.\n")
    }

    // ---- kappa ----
    printf("\n{txt}kappa-based sensitivity (MCA diagnostics)\n")
    printf("{txt}{hline 78}\n")
    if (G.O.weighted) printf("{txt}  (Hajek-weighted moments, %s)\n", st_local("wtype"))
    printf("{txt}  c_TM = Cov(eps_T, eps_M)   = {res}%s\n", _ivm_num(R.c_TM, 12))
    printf("{txt}  kappa* (breakeven)         = {res}%s\n", _ivm_num(R.kappa_star, 12))
    printf("{txt}  kappa*_CI (inf. breakeven) = {res}%s{txt}%s\n", _ivm_num(R.kappa_star_CI, 12),
           (R.ks_inf == 1 ? "" : "  (uninformative: the CI already contains 0)"))
    printf("{txt}  rho*_TY (as correlation)   = {res}%s{txt}%s\n", _ivm_num(R.rho_TY_star, 12),
           ((R.rho_TY_star < . & abs(R.rho_TY_star) > 1) ? "  (|rho*| > 1: no feasible error correlation overturns theta)" : ""))
    if (kappatable) {
        printf("\n{txt}  kappa-adjusted theta (kappa0 swept as a fraction of kappa*):\n")
        printf("{txt}  %8s %12s %12s %12s %12s\n", "% kappa*", "kappa0", "theta_adj", "CI_lo", "CI_hi")
        for (i = 1; i <= rows(R.kadj); i++) {
            printf("{res}  %7.0f%% %s %s %s %s\n", 100 * R.kadj[i, 2], _ivm_num(R.kadj[i, 1] + 0, 12),
                   _ivm_num(R.kadj[i, 3], 12), _ivm_num(R.kadj[i, 4], 12), _ivm_num(R.kadj[i, 5], 12))
        }
        // the conditional first stage below its OP critical value: one note under the table
        if (R.kt_rel == 0) {
            s = "Note: " + (K > 1 ? "F(M, Z|T,X), OP effective" : "F(M, Z|T,X)") + " = " + _ivm_f2(R.F_MZ_T_op) +
                " is below the OP critical value (" + _ivm_f2(R.F_MZ_T_crit) + ", tau = " + _ivm_p(G.O.optau) +
                "), so c_TM is estimated imprecisely" +
                (R.ar_done ? "; the AR set for theta above is robust to weak instruments." : ".")
            printf("{p 2 2 2}{txt}%s{p_end}\n", s)
        }
    }

    // ---- R^2 (with r2sensitivity) ----
    if (r2sens) {
        printf("\n{txt}R^2-based sensitivity (reduced-form robustness)\n")
        printf("{txt}{hline 78}\n")
        printf("{txt}  (t-statistics under %s SEs%s)\n", se, (K > 1 ? "; with K > 1, RV/XRV use the first instrument, " + zn[1] : ""))
        printf("{txt}  %-23s = {res}%10s\n", "F (reduced form, = t^2)", _ivm_f2(R.F_rf))
        printf("{txt}  %-23s = {res}%s\n", "RV(theta)", _ivm_num(R.RV, 10))
        printf("{txt}  %-23s = {res}%s\n", "XRV(theta)", _ivm_num(R.XRV, 10))
        printf("{txt}    (RV_rf = %s, RV_fs = %s)\n", _ivm_numc(R.rv_rf), _ivm_numc(R.rv_fs))
    }

    // ---- LASSO ----
    if (G.O.lasso) {
        printf("\n{txt}Post-double-selection LASSO (plugin penalty)\n")
        printf("{txt}{hline 78}\n")
        printf("{txt}  selected {res}%g{txt} of {res}%g{txt} expansion terms%s\n", R.lasso_nsel, R.lasso_p,
               (R.lasso_nkeep > 0 ? ", " + strofreal(R.lasso_nkeep) + " forced via lkeep()" : ""))
        if (R.lasso_nsel > 0) printf("{txt}  terms: {res}%s\n", _ivm_join(select(R.rich_names, R.lasso_union), " "))
        else                  printf("{txt}  terms: {res}(none)\n")
        if (showsel & R.lasso_p > 0) {
            lab = ("Y", "T", "M")
            for (j = 1; j <= K; j++) lab = (lab, "Z" + strofreal(j))
            for (j = 1; j <= 3 + K; j++) {
                _ivm_selline(lab[j], R.lasso_S[., j], R.lasso_lambda[j], R.lasso_niter[j],
                             R.lasso_conv[j], R.rich_names)
            }
            printf("{txt}  Union S = S_Y + S_T + S_M + S_Z: %g of %g candidates\n", R.lasso_nsel, R.lasso_p)
        }
    }

    // ---- bootstrap ----
    if (G.B > 0) {
        P = cols(G.pnames)
        printf("\n{txt}Bootstrap  (method = %s, B = %g, %s)\n", G.bmethod, G.B,
               (G.has_idx ? "replicates from bootidx()" : "seed = " + st_local("bseed")))
        printf("{txt}{hline 78}\n")
        printf("{txt}  CI type: %s (%s);  %g successful replicates%s\n", G.bcitype, lvl, G.B - G.nfailed,
               (G.nfailed > 0 ? ", " + strofreal(G.nfailed) + " failed" : ""))
        printf("{txt}  %-24s%10s %10s  %21s\n", "Parameter", "Estimate", "Boot SE", "[Boot " + lvl + " CI]")
        lab = ("tau_total (Wald)", "tau (direct T->Y)", "theta (M->Y | T)", "gamma (T->M)",
               "theta*gamma (indirect)", "kappa*", "rho*_TY", "RV(theta)", "XRV(theta)")
        for (j = 1; j <= P; j++) {
            printf("{txt}  %-24s{res}%s %s  %s %s {txt}%s\n", lab[j], _ivm_num(G.bpoint[j], 10),
                   _ivm_num(G.bse[j], 10), _ivm_num(G.bci[j, 1], 10), _ivm_num(G.bci[j, 2], 10),
                   ((dostars & j <= 5) ? _ivm_stars(G.bpoint[j], G.bse[j], cut) : ""))
        }
        ks = (R.theta - sign(R.theta) * zc * G.bse[3]) * R.c_TM
        inf = (G.bse[3] < . ? abs(R.theta) > zc * G.bse[3] : 0)
        printf("{txt}  kappa*_CI with the bootstrap SE of theta = {res}%s{txt}%s\n", _ivm_numc(ks),
               (inf ? "" : "  (uninformative: the CI already contains 0)"))
        printf("{txt}  The mediation estimates table, e(b) and e(V) keep the analytic SEs.\n")
    }

    printf("{txt}{hline 78}\n")
    printf("{txt}Excluded instrument(s): {res}%s\n", _ivm_join(zn, " "))
}

end
