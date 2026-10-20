<table>
<thead>
<tr><th rowspan="2">N</th><th colspan="4">Defaults</th><th colspan="4">Defaults + post-double-selection LASSO</th></tr>
<tr><th>Python</th><th>R</th><th>Stata</th><th>MATLAB</th><th>Python</th><th>R</th><th>Stata</th><th>MATLAB</th></tr>
</thead>
<tbody>
<tr><td align="right">1,000</td><td align="right">0.08</td><td align="right">0.11</td><td align="right">0.07</td><td align="right">0.07</td><td align="right">0.08</td><td align="right">0.11</td><td align="right">0.11</td><td align="right">0.05</td></tr>
<tr><td align="right">10,000</td><td align="right">0.74</td><td align="right">1.12</td><td align="right">0.78</td><td align="right">0.32</td><td align="right">2.54</td><td align="right">6.33</td><td align="right">4.37</td><td align="right">1.19</td></tr>
<tr><td align="right">100,000</td><td align="right">5.36</td><td align="right">15.69</td><td align="right">10.94</td><td align="right">1.84</td><td align="right">20.55</td><td align="right">92.86</td><td align="right">127</td><td align="right">10.39</td></tr>
</tbody>
</table>

*Median wall-clock time in seconds of three runs of the estimation call, after one warm-up run; loading the data and printing are excluded. Simulated data (benchmarks/gen_data.py): one instrument, 50 covariates and a fixed effect with 100 groups, absorbed. Defaults: robust standard errors, Anderson–Rubin sets and the κ diagnostics. The LASSO selects the controls from the 50 covariates and their squares and cubes (200 candidates, plug-in penalty).*

*Machine: Intel(R) Core(TM) Ultra 9 285K (24 cores, 24 threads), 63 GB RAM, Microsoft Windows 11 Home. Software: Python 3.13.9; ivmediate 2.0.0; numpy 2.3.5 (mkl-sdl 2025); pandas 2.3.3; R version 4.5.2 (2025-10-31 ucrt); ivmediate 2.0.0; BLAS (R's reference BLAS); LAPACK 3.12.1; Stata 19.5 SE (processors used: 1); ivmediate 2.0.0; MATLAB 26.2.0.3386108 (R2026b); ivmediate 2.0.0; Intel(R) oneAPI Math Kernel Library Version 2025.3-Product Build 20251007 for Intel(R) 64 architecture applications (CNR branch AVX2); computational threads 24.*

*The four languages give the same point estimates (largest relative difference 1e-12):*

*Defaults:* N = 1,000: θ̂ = 0.8823, γ̂ = 0.4723, θ̂γ̂ = 0.4167, τ̂ = -0.0019, τ̂_total = 0.4148; N = 10,000: θ̂ = 0.5209, γ̂ = 0.3933, θ̂γ̂ = 0.2048, τ̂ = 0.2951, τ̂_total = 0.5000; N = 100,000: θ̂ = 0.4857, γ̂ = 0.4026, θ̂γ̂ = 0.1955, τ̂ = 0.3077, τ̂_total = 0.5032.
*Defaults + post-double-selection LASSO:* N = 1,000: θ̂ = 1.1975, γ̂ = 0.4299, θ̂γ̂ = 0.5148, τ̂ = -0.1451, τ̂_total = 0.3697; N = 10,000: θ̂ = 0.5214, γ̂ = 0.3930, θ̂γ̂ = 0.2049, τ̂ = 0.2946, τ̂_total = 0.4995; N = 100,000: θ̂ = 0.4856, γ̂ = 0.4025, θ̂γ̂ = 0.1955, τ̂ = 0.3077, τ̂_total = 0.5032.
