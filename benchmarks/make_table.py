"""Speed-comparison table (benchmarks/README.md).

Reads results/{python,r,stata,matlab}.csv, checks that the four languages give
the same point estimates (relative difference below 1e-8) for every size and
specification, and writes results/table.md: an HTML table with the median run
time in seconds (rows N; the defaults, then the defaults plus the
post-double-selection LASSO) and a note with the machine, the program versions
and the point estimates.

Usage:  python make_table.py
"""
import os
import platform
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
LANGS = [("python", "Python"), ("r", "R"), ("stata", "Stata"), ("matlab", "MATLAB")]
SPECS = [("default", "Defaults"), ("lasso", "Defaults + post-double-selection LASSO")]
EST = ["theta", "gamma", "tau", "total"]
TOL = 1e-8


def machine():
    """CPU, cores and memory (Windows via CIM; elsewhere what platform knows)."""
    if platform.system() == "Windows":
        ps = ("$c = Get-CimInstance Win32_Processor | Select-Object -First 1; "
              "$m = (Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory; "
              "$o = (Get-CimInstance Win32_OperatingSystem).Caption; "
              "\"$($c.Name.Trim())|$($c.NumberOfCores)|$($c.NumberOfLogicalProcessors)|$m|$o\"")
        try:
            out = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True,
                                 text=True, timeout=60).stdout.strip()
            name, cores, logical, mem, os_name = out.split("|")
            return (f"{name} ({cores} cores, {logical} threads), "
                    f"{round(int(mem) / 2**30)} GB RAM, {os_name}")
        except Exception:
            pass
    return f"{platform.processor()}, {platform.platform()}"


def main():
    frames = []
    for key, _ in LANGS:
        path = os.path.join(RES, f"{key}.csv")
        if not os.path.exists(path):
            sys.exit(f"missing {path}: run the {key} benchmark first")
        frames.append(pd.read_csv(path))
    df = pd.concat(frames, ignore_index=True)
    df["language"] = df["language"].astype(str)
    sizes = sorted(df["n"].unique())

    # identical point estimates across the four languages
    worst = 0.0
    for (n, spec), g in df.groupby(["n", "spec"]):
        if len(g) != len(LANGS):
            sys.exit(f"n={n} {spec}: {len(g)} languages instead of {len(LANGS)}")
        for c in EST:
            v = g[c].astype(float)
            ref = v.iloc[0]
            rel = ((v - ref).abs() / max(abs(ref), 1e-300)).max()
            worst = max(worst, rel)
            if rel > TOL:
                sys.exit(f"n={n} {spec} {c}: languages differ (max relative difference {rel:.2e})")

    def cell(n, spec, lang):
        v = df[(df.n == n) & (df.spec == spec) & (df.language == lang)]["median"].iloc[0]
        return f"{v:.2f}" if v < 100 else f"{v:.0f}"

    h = ["<table>", "<thead>",
         "<tr><th rowspan=\"2\">N</th>"
         + "".join(f"<th colspan=\"{len(LANGS)}\">{lab}</th>" for _, lab in SPECS) + "</tr>",
         "<tr>" + "".join(f"<th>{lab}</th>" for _ in SPECS for _, lab in LANGS) + "</tr>",
         "</thead>", "<tbody>"]
    for n in sizes:
        h.append(f"<tr><td align=\"right\">{n:,}</td>"
                 + "".join(f"<td align=\"right\">{cell(n, s, lab)}</td>"
                           for s, _ in SPECS for _, lab in LANGS) + "</tr>")
    h += ["</tbody>", "</table>"]

    env = []
    for key, lab in LANGS:
        p = os.path.join(RES, f"env_{key}.txt")
        env.append(open(p).read().strip() if os.path.exists(p) else f"{lab}: (no environment file)")

    est = []
    for spec, lab in SPECS:
        parts = []
        for n in sizes:
            g = df[(df.n == n) & (df.spec == spec)].iloc[0]
            parts.append(f"N = {n:,}: θ̂ = {g.theta:.4f}, γ̂ = {g.gamma:.4f}, θ̂γ̂ = {g.theta * g.gamma:.4f}, "
                         f"τ̂ = {g.tau:.4f}, τ̂_total = {g.total:.4f}")
        est.append(f"*{lab}:* " + "; ".join(parts) + ".")

    note = [
        "",
        "*Median wall-clock time in seconds of three runs of the estimation call, after one warm-up run; "
        "loading the data and printing are excluded. Simulated data (benchmarks/gen_data.py): one instrument, "
        "50 covariates and a fixed effect with 100 groups, absorbed. Defaults: robust standard errors, "
        "Anderson–Rubin sets and the κ diagnostics. The LASSO selects the controls from the 50 covariates and "
        "their squares and cubes (200 candidates, plug-in penalty).*",
        "",
        f"*Machine: {machine()}. Software: " + "; ".join(env) + ".*",
        "",
        f"*The four languages give the same point estimates (largest relative difference {worst:.0e}):*",
        "",
    ] + est
    with open(os.path.join(RES, "table.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(h + note) + "\n")
    print("\n".join(h + note))


if __name__ == "__main__":
    main()
