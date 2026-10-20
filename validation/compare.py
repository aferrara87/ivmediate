"""Compare out_py.json (the Python reference) with another language's output.

Usage: python compare.py [r | matlab | stata | <path to a JSON file>]   (default: r)
  r       out_r.json       (run_r.R)
  matlab  out_matlab.json  (the MATLAB runner)
  stata   out_stata.json   (stata/run_stata.do)

Numbers agree within rel 1e-8 or abs 1e-10; text must match exactly. A config
that stops with an error must stop with an error in both languages, with the
same "_error_class" (the messages may differ). Keys starting with "_" are
metadata and are not compared key by key. Exit code 0 iff every config agrees.
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = {"r": "out_r.json", "matlab": "out_matlab.json", "stata": "out_stata.json"}
arg = sys.argv[1] if len(sys.argv) > 1 else "r"
if arg.lower() in TARGETS:
    LANG, other_path = arg.upper(), os.path.join(HERE, TARGETS[arg.lower()])
else:
    LANG, other_path = os.path.basename(arg), arg
with open(os.path.join(HERE, "out_py.json"), encoding="utf-8") as f:
    PY = json.load(f)
with open(other_path, encoding="utf-8") as f:
    OT = json.load(f)

REL = 1e-8
ABS = 1e-10
NUM_RE = re.compile(r"^[+-]?((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+)?$")


def is_num(s):
    return isinstance(s, str) and bool(NUM_RE.match(s))


def close(a, b):
    if a == b:
        return True
    if a in ("NaN",) and b in ("NaN",):
        return True
    if is_num(a) and is_num(b):
        fa, fb = float(a), float(b)
        if math.isnan(fa) and math.isnan(fb):
            return True
        return abs(fa - fb) <= max(ABS, REL * max(abs(fa), abs(fb)))
    return False


def reldiff(a, b):
    try:
        fa, fb = float(a), float(b)
        d = abs(fa - fb)
        m = max(abs(fa), abs(fb))
        return f"absdiff={d:.3e} reldiff={(d / m if m else 0):.3e}"
    except Exception:
        return ""


def status(out):
    return "error" if ("_ERROR" in out or out.get("_status") == "error") else "ok"


def describe(out):
    if status(out) == "ok":
        return "ok"
    return f"error [{out.get('_error_class', 'no class')}] {out.get('_ERROR', '')}"


n_cfg_ok = 0
n_cfg_bad = 0
report = []
for cfg in sorted(set(PY) | set(OT)):
    if cfg not in PY:
        report.append(f"### {cfg}: MISSING IN PYTHON OUTPUT"); n_cfg_bad += 1; continue
    if cfg not in OT:
        report.append(f"### {cfg}: MISSING IN {LANG} OUTPUT"); n_cfg_bad += 1; continue
    p, o = PY[cfg], OT[cfg]
    if status(p) == "error" or status(o) == "error":
        cls = p.get("_error_class")
        if status(p) == status(o) == "error" and cls is not None and cls == o.get("_error_class"):
            n_cfg_ok += 1
            report.append(f"### {cfg}: OK (both stop with an error of class {cls})")
        else:
            n_cfg_bad += 1
            report.append(f"### {cfg}: STATUS DIFFERS: PY {describe(p)} | {LANG} {describe(o)}")
        continue
    diffs = []
    for k in sorted(set(p) | set(o)):
        if k.startswith("_"):
            continue
        pv, ov = p.get(k, "<missing>"), o.get(k, "<missing>")
        if not close(pv, ov):
            extra = reldiff(pv, ov) if (is_num(pv) and is_num(ov)) else ""
            diffs.append(f"  {k}: py={pv!r}  {LANG.lower()}={ov!r}  {extra}")
    if diffs:
        n_cfg_bad += 1
        report.append(f"### {cfg}: {len(diffs)} mismatched keys")
        report.extend(diffs[:80])
        if len(diffs) > 80:
            report.append(f"  ... (+{len(diffs) - 80} more)")
    else:
        n_cfg_ok += 1
        report.append(f"### {cfg}: OK ({len([k for k in p if not k.startswith('_')])} keys agree)")

print("\n".join(report))
print(f"\n==== SUMMARY (Python vs {LANG}): {n_cfg_ok} configs fully agree, "
      f"{n_cfg_bad} configs with differences ====")
sys.exit(0 if n_cfg_bad == 0 else 1)
