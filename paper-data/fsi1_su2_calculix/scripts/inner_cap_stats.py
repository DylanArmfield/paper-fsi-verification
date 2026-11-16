#!/usr/bin/env python3
"""Inner-iteration statistics of a coupled SU2 run (history_0*.csv, one row per SU2 solve, i.e. per
coupling iteration): solves that hit the inner cap (Inner_Iter + 1 >= CAP), the last time step with a
capped solve, and the share of capped solves in the closing 20 % of the run.
Usage: inner_cap_stats.py CASEDIR [CAP=200]"""
import csv, glob, sys
c = sys.argv[1]; cap = int(sys.argv[2]) if len(sys.argv) > 2 else 200
f = sorted(glob.glob(f"{c}/fluid-su2/history_0*.csv"))[0]
R = [{k.strip().strip('"'): v.strip() for k, v in r.items()} for r in csv.DictReader(open(f))]
t = [float(r["Cur_Time"]) + 2.0 for r in R]  # SU2 Cur_Time counts from the coupled start (t = 2 s)
it = [int(r["Inner_Iter"]) + 1 for r in R]
capped = [ti for ti, n in zip(t, it) if n >= cap]
res = [max(float(r["rms[P]"]), float(r["rms[U]"]), float(r["rms[V]"])) for r in R]
cres = [x for x, n in zip(res, it) if n >= cap]
t0, t1 = min(t), max(t); tc = 30.0 - 0.2 * 28.0  # closing 20 % of t = 2-30 s, as evaluate_fsi1.py
late = [n for ti, n in zip(t, it) if ti >= tc]
print(f"capped solves: max rms at stop {max(cres) if cres else 0:.2f}, median {sorted(cres)[len(cres)//2] if cres else 0:.2f} (log10)")
print(f"{c}: solves {len(R)}, capped {len(capped)} ({100*len(capped)/len(R):.1f} %), last capped at t = {max(capped) if capped else '-'} s "
      f"(run {t0:g}-{t1:g} s); closing 20 % (t >= {tc:.2f}): {len(late)} solves, capped {sum(n >= cap for n in late)}, "
      f"inner its mean {sum(late)/len(late):.1f} max {max(late)}")
