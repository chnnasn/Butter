"""Summarize measured configurations, retaining failed quality gates.

Usage: python benchmarks/summarize_comparison.py build/comparison
No speedup ratio is inferred between unequal physics trajectories.
"""
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

root = Path(sys.argv[1])
groups = defaultdict(list)
for path in sorted(root.glob("*.csv")):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            key = tuple(row[k] for k in ("engine", "mode", "scenario", "iterations", "substeps", "activity"))
            groups[key].append(row)

writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(["engine", "mode", "scenario", "iterations", "substeps", "activity",
                 "median_mean_ms", "median_p95_ms", "worst_penetration_m", "worst_joint_gap_m",
                 "worst_top_error_m", "all_quality_gates", "median_cpp_new_calls", "samples"])
for key, rows in sorted(groups.items()):
    median = lambda field: statistics.median(float(r[field]) for r in rows)
    worst = lambda field: max(float(r[field]) for r in rows)
    writer.writerow([*key, median("mean_ms"), median("p95_ms"), worst("max_penetration_m"),
                     worst("max_joint_gap_m"), worst("final_top_error_m"),
                     int(all(r["quality_pass"] == "1" for r in rows)), median("cpp_new_calls"), len(rows)])
