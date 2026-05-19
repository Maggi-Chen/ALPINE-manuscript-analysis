"""Parse GNU time -v reports into a CSV summary for the 4-quadrant simulation benchmark."""
from __future__ import annotations

import csv
import re
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
LOGS = BENCH / "timing_logs"

SAMPLES = [
    ("hifi_trac", "PacBio HiFi TRAC", 60000, "pacbio"),
    ("hifi_trbc", "PacBio HiFi TRBC", 60000, "pacbio"),
    ("ont_trac",  "Nanopore TRAC",    60000, "nanopore"),
    ("ont_trbc",  "Nanopore TRBC",    60000, "nanopore"),
]


def parse_elapsed(s: str) -> float:
    parts = s.split(":")
    if len(parts) == 3:
        h, m, sec = parts
        return int(h) * 3600 + int(m) * 60 + float(sec)
    if len(parts) == 2:
        m, sec = parts
        return int(m) * 60 + float(sec)
    return float(s)


def parse_report(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if ":" not in line:
            continue
        key, _, val = line.strip().partition(":")
        out[key.strip()] = val.strip()
    return out


def main() -> None:
    rows = []
    for label, alias, reads, platform in SAMPLES:
        report = LOGS / f"time_{label}.txt"
        if not report.exists() or report.stat().st_size == 0:
            print(f"MISSING: {report}")
            continue
        r = parse_report(report)
        m = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", report.read_text())
        elapsed_s = parse_elapsed(m.group(1)) if m else float("nan")
        rss_kb = int(r["Maximum resident set size (kbytes)"])
        rows.append({
            "label": label,
            "alias": alias,
            "platform": platform,
            "reads_input": reads,
            "elapsed_seconds": round(elapsed_s, 2),
            "elapsed_minutes": round(elapsed_s / 60, 2),
            "max_rss_kb": rss_kb,
            "max_rss_gb": round(rss_kb / 1024 / 1024, 2),
            "user_seconds": float(r["User time (seconds)"]),
            "sys_seconds": float(r["System time (seconds)"]),
            "cpu_percent": r["Percent of CPU this job got"].rstrip("%"),
            "reads_per_second": round(reads / elapsed_s, 1) if elapsed_s > 0 else None,
        })

    if not rows:
        print("No timing reports found.")
        return

    out_path = BENCH / "benchmark_summary.csv"
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_path}")
    print()
    print(f"{'Quadrant':22} {'reads':>6}  {'time':>8}  {'RAM':>6}  {'CPU%':>5}  {'reads/s':>8}")
    for r in rows:
        print(
            f"{r['alias']:22} "
            f"{r['reads_input']:>6d}  "
            f"{r['elapsed_minutes']:>6.2f}min  "
            f"{r['max_rss_gb']:>4.2f}GB  "
            f"{r['cpu_percent']:>4s}%  "
            f"{r['reads_per_second']:>6}/s"
        )


if __name__ == "__main__":
    main()
