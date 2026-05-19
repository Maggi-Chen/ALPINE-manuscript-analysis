"""
Compute per-class benchmark metrics for knock-knock vs ALPINE on Maggie's v1.6
PBSIM3 simulation (HiFi + ONT × TRAC + TRBC, 60,000 reads per quadrant).

Truth labels are embedded in read names:
  simulated_<class>.<NNNNN>
where <class> is one of 12 PBSIM3-generated variants (e.g. trac_cd19_hdr,
trac_del_large, trbc_wt, etc).

knock-knock's per-read calls live in:
  run/results/<batch>/<sample>/outcome_list.txt
with (tab-separated) columns: read_name, UMI, length_ratio, quality,
category, subcategory, details...

This script:
  1. Parses the 4 knock-knock outcome_list.txt files.
  2. Maps each (truth_class, knock-knock category/subcategory) pair to an
     ALPINE-compatible class (the same 13 classes in alpine_v1_6_per_class_metrics.csv).
  3. Emits a knock-knock row set in ALPINE's table format and writes
     knock_knock_v0_8_2_per_class_metrics.csv.

Note on mapping: knock-knock was designed for amplicons where primers sit
OUTSIDE the homology arms. Maggie's simulation places primers AT the HA edges
(primer-inside-HA design), so knock-knock calls most HDR reads "complex
misintegration". We therefore include explicit mapping rows that faithfully
represent knock-knock's native labels rather than inventing a correspondence
that doesn't exist.
"""

import argparse
import collections
import pathlib
import sys


# Map (alpine-schema class, knock-knock category, knock-knock subcategory)
# Knock-knock's classifier is hierarchical; we collapse (category, subcategory)
# onto the closest ALPINE class. Entries later in the list do NOT override
# earlier ones — the match is first-wins using a priority order in the mapper.

# ALPINE's canonical class schema (13 classes + Unclassified).
ALPINE_CLASSES = [
    "Unmodified",
    "Unmodified-with-SNP",
    "DEL-small",
    "DEL-large",
    "INS-small",
    "INS-large",
    "HDR-CD19",
    "HDR-CD22",
    "Non-HDR-with-ITR-CD19",
    "Non-HDR-with-ITR-CD22",
    "Non-HDR-without-ITR-CD19",
    "Non-HDR-without-ITR-CD22",
    "Unclassified",
]


# Truth-class token → ALPINE-class (matching Maggie's canonical v1.6 schema).
# Key subtlety: `truncated_hdr` reads are mapped to Non-HDR-without-ITR-<AAV>
# in Maggie's schema, not to HDR-<AAV>, because biologically a truncated
# integration lacks full HDR architecture and has no ITR — so ALPINE labels
# it as a "Non-HDR-without-ITR" insertion of the donor payload.
TRUTH_TO_ALPINE = {
    "wt":                 "Unmodified",
    "snp":                "Unmodified-with-SNP",
    "del_small":          "DEL-small",
    "del_large":          "DEL-large",
    "ins_small":          "INS-small",
    "ins_large":          "INS-large",
    "cd19_hdr":           "HDR-CD19",
    "cd19_truncated_hdr": "Non-HDR-without-ITR-CD19",
    "cd22_hdr":           "HDR-CD22",
    "cd22_truncated_hdr": "Non-HDR-without-ITR-CD22",
    "cd19_itr":           "Non-HDR-with-ITR-CD19",
    "cd19_itr_oneside":   "Non-HDR-with-ITR-CD19",
    "cd22_itr":           "Non-HDR-with-ITR-CD22",
    "cd22_itr_oneside":   "Non-HDR-with-ITR-CD22",
}


def _parse_insertion_seq(details):
    """Extract inserted sequence from knock-knock details field.
    Format: insertions=I:<pos>,<SEQ> (may have multiple; take first)."""
    if "insertions=" not in details:
        return None
    try:
        # insertions=I:309,ACGGTCT  -> SEQ = 'ACGGTCT'
        chunk = details.split("insertions=", 1)[1]
        first = chunk.split(";")[0]  # handle multi-insert case
        # first looks like 'I:309,ACGGTCT'
        if "," in first:
            return first.split(",", 1)[1].split()[0]
    except Exception:
        return None
    return None


def parse_truth_class(read_name):
    """Read name: simulated_trac_cd19_hdr.00001  ->  truth class 'HDR-CD19'.

    NanoporeExperiment prepends a 10-digit index (`NNNNNNNNNN_`) to each read
    name during preprocessing, so we also strip that prefix if present.
    """
    # Strip NanoporeExperiment's 10-digit prefix if present
    if len(read_name) > 11 and read_name[10] == "_" and read_name[:10].isdigit():
        read_name = read_name[11:]
    # Strip 'simulated_' and trailing '.NNNNN'
    stripped = read_name.replace("simulated_", "", 1).rsplit(".", 1)[0]
    # First token is site (trac/trbc), rest is class token
    parts = stripped.split("_", 1)
    if len(parts) != 2:
        return None
    site_token, class_token = parts
    return TRUTH_TO_ALPINE.get(class_token, "Unclassified")


def map_kk_to_alpine(category, subcategory, site, details):
    """
    Map a knock-knock (category, subcategory, details, site) call to an
    ALPINE-compatible class. 'site' is 'TRAC' or 'TRBC'. Unmapped combos
    fall through to 'Unclassified' — this is intentional and honest, not
    a failure: knock-knock's scope is narrower than ALPINE's on this data.
    """
    cat = category.lower()
    sub = (subcategory or "").lower()

    if cat == "wt":
        # knock-knock's WT bucket is closest to ALPINE 'Unmodified'
        return "Unmodified"

    if cat == "simple indel":
        # subcategory tells us insertion/deletion; for insertions, parse size
        # from details field (format: insertions=I:pos,SEQ).
        if "insertion" in sub:
            ins_seq = _parse_insertion_seq(details)
            if ins_seq is not None and len(ins_seq) >= 50:
                return "INS-large"
            return "INS-small"
        if "deletion <50" in sub:
            return "DEL-small"
        if "deletion >=50" in sub:
            return "DEL-large"
        return "Unclassified"

    if cat == "hdr":
        # For genuine HDR calls — rare on this simulation due to amplicon-design
        # mismatch. Map by inferred site.
        return f"HDR-CD19" if site == "TRAC" else f"HDR-CD22"

    # Everything else (complex misintegration, blunt misintegration,
    # uncategorized, malformed layout, etc) lands in Unclassified. This
    # captures the scope limitation faithfully.
    return "Unclassified"


def site_from_batch(batch_name):
    b = batch_name.lower()
    if "trac" in b:
        return "TRAC"
    if "trbc" in b:
        return "TRBC"
    raise ValueError(f"cannot infer site from batch {batch_name!r}")


def platform_from_batch(batch_name):
    b = batch_name.lower()
    if b.startswith("hifi"):
        return "HiFi"
    if b.startswith("ont"):
        return "ONT"
    raise ValueError(f"cannot infer platform from batch {batch_name!r}")


def parse_outcome_list(path):
    """Yield (read_name, category, subcategory, details) per read."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            # Expected: name, UMI, len_ratio, qual, category, subcategory, details...
            if len(parts) < 6:
                continue
            name = parts[0]
            category = parts[4]
            subcategory = parts[5]
            details = "\t".join(parts[6:]) if len(parts) > 6 else ""
            yield name, category, subcategory, details


def compute_per_class_metrics(truth_labels, predicted_labels, classes):
    """
    For each class in `classes`, compute TP/FP/FN/precision/recall/F1.
    truth_labels and predicted_labels are parallel lists of ALPINE-schema
    class names.
    """
    rows = []
    for cls in classes:
        tp = sum(1 for t, p in zip(truth_labels, predicted_labels) if t == cls and p == cls)
        fp = sum(1 for t, p in zip(truth_labels, predicted_labels) if t != cls and p == cls)
        fn = sum(1 for t, p in zip(truth_labels, predicted_labels) if t == cls and p != cls)
        support = sum(1 for t in truth_labels if t == cls)
        prec = 100.0 * tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec = 100.0 * tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
        rows.append({
            "class": cls,
            "TP": tp, "FP": fp, "FN": fn, "support": support,
            "precision": round(prec, 2),
            "recall": round(rec, 2),
            "f1": round(f1, 2),
        })
    return rows


def process_sample(batch_name, sample_name, base_dir, tool_name):
    outcome_path = base_dir / "results" / batch_name / sample_name / "outcome_list.txt"
    if not outcome_path.exists():
        print(f"  SKIPPED: {outcome_path} not found", file=sys.stderr)
        return []

    site = site_from_batch(batch_name)
    platform = platform_from_batch(batch_name)

    truths, preds = [], []
    for name, cat, sub, details in parse_outcome_list(outcome_path):
        truth = parse_truth_class(name)
        if truth is None:
            continue
        pred = map_kk_to_alpine(cat, sub, site, details)
        truths.append(truth)
        preds.append(pred)

    per_class = compute_per_class_metrics(truths, preds, ALPINE_CLASSES)

    # Add a 'Total' row = overall accuracy (reads correct / reads total)
    total_reads = len(truths)
    correct = sum(1 for t, p in zip(truths, preds) if t == p)
    total_acc = round(100.0 * correct / total_reads, 2) if total_reads else 0.0
    per_class.append({
        "class": "Total",
        "TP": correct,
        "FP": total_reads - correct,
        "FN": total_reads - correct,
        "support": total_reads,
        "precision": total_acc,
        "recall": total_acc,
        "f1": total_acc,
    })

    for row in per_class:
        row["tool"] = tool_name
        row["platform"] = platform
        row["site"] = site
    return per_class


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        default="/home/gaox36/bioinfo_tools/on_target_LRS_pipeline/alpine/manuscript_submission/revision/external_datasets/analysis/knock-knock/simulations_v1.6/run",
        help="knock-knock run directory (contains results/)",
    )
    parser.add_argument(
        "--out",
        default="knock_knock_v0_8_2_per_class_metrics.csv",
        help="output CSV file (ALPINE-compatible schema)",
    )
    parser.add_argument(
        "--alpine-csv",
        default="alpine_v1_6_per_class_metrics.csv",
        help="existing ALPINE per-class metrics CSV to merge with (optional)",
    )
    parser.add_argument(
        "--merged-out",
        default="alpine_vs_knock_knock_per_class_metrics.csv",
        help="output CSV file for merged ALPINE + knock-knock metrics",
    )
    parser.add_argument("--tool-name", default="knock-knock")
    args = parser.parse_args()

    base_dir = pathlib.Path(args.base_dir)

    samples = [
        ("hifi_trac", "simulation_hifi_trac"),
        ("hifi_trbc", "simulation_hifi_trbc"),
        ("ont_trac", "simulation_ont_trac"),
        ("ont_trbc", "simulation_ont_trbc"),
    ]

    all_rows = []
    for batch, sample in samples:
        print(f"Processing {batch}/{sample}...", file=sys.stderr)
        rows = process_sample(batch, sample, base_dir, args.tool_name)
        all_rows.extend(rows)

    # Write knock-knock-only metrics CSV
    header = ["tool", "platform", "site", "class", "TP", "FP", "FN", "support", "precision", "recall", "f1"]
    out_path = pathlib.Path(args.out)
    with out_path.open("w") as fh:
        fh.write(",".join(header) + "\n")
        for row in all_rows:
            fh.write(",".join(str(row[col]) for col in header) + "\n")
    print(f"\nWrote {len(all_rows)} knock-knock rows to {out_path}", file=sys.stderr)

    # Optionally merge with the ALPINE CSV
    alpine_path = pathlib.Path(args.alpine_csv)
    if alpine_path.exists():
        merged_path = pathlib.Path(args.merged_out)
        with alpine_path.open() as src, merged_path.open("w") as dst:
            alpine_content = src.read()
            dst.write(alpine_content)
            if not alpine_content.endswith("\n"):
                dst.write("\n")
            for row in all_rows:
                dst.write(",".join(str(row[col]) for col in header) + "\n")
        print(f"Wrote merged ALPINE + knock-knock table to {merged_path}", file=sys.stderr)
    else:
        print(f"NOTE: ALPINE CSV {alpine_path} not found; skipping merge.", file=sys.stderr)


if __name__ == "__main__":
    main()
