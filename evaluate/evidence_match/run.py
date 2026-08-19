from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from evaluate.common import add_project_paths, read_tsv, write_tsv
from evaluate.evidence_match.config import (
    GENIE_NORMALIZED_PATH,
    HIERARCHY_PATH,
    INTERMEDIATE_DIR,
    MSK_NORMALIZED_PATH,
    RESULT_DIR,
    SOURCES,
    VARIANTS_PATH,
)
from evaluate.evidence_match.patient_genie import write_normalized_records as write_genie_normalized_records
from evaluate.evidence_match.patient_msk import write_normalized_records as write_msk_normalized_records


add_project_paths()
from hierarchy_class import VariationHierarchicalRelation, load_hierarchy  # noqa: E402


SMOKE_INTERMEDIATE_DIR = INTERMEDIATE_DIR / "smoke"
SMOKE_GENIE_NORMALIZED_PATH = SMOKE_INTERMEDIATE_DIR / GENIE_NORMALIZED_PATH.name
SMOKE_MSK_NORMALIZED_PATH = SMOKE_INTERMEDIATE_DIR / MSK_NORMALIZED_PATH.name


def load_patient_records(cohort: str, sample_limit: int | None = None) -> list[dict[str, str]]:
    paths = []
    if cohort in ("genie", "all"):
        path = SMOKE_GENIE_NORMALIZED_PATH if sample_limit else GENIE_NORMALIZED_PATH
        if not path.exists():
            print(f"Building GENIE normalized intermediate: {path}", flush=True)
            write_genie_normalized_records(path, sample_limit=sample_limit)
        paths.append(path)
    if cohort in ("msk", "all"):
        path = SMOKE_MSK_NORMALIZED_PATH if sample_limit else MSK_NORMALIZED_PATH
        if not path.exists():
            print(f"Building MSK normalized intermediate: {path}", flush=True)
            write_msk_normalized_records(path, sample_limit=sample_limit)
        paths.append(path)

    records = []
    for path in paths:
        records.extend(read_tsv(path))
    return records


def load_evidence_sets() -> dict[str, set[str]]:
    with VARIANTS_PATH.open("r", encoding="utf-8") as f:
        variants = json.load(f)

    evidence_sets = {source: set() for source in SOURCES}
    evidence_sets["All"] = set()
    for variant in variants:
        name = variant.get("display_name", "")
        if not name:
            continue
        evidence_sets["All"].add(name)
        for source in variant.get("source", []):
            if source in evidence_sets:
                evidence_sets[source].add(name)
    return evidence_sets


def load_relation_graph() -> VariationHierarchicalRelation:
    if not HIERARCHY_PATH.exists():
        return VariationHierarchicalRelation([])
    return VariationHierarchicalRelation(load_hierarchy(HIERARCHY_PATH), reduce=False)


def evaluate_records(records: list[dict[str, str]]) -> list[dict[str, Any]]:
    print(f"Evaluating patient variant records: {len(records)}", flush=True)
    evidence_sets = load_evidence_sets()
    graph = load_relation_graph()

    rows = []
    for index, row in enumerate(records, start=1):
        if index % 50000 == 0:
            print(f"Evaluated {index}/{len(records)} records", flush=True)

        normalized = row.get("normalized_variant", "")
        ancestors = graph.find_ancestor(normalized) if normalized else []
        exact_sources = [source for source in SOURCES if normalized in evidence_sets[source]]
        hierarchy_sources = [
            source
            for source in SOURCES
            if normalized in evidence_sets[source] or any(parent in evidence_sets[source] for parent in ancestors)
        ]
        matched_ancestors = sorted(parent for parent in ancestors if parent in evidence_sets["All"])

        rows.append(
            {
                **row,
                "normalized_variant": normalized,
                "parse_status": "parsed" if normalized else "unparsed",
                "exact_sources": ";".join(exact_sources),
                "hierarchy_sources": ";".join(hierarchy_sources),
                "is_exact_matched": int(bool(exact_sources)),
                "is_hierarchy_matched": int(bool(hierarchy_sources)),
                "is_hierarchy_added": int(bool(hierarchy_sources) and not bool(exact_sources)),
                "matched_ancestors": ";".join(matched_ancestors),
            }
        )
    return rows


def summarize_by_cohort_source(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary = []
    for cohort in sorted({row["cohort"] for row in rows}):
        cohort_rows = [row for row in rows if row["cohort"] == cohort and row["parse_status"] == "parsed"]
        for source in ["All"] + SOURCES:
            if source == "All":
                exact = sum(int(row["is_exact_matched"]) for row in cohort_rows)
                hierarchy = sum(int(row["is_hierarchy_matched"]) for row in cohort_rows)
            else:
                exact = sum(source in row["exact_sources"].split(";") for row in cohort_rows)
                hierarchy = sum(source in row["hierarchy_sources"].split(";") for row in cohort_rows)
            total = len(cohort_rows)
            summary.append(
                {
                    "cohort": cohort,
                    "evidence_source": source,
                    "parsed_variant_occurrences": total,
                    "exact_matched": exact,
                    "hierarchy_matched": hierarchy,
                    "hierarchy_added": hierarchy - exact,
                    "exact_rate": f"{exact / total:.4f}" if total else "0.0000",
                    "hierarchy_rate": f"{hierarchy / total:.4f}" if total else "0.0000",
                }
            )
    return summary


def summarize_by_cohort(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary = []
    for cohort in sorted({row["cohort"] for row in rows}):
        cohort_rows = [row for row in rows if row["cohort"] == cohort]
        parsed_rows = [row for row in cohort_rows if row["parse_status"] == "parsed"]
        source_exact = sum(len([s for s in row["exact_sources"].split(";") if s]) for row in parsed_rows)
        source_hierarchy = sum(len([s for s in row["hierarchy_sources"].split(";") if s]) for row in parsed_rows)
        summary.append(
            {
                "cohort": cohort,
                "variant_occurrences": len(cohort_rows),
                "parsed_variant_occurrences": len(parsed_rows),
                "unique_normalized_variants": len({row["normalized_variant"] for row in parsed_rows}),
                "exact_matched": sum(int(row["is_exact_matched"]) for row in parsed_rows),
                "hierarchy_matched": sum(int(row["is_hierarchy_matched"]) for row in parsed_rows),
                "hierarchy_added": sum(int(row["is_hierarchy_added"]) for row in parsed_rows),
                "source_level_exact_matches": source_exact,
                "source_level_hierarchy_matches": source_hierarchy,
                "source_level_hierarchy_added": source_hierarchy - source_exact,
            }
        )
    return summary


def summarize_by_sample(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["cohort"], row["sample_id"])].append(row)

    summary = []
    for (cohort, sample_id), sample_rows in sorted(grouped.items()):
        parsed_rows = [row for row in sample_rows if row["parse_status"] == "parsed"]
        summary.append(
            {
                "cohort": cohort,
                "sample_id": sample_id,
                "patient_id": sample_rows[0].get("patient_id", ""),
                "cancer_type": sample_rows[0].get("cancer_type", ""),
                "parsed_variant_occurrences": len(parsed_rows),
                "has_exact_match": int(any(int(row["is_exact_matched"]) for row in parsed_rows)),
                "has_hierarchy_match": int(any(int(row["is_hierarchy_matched"]) for row in parsed_rows)),
            }
        )
    return summary


def write_svg_chart(summary: list[dict[str, Any]], output_dir: Path) -> None:
    path = output_dir / "evidence_match_summary.svg"
    width = 900
    height = 300 + 80 * len(summary)
    max_value = max([int(row["source_level_hierarchy_matches"]) for row in summary] + [1])
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="24" y="34" font-family="Arial" font-size="20">Evidence matching: exact vs hierarchy-aware</text>',
    ]
    y = 70
    for row in summary:
        exact = int(row["source_level_exact_matches"])
        hierarchy = int(row["source_level_hierarchy_matches"])
        exact_w = int(620 * exact / max_value)
        hierarchy_w = int(620 * hierarchy / max_value)
        lines.append(f'<text x="24" y="{y + 18}" font-family="Arial" font-size="14">{row["cohort"]}</text>')
        lines.append(f'<rect x="220" y="{y}" width="{exact_w}" height="22" fill="#6b7280"/>')
        lines.append(f'<rect x="220" y="{y + 28}" width="{hierarchy_w}" height="22" fill="#2563eb"/>')
        lines.append(f'<text x="{230 + exact_w}" y="{y + 17}" font-family="Arial" font-size="13">{exact}</text>')
        lines.append(f'<text x="{230 + hierarchy_w}" y="{y + 45}" font-family="Arial" font-size="13">{hierarchy}</text>')
        y += 80
    lines.append('<text x="220" y="250" font-family="Arial" font-size="13" fill="#6b7280">gray: exact</text>')
    lines.append('<text x="330" y="250" font-family="Arial" font-size="13" fill="#2563eb">blue: hierarchy-aware</text>')
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def write_results(result_rows: list[dict[str, Any]], output_dir: Path) -> None:
    source_summary = summarize_by_cohort_source(result_rows)
    cohort_summary = summarize_by_cohort(result_rows)
    sample_summary = summarize_by_sample(result_rows)

    write_tsv(
        output_dir / "patient_variant_matches.tsv",
        result_rows,
        [
            "cohort",
            "patient_id",
            "sample_id",
            "cancer_type",
            "variant_source",
            "raw_variant",
            "normalized_variant",
            "parse_status",
            "exact_sources",
            "hierarchy_sources",
            "is_exact_matched",
            "is_hierarchy_matched",
            "is_hierarchy_added",
            "matched_ancestors",
        ],
    )
    write_tsv(
        output_dir / "summary_by_cohort_source.tsv",
        source_summary,
        [
            "cohort",
            "evidence_source",
            "parsed_variant_occurrences",
            "exact_matched",
            "hierarchy_matched",
            "hierarchy_added",
            "exact_rate",
            "hierarchy_rate",
        ],
    )
    write_tsv(
        output_dir / "summary_by_cohort.tsv",
        cohort_summary,
        [
            "cohort",
            "variant_occurrences",
            "parsed_variant_occurrences",
            "unique_normalized_variants",
            "exact_matched",
            "hierarchy_matched",
            "hierarchy_added",
            "source_level_exact_matches",
            "source_level_hierarchy_matches",
            "source_level_hierarchy_added",
        ],
    )
    write_tsv(
        output_dir / "summary_by_sample.tsv",
        sample_summary,
        [
            "cohort",
            "sample_id",
            "patient_id",
            "cancer_type",
            "parsed_variant_occurrences",
            "has_exact_match",
            "has_hierarchy_match",
        ],
    )
    write_svg_chart(cohort_summary, output_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate patient-level evidence matching.")
    parser.add_argument("--cohort", choices=["all", "genie", "msk"], default="all")
    parser.add_argument("--sample-limit", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(f"Loading patient records: {args.cohort}", flush=True)
    records = load_patient_records(args.cohort, sample_limit=args.sample_limit or None)
    print(f"Loaded unique patient records: {len(records)}", flush=True)

    result_rows = evaluate_records(records)
    output_dir = RESULT_DIR / "smoke" if args.sample_limit else RESULT_DIR

    print("Writing result files", flush=True)
    write_results(result_rows, output_dir)
    print(f"Patient variant records: {len(result_rows)}")
    print(f"Results written to: {output_dir}")


if __name__ == "__main__":
    main()
