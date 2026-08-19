from __future__ import annotations

import argparse
import re
from collections import Counter

from evaluate.common import add_project_paths
from evaluate.common import read_tsv
from evaluate.common import write_tsv
from evaluate.evidence_match.config import GENIE_DIR, GENIE_NORMALIZED_PATH, INTERMEDIATE_DIR


add_project_paths()
from variant_hierarchy.normalize.variant_parser import parse_variant  # noqa: E402


def clean_genie_variant_text(text: str) -> str:
    text = str(text or "").strip()
    text = re.sub(r"\s*-\s*Archer$", "", text)
    return re.sub(r"\s+", " ", text)


def split_genie_mutations(text: str) -> list[str]:
    return [clean_genie_variant_text(item) for item in str(text or "").split(";") if item.strip()]


def iter_genie_records(year: str = "all") -> list[dict[str, str]]:
    records = []
    years = ("2017", "2022") if year == "all" else (year,)
    for item_year in years:
        path = GENIE_DIR / f"genie_mskcc_samples_with_{item_year}_oncokb_annotation.txt"
        for row in read_tsv(path):
            for raw_variant in split_genie_mutations(row.get("ONCOGENIC_MUTATIONS", "")):
                records.append(
                    {
                        "cohort": f"GENIE_MSKCC_{item_year}",
                        "patient_id": row.get("PATIENT_ID", ""),
                        "sample_id": row.get("SAMPLE_ID", ""),
                        "cancer_type": row.get("CANCER_TYPE", ""),
                        "variant_source": "oncogenic_mutations",
                        "raw_variant": raw_variant,
                    }
                )
    return records


def filter_sample_limit(records: list[dict[str, str]], sample_limit: int | None = None) -> list[dict[str, str]]:
    if not sample_limit:
        return records

    selected = []
    sample_ids = []
    seen = set()
    for row in records:
        sample_id = row["sample_id"]
        if sample_id not in seen:
            if len(sample_ids) >= sample_limit:
                break
            seen.add(sample_id)
            sample_ids.append(sample_id)
        selected.append(row)
    return selected


def normalize_records(records: list[dict[str, str]]) -> list[dict[str, str]]:
    normalized_by_raw = {}
    for raw_variant in sorted({row["raw_variant"] for row in records}):
        parsed = parse_variant(raw_variant)
        normalized_by_raw[raw_variant] = parsed.get("display_name", "") if parsed else ""

    rows = []
    for row in records:
        normalized = normalized_by_raw[row["raw_variant"]]
        rows.append(
            {
                **row,
                "normalized_variant": normalized,
                "parse_status": "parsed" if normalized else "unparsed",
            }
        )
    return rows


def write_normalized_records(
    output_path=GENIE_NORMALIZED_PATH,
    year: str = "all",
    sample_limit: int | None = None,
) -> list[dict[str, str]]:
    records = filter_sample_limit(iter_genie_records(year), sample_limit)
    rows = normalize_records(records)
    write_tsv(
        output_path,
        rows,
        [
            "cohort",
            "patient_id",
            "sample_id",
            "cancer_type",
            "variant_source",
            "raw_variant",
            "normalized_variant",
            "parse_status",
        ],
    )
    return rows


def report_parse_status(year: str = "all", limit: int | None = None) -> None:
    records = iter_genie_records(year)
    raw_variants = sorted({row["raw_variant"] for row in records})
    if limit:
        raw_variants = raw_variants[:limit]

    counts = Counter()
    for raw_variant in raw_variants:
        parsed = parse_variant(raw_variant)
        if parsed:
            counts["parsed"] += 1
            print(f"parsed\t{raw_variant}\t{parsed.get('display_name', '')}")
        else:
            counts["unparsed"] += 1
            print(f"unparsed\t{raw_variant}\t")
    print(f"summary\tparsed={counts['parsed']}\tunparsed={counts['unparsed']}\ttotal={sum(counts.values())}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize GENIE/MSKCC patient variants.")
    parser.add_argument("--year", choices=["all", "2017", "2022"], default="2022")
    parser.add_argument("--sample-limit", type=int, default=0)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.report_only:
        report_parse_status(year=args.year, limit=args.limit or None)
    else:
        output_path = (
            INTERMEDIATE_DIR / "smoke" / GENIE_NORMALIZED_PATH.name
            if args.sample_limit
            else GENIE_NORMALIZED_PATH
        )
        rows = write_normalized_records(output_path, year=args.year, sample_limit=args.sample_limit or None)
        counts = Counter(row["parse_status"] for row in rows)
        print(
            f"written={output_path}\t"
            f"rows={len(rows)}\tparsed={counts['parsed']}\tunparsed={counts['unparsed']}"
        )
