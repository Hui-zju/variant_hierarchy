from __future__ import annotations

import csv
import argparse
import json
import re
import urllib.request
from collections import Counter
from pathlib import Path

from evaluate.common import add_project_paths
from evaluate.common import read_tsv
from evaluate.common import write_tsv
from evaluate.evidence_match.config import (
    CBIOPORTAL_MUTATION_URL,
    INTERMEDIATE_DIR,
    MSK_DIR,
    MSK_NORMALIZED_PATH,
    MSK_MUTATION_API_JSON_PATH,
    MSK_MUTATION_API_PATH,
)


add_project_paths()
from variant_hierarchy.normalize.variant_parser import parse_variant  # noqa: E402


def clean_msk_variant_text(text: str) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip())


def load_msk_sample_info() -> dict[str, dict[str, str]]:
    sample_path = MSK_DIR / "data_clinical_sample.txt"
    return {row["SAMPLE_ID"]: row for row in read_tsv(sample_path, skip_comments=True)}


def download_msk_mutations_from_api(sample_ids: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["sample_id", "patient_id", "gene", "protein_change", "mutation_type"]

    with output_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i in range(0, len(sample_ids), 200):
            body = json.dumps({"sampleIds": sample_ids[i : i + 200]}).encode("utf-8")
            request = urllib.request.Request(
                CBIOPORTAL_MUTATION_URL,
                data=body,
                headers={"Accept": "application/json", "Content-Type": "application/json"},
                method="POST",
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                mutations = json.loads(response.read().decode("utf-8"))
            for item in mutations:
                writer.writerow(
                    {
                        "sample_id": item.get("sampleId", ""),
                        "patient_id": item.get("patientId", ""),
                        "gene": item.get("gene", {}).get("hugoGeneSymbol", ""),
                        "protein_change": item.get("proteinChange", ""),
                        "mutation_type": item.get("mutationType", ""),
                    }
                )
            print(f"Downloaded MSK mutations: {min(i + 200, len(sample_ids))}/{len(sample_ids)} samples")


def mutation_type_to_text(mutation_type: str) -> str:
    mutation_type = mutation_type.lower()
    if "splice" in mutation_type:
        return "splice mutation"
    if "frame_shift" in mutation_type or "nonsense" in mutation_type or "nonstop" in mutation_type:
        return "truncating mutation"
    if "missense" in mutation_type:
        return "missense mutation"
    if "in_frame_del" in mutation_type:
        return "inframe deletion"
    if "in_frame_ins" in mutation_type:
        return "inframe insertion"
    return "mutation"


def load_msk_mutation_rows(sample_info: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    if MSK_MUTATION_API_JSON_PATH.exists():
        with MSK_MUTATION_API_JSON_PATH.open("r", encoding="utf-8") as f:
            api_rows = json.load(f)
        return [
            {
                "sample_id": item.get("sampleId", ""),
                "patient_id": item.get("patientId", ""),
                "gene": item.get("gene", {}).get("hugoGeneSymbol", ""),
                "protein_change": item.get("proteinChange", ""),
                "mutation_type": item.get("mutationType", ""),
            }
            for item in api_rows
        ]

    if not MSK_MUTATION_API_PATH.exists():
        download_msk_mutations_from_api(sorted(sample_info), MSK_MUTATION_API_PATH)
    return read_tsv(MSK_MUTATION_API_PATH)


def iter_msk_mutation_records(sample_info: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    records = []
    for row in load_msk_mutation_rows(sample_info):
        sample_id = row.get("sample_id", "")
        info = sample_info.get(sample_id, {})
        gene = row.get("gene", "")
        protein_change = clean_msk_variant_text(row.get("protein_change", ""))
        if protein_change:
            raw_variant = f"{gene} {protein_change}"
        else:
            raw_variant = f"{gene} {mutation_type_to_text(row.get('mutation_type', ''))}"
        records.append(
            {
                "cohort": "MSK_IMPACT_2017",
                "patient_id": row.get("patient_id", "") or info.get("PATIENT_ID", ""),
                "sample_id": sample_id,
                "cancer_type": info.get("CANCER_TYPE", ""),
                "variant_source": "mutation",
                "raw_variant": raw_variant,
            }
        )
    return records


def iter_msk_cna_records(sample_info: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    path = MSK_DIR / "data_cna.txt"
    records = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        sample_ids = header[1:]
        for row in reader:
            gene = row[0]
            for sample_id, value in zip(sample_ids, row[1:]):
                if value not in {"2", "-2"}:
                    continue
                info = sample_info.get(sample_id, {})
                raw_variant = f"{gene} amplification" if value == "2" else f"{gene} copy number deletion"
                records.append(
                    {
                        "cohort": "MSK_IMPACT_2017",
                        "patient_id": info.get("PATIENT_ID", ""),
                        "sample_id": sample_id,
                        "cancer_type": info.get("CANCER_TYPE", ""),
                        "variant_source": "cna",
                        "raw_variant": raw_variant,
                    }
                )
    return records


def iter_msk_sv_records(sample_info: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    records = []
    for row in read_tsv(MSK_DIR / "data_sv.txt"):
        sample_id = row.get("Sample_Id", "")
        info = sample_info.get(sample_id, {})
        g1 = row.get("Site1_Hugo_Symbol", "")
        g2 = row.get("Site2_Hugo_Symbol", "")
        event_info = row.get("Event_Info", "")
        if g1 and g2 and "fusion" in event_info.lower():
            raw_variant = f"{g1}::{g2} fusion"
        elif event_info:
            raw_variant = event_info
        else:
            continue
        records.append(
            {
                "cohort": "MSK_IMPACT_2017",
                "patient_id": info.get("PATIENT_ID", ""),
                "sample_id": sample_id,
                "cancer_type": info.get("CANCER_TYPE", ""),
                "variant_source": "sv",
                "raw_variant": clean_msk_variant_text(raw_variant),
            }
        )
    return records


def iter_msk_records() -> list[dict[str, str]]:
    if not (MSK_DIR / "data_clinical_sample.txt").exists():
        return []
    sample_info = load_msk_sample_info()
    return (
        iter_msk_mutation_records(sample_info)
        + iter_msk_cna_records(sample_info)
        + iter_msk_sv_records(sample_info)
    )


def filter_sample_limit(records: list[dict[str, str]], sample_limit: int | None = None) -> list[dict[str, str]]:
    if not sample_limit:
        return records

    selected_sample_ids = []
    seen = set()
    for row in records:
        sample_id = row["sample_id"]
        if sample_id not in seen:
            if len(selected_sample_ids) >= sample_limit:
                break
            seen.add(sample_id)
            selected_sample_ids.append(sample_id)
    selected = set(selected_sample_ids)
    return [row for row in records if row["sample_id"] in selected]


def normalize_records(records: list[dict[str, str]]) -> list[dict[str, str]]:
    normalized_by_raw = {}
    for raw_variant in sorted({row["raw_variant"] for row in records}):
        parsed = parse_variant(raw_variant)
        normalized_by_raw[raw_variant] = parsed.get("display_name", "") if parsed else ""

    rows = []
    for index, row in enumerate(records):
        print(f"{index} finished!")
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
    output_path=MSK_NORMALIZED_PATH,
    sample_limit: int | None = None,
) -> list[dict[str, str]]:
    records = filter_sample_limit(iter_msk_records(), sample_limit)
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


def report_parse_status(limit: int | None = None) -> None:
    records = iter_msk_records()
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
    parser = argparse.ArgumentParser(description="Normalize MSK-IMPACT patient variants.")
    parser.add_argument("--sample-limit", type=int, default=0)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.report_only:
        report_parse_status(limit=args.limit or None)
    else:
        output_path = (
            INTERMEDIATE_DIR / "smoke" / MSK_NORMALIZED_PATH.name
            if args.sample_limit
            else MSK_NORMALIZED_PATH
        )
        rows = write_normalized_records(output_path, sample_limit=args.sample_limit or None)
        counts = Counter(row["parse_status"] for row in rows)
        print(
            f"written={output_path}\t"
            f"rows={len(rows)}\tparsed={counts['parsed']}\tunparsed={counts['unparsed']}"
        )
