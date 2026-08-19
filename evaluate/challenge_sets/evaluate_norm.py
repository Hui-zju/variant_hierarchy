from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

from evaluate.common import PROJECT_ROOT, add_project_paths, read_tsv, write_tsv


CHALLENGE_PATH = PROJECT_ROOT / "data" / "challenge_sets" / "normalization_challenge_set.tsv"
RESULT_DIR = Path(__file__).resolve().parent / "results"
RESULT_PATH = RESULT_DIR / "normalization_results.tsv"
SUMMARY_PATH = RESULT_DIR / "normalization_summary_by_stratum.tsv"

add_project_paths()
from variant_hierarchy.normalize.variant_parser import parse_variant  # noqa: E402


def evaluate() -> list[dict[str, Any]]:
    results = []
    for row in read_tsv(CHALLENGE_PATH):
        expected = row.get("expected_variant_name") or ""
        input_text = row.get("input_text") or ""

        try:
            parsed = parse_variant(input_text)
            predicted = parsed.get("display_name", "") if parsed else ""
            error = ""
        except Exception as exc:
            predicted = ""
            error = f"{type(exc).__name__}: {exc}"

        if error:
            status = "exception"
        elif not predicted:
            status = "unparsed"
        elif predicted == expected:
            status = "correct"
        else:
            status = "incorrect"

        results.append(
            {
                "id": row.get("id", ""),
                "stratum_id": row.get("stratum_id", ""),
                "input_text": input_text,
                "expected_variant_name": expected,
                "predicted_variant_name": predicted,
                "status": status,
                "is_correct": int(status == "correct"),
                "error": error,
            }
        )
    return results


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[row["stratum_id"]].append(row)

    summary = []
    for stratum_id in sorted(groups):
        total = len(groups[stratum_id])
        correct = sum(int(row["is_correct"]) for row in groups[stratum_id])
        summary.append(
            {
                "stratum_id": stratum_id,
                "total": total,
                "correct": correct,
                "accuracy": f"{correct / total:.4f}" if total else "0.0000",
            }
        )
    return summary


def main() -> None:
    rows = evaluate()
    summary = summarize(rows)

    write_tsv(
        RESULT_PATH,
        rows,
        [
            "id",
            "stratum_id",
            "input_text",
            "expected_variant_name",
            "predicted_variant_name",
            "status",
            "is_correct",
            "error",
        ],
    )
    write_tsv(SUMMARY_PATH, summary, ["stratum_id", "total", "correct", "accuracy"])

    total = len(rows)
    correct = sum(int(row["is_correct"]) for row in rows)
    accuracy = correct / total if total else 0
    print(f"Normalization accuracy: {correct}/{total} ({accuracy:.1%})")
    print(f"Results written to: {RESULT_PATH}")
    print(f"Summary written to: {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
