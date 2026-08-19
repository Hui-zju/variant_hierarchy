from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from evaluate.common import PROJECT_ROOT, add_project_paths, read_tsv, write_tsv


CHALLENGE_PATH = PROJECT_ROOT / "data" / "challenge_sets" / "hierarchy_relation_challenge_set.tsv"
HIERARCHY_PATH = PROJECT_ROOT / "data" / "hierarchy" / "hierarchy.json"
RESULT_DIR = Path(__file__).resolve().parent / "results"
RESULT_PATH = RESULT_DIR / "hierarchy_relation_results.tsv"
SUMMARY_PATH = RESULT_DIR / "hierarchy_summary_by_stratum.tsv"
CONFUSION_PATH = RESULT_DIR / "hierarchy_confusion_matrix.tsv"

add_project_paths()
from hierarchy_class import VariationHierarchicalRelation, load_hierarchy  # noqa: E402
from variant_hierarchy.contains import GenomicBackendUnavailable, contains_variant  # noqa: E402
from variant_hierarchy.normalize.variant_parser import parse_variant  # noqa: E402


def load_relation_graph() -> VariationHierarchicalRelation:
    if not HIERARCHY_PATH.exists():
        return VariationHierarchicalRelation([])
    return VariationHierarchicalRelation(load_hierarchy(HIERARCHY_PATH), reduce=False)


def parse_challenge_variant(name: str) -> tuple[dict[str, Any] | None, str, str]:
    variant = parse_variant(name)
    if variant:
        return variant, variant.get("display_name", ""), "parsed"
    return None, "", "unresolved"


def evaluate() -> list[dict[str, Any]]:
    results = []
    graph = load_relation_graph()
    for row in read_tsv(CHALLENGE_PATH):
        parent_variant, parent, parent_method = parse_challenge_variant(row.get("parent_variant", ""))
        child_variant, child, child_method = parse_challenge_variant(row.get("child_variant", ""))

        if parent_variant and child_variant:
            graph_contains = parent in graph.find_ancestor(child)
            try:
                logic_contains = contains_variant(parent_variant, child_variant)
            except GenomicBackendUnavailable:
                logic_contains = False
            predicted = "1" if graph_contains or logic_contains else "0"
        else:
            graph_contains = False
            logic_contains = False
            predicted = "unresolved"

        expected = row.get("expected_contains", "")
        results.append(
            {
                "id": row.get("id", ""),
                "stratum_id": row.get("stratum_id", ""),
                "parent_variant": row.get("parent_variant", ""),
                "child_variant": row.get("child_variant", ""),
                "expected_contains": expected,
                "predicted_contains": predicted,
                "is_correct": int(predicted == expected),
                "parent_resolved": parent,
                "child_resolved": child,
                "parent_resolution_method": parent_method,
                "child_resolution_method": child_method,
                "matched_by_graph": int(graph_contains),
                "matched_by_logic": int(logic_contains),
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


def confusion_matrix(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    counts = Counter((row["expected_contains"], row["predicted_contains"]) for row in rows)
    return [
        {"expected_contains": expected, "predicted_contains": predicted, "count": count}
        for (expected, predicted), count in sorted(counts.items())
    ]


def main() -> None:
    rows = evaluate()
    summary = summarize(rows)
    confusion = confusion_matrix(rows)

    write_tsv(
        RESULT_PATH,
        rows,
        [
            "id",
            "stratum_id",
            "parent_variant",
            "child_variant",
            "expected_contains",
            "predicted_contains",
            "is_correct",
            "parent_resolved",
            "child_resolved",
            "parent_resolution_method",
            "child_resolution_method",
            "matched_by_graph",
            "matched_by_logic",
        ],
    )
    write_tsv(SUMMARY_PATH, summary, ["stratum_id", "total", "correct", "accuracy"])
    write_tsv(CONFUSION_PATH, confusion, ["expected_contains", "predicted_contains", "count"])

    total = len(rows)
    correct = sum(int(row["is_correct"]) for row in rows)
    accuracy = correct / total if total else 0
    print(f"Hierarchy accuracy: {correct}/{total} ({accuracy:.1%})")
    print(f"Results written to: {RESULT_PATH}")
    print(f"Summary written to: {SUMMARY_PATH}")
    print(f"Confusion matrix written to: {CONFUSION_PATH}")


if __name__ == "__main__":
    main()
