from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

from evaluate.common import add_project_paths, write_tsv
from evaluate.hierarchy_construction.config import (
    FULL_HIERARCHY_PATH,
    NON_REDUNDANT_HIERARCHY_PATH,
    RELATION_SOURCE_GROUPS,
    RESULT_DIR,
    SOURCES,
    VARIANTS_PATH,
    manuscript_class_label,
    source_set_sort_key,
)


add_project_paths()
from variant_hierarchy.normalize.variant_parser import parse_variant  # noqa: E402


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def edge_name(edge: dict[str, Any], role: str) -> str:
    value = edge.get(role)
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return value.get("display_name") or value.get("displayName") or ""
    return ""


def build_variant_info(variants: list[dict[str, Any]], full_edges: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    info = {
        variant["display_name"]: variant
        for variant in variants
        if variant.get("display_name")
    }
    for edge in full_edges:
        for role in ("parent", "child"):
            value = edge.get(role)
            if isinstance(value, dict):
                name = value.get("display_name") or value.get("displayName")
                if name:
                    info.setdefault(name, value)
    return info


def get_variant_info(name: str, info: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    if name not in info:
        parsed = parse_variant(name)
        if parsed:
            info[name] = parsed
    return info.get(name)


def variant_class(name: str, info: dict[str, dict[str, Any]]) -> str:
    variant = get_variant_info(name, info)
    return (variant or {}).get("variant_class", "Unclassified")


def count_edge_categories(
    edges: list[dict[str, Any]],
    info: dict[str, dict[str, Any]],
) -> Counter[tuple[str, str]]:
    counts: Counter[tuple[str, str]] = Counter()
    for edge in edges:
        parent = variant_class(edge_name(edge, "parent"), info)
        child = variant_class(edge_name(edge, "child"), info)
        counts[(parent, child)] += 1
    return counts


def count_relation_sources(edges: list[dict[str, Any]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for edge in edges:
        source = edge.get("source", [])
        if isinstance(source, list):
            for item in source:
                counts[str(item)] += 1
        elif isinstance(source, dict):
            counts["dict_source"] += 1
        else:
            counts[str(source)] += 1
    return counts


def group_relation_sources(source_counts: Counter[str]) -> Counter[str]:
    grouped: Counter[str] = Counter()
    assigned_sources = set()
    for group_key, _, source_keys in RELATION_SOURCE_GROUPS:
        grouped[group_key] = sum(source_counts.get(source_key, 0) for source_key in source_keys)
        assigned_sources.update(source_keys)
    for source_key, count in source_counts.items():
        if source_key not in assigned_sources:
            grouped[source_key] += count
    return grouped


def pct(value: int, total: int) -> str:
    return f"{value / total * 100:.1f}" if total else "0.0"


def collect_stats() -> dict[str, Any]:
    variants = load_json(VARIANTS_PATH)
    full_edges = load_json(FULL_HIERARCHY_PATH)
    non_edges = load_json(NON_REDUNDANT_HIERARCHY_PATH)
    info = build_variant_info(variants, full_edges)

    variant_class_counts = Counter(variant.get("variant_class", "") for variant in variants)
    source_counts = Counter(source for variant in variants for source in variant.get("source", []))
    source_set_counts = Counter(tuple(sorted(variant.get("source", []))) for variant in variants)

    pairwise_counts = {}
    for left, right in combinations(SOURCES, 2):
        pairwise_counts[(left, right)] = sum(
            1
            for variant in variants
            if left in variant.get("source", []) and right in variant.get("source", [])
        )

    full_edge_categories = count_edge_categories(full_edges, info)
    non_edge_categories = count_edge_categories(non_edges, info)
    full_relation_sources = count_relation_sources(full_edges)
    non_relation_sources = count_relation_sources(non_edges)
    full_relation_source_groups = group_relation_sources(full_relation_sources)
    non_relation_source_groups = group_relation_sources(non_relation_sources)

    def nodes(edges: list[dict[str, Any]]) -> set[str]:
        return {
            name
            for edge in edges
            for name in (edge_name(edge, "parent"), edge_name(edge, "child"))
            if name
        }

    return {
        "variants": variants,
        "full_edges": full_edges,
        "non_edges": non_edges,
        "variant_class_counts": variant_class_counts,
        "source_counts": source_counts,
        "source_set_counts": source_set_counts,
        "pairwise_counts": pairwise_counts,
        "full_edge_categories": full_edge_categories,
        "non_edge_categories": non_edge_categories,
        "full_relation_sources": full_relation_sources,
        "non_relation_sources": non_relation_sources,
        "full_relation_source_groups": full_relation_source_groups,
        "non_relation_source_groups": non_relation_source_groups,
        "full_nodes": nodes(full_edges),
        "non_nodes": nodes(non_edges),
    }


def write_tables(stats: dict[str, Any]) -> None:
    total_variants = len(stats["variants"])
    full_edges = stats["full_edges"]
    non_edges = stats["non_edges"]
    removed_edges = len(full_edges) - len(non_edges)

    write_tsv(
        RESULT_DIR / "variant_class_distribution.tsv",
        [
            {
                "semantic_category": manuscript_class_label(key),
                "internal_variant_class": key,
                "variant_count": value,
                "percentage": pct(value, total_variants),
            }
            for key, value in stats["variant_class_counts"].most_common()
        ],
        ["semantic_category", "internal_variant_class", "variant_count", "percentage"],
    )

    write_tsv(
        RESULT_DIR / "source_distribution.tsv",
        [
            {"source": key, "variant_count": stats["source_counts"].get(key, 0)}
            for key in SOURCES
        ],
        ["source", "variant_count"],
    )

    write_tsv(
        RESULT_DIR / "source_intersections.tsv",
        [
            {
                "source_set": "+".join(key) if key else "No source",
                "variant_count": value,
                "source_count": len(key),
            }
            for key, value in sorted(stats["source_set_counts"].items(), key=source_set_sort_key)
        ],
        ["source_set", "source_count", "variant_count"],
    )

    write_tsv(
        RESULT_DIR / "pairwise_source_overlap.tsv",
        [
            {"source_1": left, "source_2": right, "shared_variants": value}
            for (left, right), value in stats["pairwise_counts"].items()
        ],
        ["source_1", "source_2", "shared_variants"],
    )

    write_tsv(
        RESULT_DIR / "graph_summary.tsv",
        [
            {"statistic": "Full hierarchy nodes", "value": len(stats["full_nodes"])},
            {"statistic": "Full hierarchy edges", "value": len(full_edges)},
            {"statistic": "Non-redundant hierarchy nodes", "value": len(stats["non_nodes"])},
            {"statistic": "Non-redundant hierarchy edges", "value": len(non_edges)},
            {"statistic": "Removed redundant edges", "value": removed_edges},
            {"statistic": "Edge reduction (%)", "value": f"{removed_edges / len(full_edges) * 100:.1f}"},
        ],
        ["statistic", "value"],
    )

    category_rows = []
    for key in sorted(
        set(stats["full_edge_categories"]) | set(stats["non_edge_categories"]),
        key=lambda item: (
            stats["non_edge_categories"].get(item, 0),
            stats["full_edge_categories"].get(item, 0),
            item,
        ),
        reverse=True,
    ):
        full_count = stats["full_edge_categories"].get(key, 0)
        non_count = stats["non_edge_categories"].get(key, 0)
        category_rows.append(
            {
                "parent_category": manuscript_class_label(key[0]),
                "child_category": manuscript_class_label(key[1]),
                "internal_parent_class": key[0],
                "internal_child_class": key[1],
                "full_hierarchy_edges": full_count,
                "non_redundant_edges": non_count,
                "removed_redundant_edges": full_count - non_count,
            }
        )
    write_tsv(
        RESULT_DIR / "edge_category_summary.tsv",
        category_rows,
        [
            "parent_category",
            "child_category",
            "internal_parent_class",
            "internal_child_class",
            "full_hierarchy_edges",
            "non_redundant_edges",
            "removed_redundant_edges",
        ],
    )

    write_tsv(
        RESULT_DIR / "edge_relation_sources.tsv",
        [
            {
                "relation_source": label,
                "full_hierarchy_edges": stats["full_relation_source_groups"].get(key, 0),
                "non_redundant_edges": stats["non_relation_source_groups"].get(key, 0),
            }
            for key, label, _ in RELATION_SOURCE_GROUPS
        ],
        ["relation_source", "full_hierarchy_edges", "non_redundant_edges"],
    )


def main() -> None:
    from evaluate.hierarchy_construction.plot import write_figures

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    stats = collect_stats()
    write_tables(stats)
    write_figures(stats)
    print(f"Variants: {len(stats['variants'])}")
    print(f"Full hierarchy edges: {len(stats['full_edges'])}")
    print(f"Non-redundant hierarchy edges: {len(stats['non_edges'])}")
    print(f"Results written to: {RESULT_DIR}")


if __name__ == "__main__":
    main()
