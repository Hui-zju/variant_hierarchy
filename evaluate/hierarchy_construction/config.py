from __future__ import annotations

from pathlib import Path

from evaluate.common import PROJECT_ROOT


DATA_DIR = PROJECT_ROOT / "data"
VARIANTS_PATH = DATA_DIR / "variants" / "variants.json"
FULL_HIERARCHY_PATH = DATA_DIR / "hierarchy" / "hierarchy.json"
NON_REDUNDANT_HIERARCHY_PATH = DATA_DIR / "hierarchy" / "hierarchy_non_redundant.json"
RESULT_DIR = Path(__file__).resolve().parent / "results"

SOURCES = ["OncoKB", "CGI", "CIVIC", "COSMIC"]
CLASS_LABELS = {
    "HGVS Variant": "p.HGVS variant",
    "Sequence Variant": "Categorical sequence variant",
    "Function Variant": "Functional assertion category",
    "Copy Number Variant": "Copy-number alteration category",
    "Structural Variant": "Structural alteration category",
    "Expression Variant": "Expression-state category",
    "Epigenetic Variant": "Methylation-state category",
    "Uncertain Variant": "Unresolved variant term",
    "Unclassified": "Unclassified term",
}
RELATION_SOURCE_GROUPS = [
    ("semantic", "Semantically inferred relationships", {"VH", "RBI"}),
    (
        "knowledge",
        "Knowledge-supported relationships",
        {"CGI", "OncoKB-SOP", "OncoKB", "Oncokb", "PMC4232638"},
    ),
]
COLORS = {
    "blue": "#2563eb",
    "blue_light": "#93c5fd",
    "green": "#059669",
    "green_light": "#86efac",
    "gray": "#6b7280",
    "gray_light": "#d1d5db",
    "orange": "#ea580c",
    "purple": "#7c3aed",
    "red": "#dc2626",
    "ink": "#111827",
    "muted": "#4b5563",
    "grid": "#e5e7eb",
    "paper": "#ffffff",
}


def manuscript_class_label(value: str) -> str:
    return CLASS_LABELS.get(value, value)


def source_set_sort_key(item: tuple[tuple[str, ...], int]) -> tuple[int, int, tuple[str, ...]]:
    source_set, value = item
    if not source_set:
        return (99, 0, source_set)
    return (len(source_set), -value, source_set)
