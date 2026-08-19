from __future__ import annotations

import html
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from evaluate.common import read_tsv, write_tsv
from evaluate.evidence_match.config import RESULT_DIR, SOURCES, VARIANTS_PATH


MATCHES_PATH = RESULT_DIR / "patient_variant_matches.tsv"
COHORT_SUMMARY_PATH = RESULT_DIR / "summary_by_cohort.tsv"
SOURCE_SUMMARY_PATH = RESULT_DIR / "summary_by_cohort_source.tsv"
FIGURE_PATH = RESULT_DIR / "figure6_evidence_retrieval.pdf"
LOCATION_CATEGORY_SUMMARY_PATH = RESULT_DIR / "figure6_location_semantic_category_added_links.tsv"

COHORT_LABELS = {
    "GENIE_MSKCC_2022": "GENIE/MSKCC",
    "MSK_IMPACT_2017": "MSK-IMPACT 2017",
}
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
CLASS_PRIORITY = {
    "HGVS Variant": 0,
    "Sequence Variant": 1,
    "Copy Number Variant": 2,
    "Structural Variant": 3,
    "Expression Variant": 4,
    "Epigenetic Variant": 5,
    "Function Variant": 6,
    "Uncertain Variant": 7,
    "Unclassified": 8,
}
COLORS = {
    "blue": "#2563eb",
    "blue_light": "#93c5fd",
    "green": "#059669",
    "green_light": "#86efac",
    "gray": "#6b7280",
    "gray_light": "#d1d5db",
    "orange": "#ea580c",
    "purple": "#7c3aed",
    "ink": "#111827",
    "muted": "#4b5563",
    "grid": "#e5e7eb",
    "paper": "#ffffff",
}
SERIES = [
    ("GENIE_MSKCC_2022", "GENIE/MSKCC", COLORS["green"]),
    ("MSK_IMPACT_2017", "MSK-IMPACT 2017", COLORS["orange"]),
]
NORMALIZATION_SERIES = [("Normalized", COLORS["blue"]), ("Unnormalized", COLORS["gray_light"])]


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def svg_tag(name: str, attrs: dict[str, object], content: str | None = None) -> str:
    attr = " ".join(f'{key}="{esc(value)}"' for key, value in attrs.items() if value is not None)
    if content is None:
        return f"<{name} {attr}/>"
    return f"<{name} {attr}>{content}</{name}>"


def text(
    x: float,
    y: float,
    value: object,
    size: float = 16,
    weight: str = "400",
    fill: str = COLORS["ink"],
    anchor: str = "start",
) -> str:
    return svg_tag(
        "text",
        {
            "x": f"{x:.1f}",
            "y": f"{y:.1f}",
            "font-family": "Arial, Helvetica, sans-serif",
            "font-size": size,
            "font-weight": weight,
            "fill": fill,
            "text-anchor": anchor,
        },
        esc(value),
    )


def multiline_text(
    x: float,
    y: float,
    lines: list[str],
    size: float = 14,
    fill: str = COLORS["muted"],
    anchor: str = "start",
    line_height: float = 18,
) -> list[str]:
    return [
        svg_tag(
            "text",
            {
                "x": f"{x:.1f}",
                "y": f"{y + i * line_height:.1f}",
                "font-family": "Arial, Helvetica, sans-serif",
                "font-size": size,
                "fill": fill,
                "text-anchor": anchor,
            },
            esc(line),
        )
        for i, line in enumerate(lines)
    ]


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str | None = None, rx: float = 0) -> str:
    return svg_tag(
        "rect",
        {
            "x": f"{x:.1f}",
            "y": f"{y:.1f}",
            "width": f"{max(w, 0):.1f}",
            "height": f"{h:.1f}",
            "fill": fill,
            "stroke": stroke,
            "stroke-width": 1 if stroke else None,
            "rx": rx or None,
        },
    )


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = COLORS["grid"], width: float = 1) -> str:
    return svg_tag(
        "line",
        {
            "x1": f"{x1:.1f}",
            "y1": f"{y1:.1f}",
            "x2": f"{x2:.1f}",
            "y2": f"{y2:.1f}",
            "stroke": stroke,
            "stroke-width": width,
        },
    )


def wrap_label(label: str, max_chars: int = 31) -> list[str]:
    words = label.split()
    lines: list[str] = []
    current = ""
    for word in words:
        if len(current) + len(word) + bool(current) <= max_chars:
            current = f"{current} {word}".strip()
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines[:3]


def label_for_cohort(cohort: str) -> str:
    return COHORT_LABELS.get(cohort, cohort)


def source_set(value: str) -> set[str]:
    return {item for item in value.split(";") if item}


def load_variant_maps() -> tuple[dict[str, dict[str, Any]], dict[str, set[str]]]:
    with VARIANTS_PATH.open("r", encoding="utf-8") as f:
        variants = json.load(f)
    info = {item["display_name"]: item for item in variants if item.get("display_name")}
    source_to_variants = {source: set() for source in SOURCES}
    for item in variants:
        name = item.get("display_name", "")
        for source in item.get("source", []):
            if source in source_to_variants:
                source_to_variants[source].add(name)
    return info, source_to_variants


def category_label(name: str, info: dict[str, dict[str, Any]]) -> str:
    variant_class = info.get(name, {}).get("variant_class", "Unclassified")
    return CLASS_LABELS.get(variant_class, variant_class)


def location_category_label(name: str, info: dict[str, dict[str, Any]]) -> str:
    variant = info.get(name, {})
    variant_class = variant.get("variant_class", "Unclassified")
    location_type = variant.get("location_type", "")
    location = {
        "Gene": "gene",
        "Exon": "exon",
        "Amino Acid": "amino-acid",
    }.get(location_type)
    if variant_class == "Sequence Variant" and location:
        return f"Categorical {location} sequence variant"
    return CLASS_LABELS.get(variant_class, variant_class)


def category_rank(name: str, info: dict[str, dict[str, Any]]) -> tuple[int, int, str]:
    variant_class = info.get(name, {}).get("variant_class", "Unclassified")
    location_type = info.get(name, {}).get("location_type", "")
    location_priority = {"Amino Acid": 0, "Exon": 1, "Gene": 2}.get(location_type, 3)
    return CLASS_PRIORITY.get(variant_class, CLASS_PRIORITY["Unclassified"]), location_priority, name


def collect_added_categories() -> list[dict[str, Any]]:
    info, source_to_variants = load_variant_maps()
    counts: Counter[tuple[str, str]] = Counter()
    rows = read_tsv(MATCHES_PATH)
    for row in rows:
        if row.get("parse_status") != "parsed":
            continue
        added_sources = source_set(row.get("hierarchy_sources", "")) - source_set(row.get("exact_sources", ""))
        if not added_sources:
            continue
        ancestors = [item for item in row.get("matched_ancestors", "").split(";") if item]
        for source in added_sources:
            candidates = [item for item in ancestors if item in source_to_variants.get(source, set())]
            if candidates:
                chosen = sorted(candidates, key=lambda item: category_rank(item, info))[0]
                category = location_category_label(chosen, info)
            else:
                category = CLASS_LABELS["Unclassified"]
            counts[(row["cohort"], category)] += 1

    cohort_order = list(COHORT_LABELS)
    category_order = sorted(
        {category for _, category in counts},
        key=lambda category: (
            -sum(counts[(cohort, category)] for cohort in cohort_order),
            category,
        ),
    )
    return [
        {
            "cohort": cohort,
            "cohort_label": label_for_cohort(cohort),
            "semantic_category": category,
            "additional_evidence_links": counts[(cohort, category)],
        }
        for category in category_order
        for cohort in cohort_order
        if counts[(cohort, category)]
    ]


def collect_figure_data() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    cohort_rows = read_tsv(COHORT_SUMMARY_PATH)
    overall = {"cohort": "Overall"}
    for key in (
        "variant_occurrences",
        "parsed_variant_occurrences",
        "unique_normalized_variants",
        "exact_matched",
        "hierarchy_matched",
        "hierarchy_added",
        "source_level_exact_matches",
        "source_level_hierarchy_matches",
        "source_level_hierarchy_added",
    ):
        overall[key] = sum(int(row[key]) for row in cohort_rows)
    cohort_rows = cohort_rows + [overall]
    source_rows = [row for row in read_tsv(SOURCE_SUMMARY_PATH) if row["evidence_source"] != "All"]
    category_rows = collect_added_categories()

    write_tsv(
        LOCATION_CATEGORY_SUMMARY_PATH,
        category_rows,
        ["cohort", "cohort_label", "semantic_category", "additional_evidence_links"],
    )
    return cohort_rows, source_rows, category_rows


def draw_panel_header(content: list[str], label: str, title: str, x: float, y: float) -> None:
    content.append(text(x, y, label, 28, "700"))
    content.append(text(x + 38, y, title, 19, "700"))


def draw_legend(content: list[str], items: list[tuple[str, str]], x: float, y: float) -> None:
    xx = x
    for label, color in items:
        content.append(rect(xx, y - 10, 21, 11, color, rx=2))
        content.append(text(xx + 29, y, label, 13, fill=COLORS["muted"]))
        xx += 185


def draw_flexible_legend(content: list[str], items: list[tuple[str, str]], x: float, y: float, gap: float = 185) -> None:
    xx = x
    for label, color in items:
        content.append(rect(xx, y - 10, 21, 11, color, rx=2))
        content.append(text(xx + 29, y, label, 13, fill=COLORS["muted"]))
        xx += gap


def axis_ticks(max_value: int, steps: int = 4) -> list[int]:
    if max_value <= 0:
        return [0]
    rough = max_value / steps
    magnitude = 10 ** (len(str(int(rough))) - 1)
    step = max(magnitude, round(rough / magnitude) * magnitude)
    top = ((max_value + step - 1) // step) * step
    return [value for value in range(0, top + step, step)]


def draw_axis(content: list[str], x: float, y: float, width: float, max_value: int) -> None:
    ticks = axis_ticks(max_value)
    top = max(ticks) or 1
    content.append(line(x, y, x + width, y, COLORS["grid"], 1))
    for tick in ticks:
        xx = x + width * tick / top
        content.append(line(xx, y - 5, xx, y + 4, COLORS["grid"], 1))
        content.append(text(xx, y + 22, f"{tick:,}", 11, fill=COLORS["muted"], anchor="middle"))


def draw_exact_hierarchy_bars(
    content: list[str],
    rows: list[dict[str, Any]],
    x: float,
    y: float,
    width: float,
    exact_key: str,
    hierarchy_key: str,
) -> None:
    label_w = 155
    chart_w = width - label_w - 105
    max_value = max(int(row[hierarchy_key]) for row in rows)
    top = max(axis_ticks(max_value)) or max_value
    draw_legend(content, [("Exact matching", COLORS["gray_light"]), ("Hierarchy-aware", COLORS["blue"])], x + label_w, y)
    yy = y + 50
    for row in rows:
        exact = int(row[exact_key])
        hierarchy = int(row[hierarchy_key])
        label = label_for_cohort(row["cohort"])
        content.append(text(x, yy + 24, label, 14, "700", COLORS["ink"]))
        exact_w = chart_w * exact / top
        hierarchy_w = chart_w * hierarchy / top
        content.append(rect(x + label_w, yy, exact_w, 17, COLORS["gray_light"], rx=2))
        content.append(rect(x + label_w, yy + 26, hierarchy_w, 17, COLORS["blue"], rx=2))
        content.append(text(x + label_w + exact_w + 8, yy + 14, f"{exact:,}", 12, fill=COLORS["muted"]))
        content.append(text(x + label_w + hierarchy_w + 8, yy + 40, f"{hierarchy:,}", 12))
        yy += 108
    draw_axis(content, x + label_w, yy - 30, chart_w, max_value)


def draw_normalization_stacked_bars(
    content: list[str],
    rows: list[dict[str, Any]],
    x: float,
    y: float,
    width: float,
) -> None:
    label_w = 155
    chart_w = width - label_w - 95
    max_value = max(int(row["variant_occurrences"]) for row in rows)
    top = max(axis_ticks(max_value)) or max_value
    draw_flexible_legend(content, NORMALIZATION_SERIES, x + label_w, y)
    yy = y + 50
    for row in rows:
        total = int(row["variant_occurrences"])
        normalized = int(row["parsed_variant_occurrences"])
        unnormalized = total - normalized
        label = label_for_cohort(row["cohort"])
        normalized_w = chart_w * normalized / top
        unnormalized_w = chart_w * unnormalized / top
        content.append(text(x, yy + 19, label, 14, "700", COLORS["ink"]))
        content.append(rect(x + label_w, yy, normalized_w, 21, COLORS["blue"], rx=2))
        content.append(rect(x + label_w + normalized_w, yy, unnormalized_w, 21, COLORS["gray_light"], rx=2))
        rate = normalized / total * 100 if total else 0
        content.append(text(x + label_w + normalized_w + unnormalized_w + 8, yy + 17, f"{normalized:,} ({rate:.1f}%)", 12))
        if unnormalized:
            content.append(text(x + label_w + normalized_w + max(unnormalized_w, 4) + 8, yy + 39, f"{unnormalized:,}", 11, fill=COLORS["muted"]))
        yy += 78
    draw_axis(content, x + label_w, yy - 18, chart_w, max_value)


def draw_grouped_added_bars(
    content: list[str],
    rows: list[tuple[str, int, int]],
    x: float,
    y: float,
    width: float,
    label_width: float,
    row_h: float,
    label_wrap: int,
) -> None:
    chart_w = width - label_width - 100
    max_value = max([max(left, right) for _, left, right in rows] + [1])
    top = max(axis_ticks(max_value)) or max_value
    draw_legend(content, [(label, color) for _, label, color in SERIES], x + label_width, y)
    yy = y + 48
    for label, left, right in rows:
        content.extend(multiline_text(x, yy + 12, wrap_label(label, label_wrap), 12, COLORS["muted"]))
        for offset, value, color in ((0, left, SERIES[0][2]), (22, right, SERIES[1][2])):
            if value:
                bar_w = max(3, chart_w * value / top)
                content.append(rect(x + label_width, yy + offset, bar_w, 15, color, rx=2))
                content.append(text(x + label_width + bar_w + 8, yy + offset + 13, f"{value:,}", 11))
        yy += row_h
    draw_axis(content, x + label_width, yy - 18, chart_w, max_value)


def write_figure(cohort_rows: list[dict[str, Any]], source_rows: list[dict[str, Any]], category_rows: list[dict[str, Any]]) -> None:
    content: list[str] = [rect(0, 0, 1700, 1050, COLORS["paper"])]

    draw_panel_header(content, "A", "Patient variant normalization", 55, 60)
    draw_normalization_stacked_bars(content, cohort_rows, 95, 105, 710)

    draw_panel_header(content, "B", "Retrieved evidence links", 885, 60)
    draw_exact_hierarchy_bars(
        content,
        cohort_rows,
        925,
        105,
        710,
        "source_level_exact_matches",
        "source_level_hierarchy_matches",
    )

    source_values = []
    for source in SOURCES:
        by_cohort = {row["cohort"]: int(row["hierarchy_added"]) for row in source_rows if row["evidence_source"] == source}
        source_values.append((source, by_cohort.get("GENIE_MSKCC_2022", 0), by_cohort.get("MSK_IMPACT_2017", 0)))

    draw_panel_header(content, "C", "Additional evidence links by source", 55, 560)
    draw_grouped_added_bars(content, source_values, 95, 605, 710, 145, 65, 24)

    category_values = []
    for category in sorted(
        {row["semantic_category"] for row in category_rows},
        key=lambda key: -sum(int(row["additional_evidence_links"]) for row in category_rows if row["semantic_category"] == key),
    ):
        by_cohort = {
            row["cohort"]: int(row["additional_evidence_links"])
            for row in category_rows
            if row["semantic_category"] == category
        }
        category_values.append((category, by_cohort.get("GENIE_MSKCC_2022", 0), by_cohort.get("MSK_IMPACT_2017", 0)))

    draw_panel_header(content, "D", "Additional evidence links by location and semantic category", 885, 560)
    draw_grouped_added_bars(content, category_values, 925, 605, 710, 285, 64, 28)

    svg = "\n".join(
        [
            '<svg xmlns="http://www.w3.org/2000/svg" width="170mm" height="105mm" viewBox="0 0 1700 1050">',
            *content,
            "</svg>",
        ]
    )
    svg_path = FIGURE_PATH.with_suffix(".svg")
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.write_text(svg, encoding="utf-8")

    import cairosvg  # type: ignore

    cairosvg.svg2pdf(url=str(svg_path), write_to=str(FIGURE_PATH))
    svg_path.unlink(missing_ok=True)


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    cohort_rows, source_rows, category_rows = collect_figure_data()
    write_figure(cohort_rows, source_rows, category_rows)
    print(f"Figure written to: {FIGURE_PATH}")
    print(f"Location-aware semantic category data written to: {LOCATION_CATEGORY_SUMMARY_PATH}")


if __name__ == "__main__":
    main()
