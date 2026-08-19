from __future__ import annotations

import html
from typing import Any

from evaluate.common import write_tsv
from evaluate.hierarchy_construction.config import (
    COLORS,
    RELATION_SOURCE_GROUPS,
    RESULT_DIR,
    SOURCES,
    manuscript_class_label,
    source_set_sort_key,
)


def svg_tag(name: str, attrs: dict[str, Any], content: str = "") -> str:
    attr_s = " ".join(f'{key}="{html.escape(str(value), quote=True)}"' for key, value in attrs.items())
    if content:
        return f"<{name} {attr_s}>{content}</{name}>"
    return f"<{name} {attr_s}/>"


def text(x: float, y: float, value: Any, size: int = 13, weight: str = "400", fill: str = COLORS["ink"], anchor: str = "start") -> str:
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
        html.escape(str(value)),
    )


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str | None = None, rx: int = 0) -> str:
    attrs = {
        "x": f"{x:.1f}",
        "y": f"{y:.1f}",
        "width": f"{max(w, 0):.1f}",
        "height": f"{max(h, 0):.1f}",
        "fill": fill,
    }
    if stroke:
        attrs["stroke"] = stroke
    if rx:
        attrs["rx"] = rx
    return svg_tag("rect", attrs)


def circle(cx: float, cy: float, r: float, fill: str, stroke: str | None = None) -> str:
    attrs = {"cx": f"{cx:.1f}", "cy": f"{cy:.1f}", "r": f"{r:.1f}", "fill": fill}
    if stroke:
        attrs["stroke"] = stroke
    return svg_tag("circle", attrs)


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


def svg_document(width: int, height: int, content: list[str]) -> str:
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            rect(0, 0, width, height, COLORS["paper"]),
            *content,
            "</svg>",
            "",
        ]
    )


def draw_horizontal_bars(
    rows: list[tuple[str, int]],
    x: int,
    y: int,
    width: int,
    bar_h: int,
    gap: int,
    color: str,
    max_value: int | None = None,
    label_width: int = 210,
) -> list[str]:
    out: list[str] = []
    max_v = max_value or max([value for _, value in rows] + [1])
    chart_x = x + label_width
    for i, (label, value) in enumerate(rows):
        yy = y + i * (bar_h + gap)
        bar_w = width * value / max_v
        out.append(text(x, yy + bar_h * 0.72, label, size=12, fill=COLORS["muted"]))
        out.append(rect(chart_x, yy, bar_w, bar_h, color, rx=2))
        out.append(text(chart_x + bar_w + 8, yy + bar_h * 0.72, f"{value:,}", size=12, fill=COLORS["ink"]))
    return out


def write_figure2(stats: dict[str, Any]) -> None:
    variant_rows = [
        (manuscript_class_label(key), value)
        for key, value in stats["variant_class_counts"].most_common()
    ]
    source_rows = [(source, stats["source_counts"].get(source, 0)) for source in SOURCES]
    intersections = [
        (key, value)
        for key, value in stats["source_set_counts"].items()
        if key
    ]
    intersections.sort(key=source_set_sort_key)

    content: list[str] = []

    content.append(text(32, 42, "A", 22, "700"))
    content.append(text(62, 42, "Source-level contribution", 16, "700"))
    content.extend(
        draw_horizontal_bars(
            source_rows,
            x=62,
            y=70,
            width=390,
            bar_h=25,
            gap=18,
            color=COLORS["green"],
            label_width=120,
        )
    )

    content.append(text(32, 278, "B", 22, "700"))
    content.append(text(62, 278, "Semantic category distribution", 16, "700"))
    content.append(text(62, 318, "Total normalized concepts", 12, fill=COLORS["muted"]))
    content.append(text(452, 318, f"{len(stats['variants']):,}", 16, "700", COLORS["ink"], anchor="end"))
    content.append(line(62, 331, 452, 331, COLORS["grid"], 1))
    content.extend(
        draw_horizontal_bars(
            variant_rows,
            x=62,
            y=342,
            width=260,
            bar_h=20,
            gap=11,
            color=COLORS["blue"],
            label_width=250,
        )
    )

    content.append(text(610, 42, "C", 22, "700"))
    content.append(text(640, 42, "Exact source-set intersections", 16, "700"))
    matrix_x = 632
    matrix_gap = 52
    bar_x = 835
    y0 = 92
    row_h = 31
    max_intersection = max([value for _, value in intersections] + [1])
    for idx, source in enumerate(SOURCES):
        content.append(text(matrix_x + idx * matrix_gap, y0 - 18, source, 11, fill=COLORS["muted"], anchor="middle"))
    content.append(circle(matrix_x, y0 + len(intersections) * row_h + 14, 5.4, COLORS["ink"]))
    content.append(text(matrix_x + 14, y0 + len(intersections) * row_h + 18, "included", 11, fill=COLORS["muted"]))
    content.append(circle(matrix_x + 86, y0 + len(intersections) * row_h + 14, 5.4, COLORS["gray_light"]))
    content.append(text(matrix_x + 100, y0 + len(intersections) * row_h + 18, "not included", 11, fill=COLORS["muted"]))
    for i, (source_set, value) in enumerate(intersections):
        yy = y0 + i * row_h
        active_positions = []
        for idx, source in enumerate(SOURCES):
            active = source in source_set
            cx = matrix_x + idx * matrix_gap
            active_positions.append(cx) if active else None
            content.append(circle(cx, yy, 5.4, COLORS["ink"] if active else COLORS["gray_light"]))
        if len(active_positions) > 1:
            content.append(line(active_positions[0], yy, active_positions[-1], yy, COLORS["ink"], 1.5))
        bar_w = 315 * value / max_intersection
        content.append(rect(bar_x, yy - 10, bar_w, 18, COLORS["orange"], rx=2))
        content.append(text(bar_x + bar_w + 8, yy - 10 + 14, f"{value:,}", 11))

    figure_path = RESULT_DIR / "figure2_variant_landscape.svg"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    figure_path.write_text(
        svg_document(1200, 580, content),
        encoding="utf-8",
    )


def write_figure3(stats: dict[str, Any]) -> None:
    full_edges = len(stats["full_edges"])
    non_edges = len(stats["non_edges"])
    removed_edges = full_edges - non_edges

    edge_rows = []
    for key in sorted(
        set(stats["full_edge_categories"]) | set(stats["non_edge_categories"]),
        key=lambda item: stats["non_edge_categories"].get(item, 0),
        reverse=True,
    ):
        edge_rows.append(
            (
                f"{manuscript_class_label(key[0])} -> {manuscript_class_label(key[1])}",
                stats["full_edge_categories"].get(key, 0),
                stats["non_edge_categories"].get(key, 0),
            )
        )

    content: list[str] = []

    content.append(text(32, 42, "A", 22, "700"))
    content.append(text(62, 42, "Graph-scale summary", 16, "700"))
    box_y = 74
    boxes = [
        ("Normalized variant concepts", f"{len(stats['variants']):,}"),
        ("Full edges", f"{full_edges:,}"),
        ("Non-redundant edges", f"{non_edges:,}"),
        ("Removed edges", f"{removed_edges:,} ({removed_edges / full_edges * 100:.1f}%)"),
    ]
    for i, (label, value) in enumerate(boxes):
        x = 62 + i * 260
        content.append(rect(x, box_y, 220, 72, "#f9fafb", COLORS["grid"], rx=4))
        content.append(text(x + 18, box_y + 28, label, 12, fill=COLORS["muted"]))
        content.append(text(x + 18, box_y + 56, value, 22, "700", COLORS["ink"]))
        if i < len(boxes) - 1:
            content.append(line(x + 228, box_y + 36, x + 252, box_y + 36, COLORS["muted"], 1.5))
            content.append(text(x + 240, box_y + 31, ">", 14, "700", COLORS["muted"], "middle"))

    content.append(text(32, 202, "B", 22, "700"))
    content.append(text(62, 202, "Relation sources", 16, "700"))
    source_label_x = 62
    source_bar_x = 430
    source_y = 234
    source_row_h = 56
    source_bar_w = 520
    max_source = max(
        [
            stats["full_relation_source_groups"].get(key, 0)
            for key, _, _ in RELATION_SOURCE_GROUPS
        ]
        + [
            stats["non_relation_source_groups"].get(key, 0)
            for key, _, _ in RELATION_SOURCE_GROUPS
        ]
        + [1]
    )
    for i, (source_key, label, _) in enumerate(RELATION_SOURCE_GROUPS):
        yy = source_y + i * source_row_h
        full_value = stats["full_relation_source_groups"].get(source_key, 0)
        non_value = stats["non_relation_source_groups"].get(source_key, 0)
        full_w = source_bar_w * full_value / max_source
        non_w = source_bar_w * non_value / max_source
        content.append(text(source_label_x, yy + 24, label, 12, fill=COLORS["muted"]))
        content.append(rect(source_bar_x, yy, full_w, 14, COLORS["gray_light"], rx=2))
        content.append(rect(source_bar_x, yy + 19, non_w, 14, COLORS["green"], rx=2))
        content.append(text(source_bar_x + full_w + 8, yy + 12, f"{full_value:,}", 11))
        content.append(text(source_bar_x + non_w + 8, yy + 31, f"{non_value:,}", 11))
    source_legend_y = source_y + len(RELATION_SOURCE_GROUPS) * source_row_h + 8
    content.append(rect(source_bar_x, source_legend_y, 18, 9, COLORS["gray_light"], rx=2))
    content.append(text(source_bar_x + 24, source_legend_y + 9, "Full", 11, fill=COLORS["muted"]))
    content.append(rect(source_bar_x + 82, source_legend_y, 18, 9, COLORS["green"], rx=2))
    content.append(text(source_bar_x + 106, source_legend_y + 9, "Non-redundant", 11, fill=COLORS["muted"]))

    content.append(text(32, 410, "C", 22, "700"))
    content.append(text(62, 410, "Semantic edge categories before and after redundancy removal", 16, "700"))
    max_edge = max([full for _, full, _ in edge_rows] + [1])
    label_x = 62
    chart_x = 515
    y0 = 438
    row_h = 29
    for i, (label, full_count, non_count) in enumerate(edge_rows):
        yy = y0 + i * row_h
        full_w = 420 * full_count / max_edge
        non_w = 420 * non_count / max_edge
        content.append(text(label_x, yy + 15, label, 10.5, fill=COLORS["muted"]))
        content.append(rect(chart_x, yy, full_w, 10, COLORS["gray_light"], rx=2))
        content.append(rect(chart_x, yy + 12, non_w, 10, COLORS["blue"], rx=2))
        content.append(text(chart_x + max(full_w, non_w) + 8, yy + 17, f"{full_count:,} / {non_count:,}", 10.5))
    legend_y = y0 + len(edge_rows) * row_h + 14
    content.append(rect(chart_x, legend_y, 18, 9, COLORS["gray_light"], rx=2))
    content.append(text(chart_x + 24, legend_y + 9, "Full", 11, fill=COLORS["muted"]))
    content.append(rect(chart_x + 82, legend_y, 18, 9, COLORS["blue"], rx=2))
    content.append(text(chart_x + 106, legend_y + 9, "Non-redundant", 11, fill=COLORS["muted"]))

    figure_path = RESULT_DIR / "figure3_hierarchy_construction.svg"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    figure_path.write_text(
        svg_document(1200, 925, content),
        encoding="utf-8",
    )


def prepare_submission_files() -> None:
    figure_paths = [
        RESULT_DIR / "figure2_variant_landscape.svg",
        RESULT_DIR / "figure3_hierarchy_construction.svg",
    ]
    notes = []
    try:
        import cairosvg  # type: ignore
    except ModuleNotFoundError:
        write_tsv(
            RESULT_DIR / "figure_submission_notes.tsv",
            [
                {
                    "item": "PDF/PNG export",
                    "status": "not generated",
                    "note": (
                        "Install cairosvg to export accepted submission formats. "
                        "Current SVG files are vector source files; the journal's listed "
                        "submission formats include PDF and PNG, but not SVG."
                    ),
                },
                {
                    "item": "Graphic title placement",
                    "status": "ok",
                    "note": "Figure titles and legends are not embedded in the graphic files.",
                },
                {
                    "item": "Composite panels",
                    "status": "ok",
                    "note": "Each multi-panel figure is generated as a single composite file.",
                },
            ],
            ["item", "status", "note"],
        )
        return

    for svg_path in figure_paths:
        pdf_path = svg_path.with_suffix(".pdf")
        cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
        svg_path.unlink(missing_ok=True)
        svg_path.with_suffix(".png").unlink(missing_ok=True)
        notes.append(
            {
                "item": svg_path.name,
                "status": "generated",
                "note": f"Exported accepted submission file: {pdf_path.name}",
            }
        )

    notes.extend(
        [
            {
                "item": "Graphic title placement",
                "status": "ok",
                "note": "Figure titles and legends are not embedded in the graphic files.",
            },
            {
                "item": "Composite panels",
                "status": "ok",
                "note": "Each multi-panel figure is generated as a single composite file.",
            },
        ]
    )
    write_tsv(RESULT_DIR / "figure_submission_notes.tsv", notes, ["item", "status", "note"])


def write_figures(stats: dict[str, Any]) -> None:
    write_figure2(stats)
    write_figure3(stats)
    prepare_submission_files()


def main() -> None:
    from evaluate.hierarchy_construction.run import collect_stats

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    stats = collect_stats()
    write_figures(stats)
    print(f"Results written to: {RESULT_DIR}")


if __name__ == "__main__":
    main()
