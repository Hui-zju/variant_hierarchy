from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def add_project_paths() -> None:
    for path in (PROJECT_ROOT / "src", PROJECT_ROOT, PROJECT_ROOT / "hierarchy"):
        path_s = str(path)
        if path_s not in sys.path:
            sys.path.insert(0, path_s)


def read_tsv(path: Path, skip_comments: bool = False) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        lines = (line for line in f if not (skip_comments and line.startswith("#")))
        return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
