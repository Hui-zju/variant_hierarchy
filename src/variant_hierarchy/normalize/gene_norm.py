import csv
import re
from functools import lru_cache
from pathlib import Path
from typing import Optional


ALIAS_PATH = Path(__file__).resolve().parent.parent / "data" / "gene_aliases.tsv"
DEFAULT_VICC_BASE_URL = "http://localhost:8001"
DEFAULT_VICC_TIMEOUT = 1.0
MIN_VICC_MATCH_TYPE = 100


def _compact_key(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def _exact_key(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip()).casefold()


def _add_index_value(index, key: str, symbol: str):
    if not key:
        return
    if key not in index:
        index[key] = symbol
        return
    existing = index[key]
    if existing is None:
        return
    if existing and existing != symbol:
        index[key] = None
        return


@lru_cache(maxsize=1)
def _load_alias_indices():
    exact_index = {}
    compact_index = {}
    symbols = set()

    with ALIAS_PATH.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            symbol = (row.get("symbol") or "").strip().upper()
            if not symbol:
                continue
            symbols.add(symbol)
            aliases = [symbol]
            aliases.extend(x.strip() for x in (row.get("aliases") or "").split("|") if x.strip())
            for alias in aliases:
                _add_index_value(exact_index, _exact_key(alias), symbol)
                _add_index_value(compact_index, _compact_key(alias), symbol)

    return symbols, exact_index, compact_index


def gene_norm_local(q) -> Optional[str]:
    if q is None:
        return None

    text = str(q).strip()
    if not text:
        return None

    symbols, exact_index, compact_index = _load_alias_indices()
    upper_text = text.upper()
    if upper_text in symbols:
        return upper_text

    exact_match = exact_index.get(_exact_key(text))
    if exact_match:
        return exact_match

    compact_match = compact_index.get(_compact_key(text))
    if compact_match:
        return compact_match

    acronym_match = re.search(r"[A-Z0-9]{2,}", text)
    if acronym_match:
        acronym = acronym_match.group()
        if acronym in symbols:
            return acronym
        exact_match = exact_index.get(_exact_key(acronym))
        if exact_match:
            return exact_match
        compact_match = compact_index.get(_compact_key(acronym))
        if compact_match:
            return compact_match

    return None


def gene_norm_vicc(q, timeout: float = DEFAULT_VICC_TIMEOUT) -> Optional[str]:
    if q is None:
        return None

    try:
        import requests

        resp = requests.get(
            f"{DEFAULT_VICC_BASE_URL}/gene/normalize",
            params={"q": q},
            timeout=timeout,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception:
        return None

    if data.get("match_type", 0) < MIN_VICC_MATCH_TYPE:
        return None

    gene = data.get("gene")
    if not gene:
        return None
    return gene.get("name")


@lru_cache(maxsize=10000)
def _gene_norm_cached(q) -> Optional[str]:
    if q is None:
        return None

    return gene_norm_vicc(q) or gene_norm_local(q)


def gene_norm(q) -> Optional[str]:
    """
    Normalize a gene name or alias to a canonical gene symbol.

    First try a VICC gene-normalizer service at DEFAULT_VICC_BASE_URL.
    If unavailable or unresolved, fall back to the packaged alias table.
    """
    return _gene_norm_cached(q)


def _clear_gene_norm_cache():
    _gene_norm_cached.cache_clear()
    _load_alias_indices.cache_clear()


gene_norm.cache_clear = _clear_gene_norm_cache
gene_norm.cache_info = _gene_norm_cached.cache_info


if __name__ == "__main__":
    print(gene_norm("EGFR"))
