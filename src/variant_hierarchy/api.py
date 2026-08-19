from typing import Optional

from .hierarchy import get_default_hierarchy
from .normalize.variant_parser import parse_variant


def normalize_variant(variant: str, reference: Optional[str] = None, **parse_kwargs) -> Optional[dict]:
    return parse_variant(variant, reference=reference, **parse_kwargs)


def parents(variant: str, reference: Optional[str] = None) -> list[str]:
    return get_default_hierarchy().parents(variant, reference=reference)


def children(variant: str, reference: Optional[str] = None) -> list[str]:
    return get_default_hierarchy().children(variant, reference=reference)


def ancestors(variant: str, reference: Optional[str] = None) -> list[str]:
    return get_default_hierarchy().ancestors(variant, reference=reference)


def descendants(variant: str, reference: Optional[str] = None) -> list[str]:
    return get_default_hierarchy().descendants(variant, reference=reference)


def gene_hierarchy(gene: str) -> Optional[dict]:
    return get_default_hierarchy().gene_hierarchy(gene)


def variant_hierarchy(variant: str, reference: Optional[str] = None) -> Optional[dict]:
    return get_default_hierarchy().variant_hierarchy(variant, reference=reference)


def hierarchy(query: str, reference: Optional[str] = None, mode: str = "auto") -> Optional[dict]:
    return get_default_hierarchy().hierarchy(query, reference=reference, mode=mode)
