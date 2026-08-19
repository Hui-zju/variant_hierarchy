import re
from typing import Dict, Any, Optional
from .ontology import VARIANT_ALTERATION_ONTOLOGY
from .ontology import VARIANT_FUNCTION_ONTOLOGY


DELETION_KEYS = frozenset({"deletion", "inframe_deletion", "deletion_cna", "exon_deletion_cna"})
INSERTION_KEYS = frozenset({"insertion", "inframe_insertion"})
EXPLICIT_INFRAME_DELETION_RE = re.compile(
    r"\b(?:in[\s\-]?frame|non[\s\-]?frameshift)\s+deletions?\b",
    flags=re.IGNORECASE,
)
COPY_NUMBER_DELETION_RE = re.compile(
    r"\b(?:homozygous|gene|copy\s+number|deep)\s+(?:deletion|del)\b",
    flags=re.IGNORECASE,
)


EXON_SHORTHAND_ALTERATIONS = {
    ("EGFR", "19", "deletion"): "inframe_deletion",
    ("MET", "14", "deletion"): "inframe_deletion",
    ("KIT", "11", "deletion"): "inframe_deletion",
    ("PDGFRA", "18", "deletion"): "inframe_deletion",
    ("EGFR", "19", "insertion"): "inframe_insertion",
    ("EGFR", "20", "insertion"): "inframe_insertion",
    ("ERBB2", "20", "insertion"): "inframe_insertion",
    ("PDGFRA", "18", "insertion"): "inframe_insertion",
}


def _exon_shorthand_alteration(gene: Optional[str], exon: Optional[str], mechanism: str) -> Optional[str]:
    if not gene or not exon:
        return None
    return EXON_SHORTHAND_ALTERATIONS.get((gene.upper(), exon, mechanism))


def _copy_number_deletion_key(keys: set[str], location: Optional[str]) -> Optional[str]:
    if location == "Exon" and "exon_deletion_cna" in keys:
        return "exon_deletion_cna"
    if "deletion_cna" in keys:
        return "deletion_cna"
    return None


def _resolve_indel_candidate(
    text: str,
    candidates: list[Dict[str, Any]],
    location: Optional[str],
    variant_class: Optional[str],
    gene: Optional[str],
    exon: Optional[str],
) -> Optional[str]:
    keys = {candidate["key"] for candidate in candidates}
    has_deletion = bool(keys & DELETION_KEYS)
    has_insertion = bool(keys & INSERTION_KEYS)
    if not has_deletion and not has_insertion:
        return None

    if has_deletion:
        if EXPLICIT_INFRAME_DELETION_RE.search(text or "") and "inframe_deletion" in keys:
            return "inframe_deletion"

        if variant_class == "Copy Number Variant":
            return _copy_number_deletion_key(keys, location)

        if variant_class in (None, "Sequence Variant") and location == "Exon":
            shorthand = _exon_shorthand_alteration(gene, exon, "deletion")
            if shorthand:
                return shorthand

        if variant_class == "Sequence Variant" and "deletion" in keys:
            return "deletion"

        if COPY_NUMBER_DELETION_RE.search(text or ""):
            cna_key = _copy_number_deletion_key(keys, None)
            if cna_key:
                return cna_key

        if location == "Amino Acid" and "inframe_deletion" in keys:
            return "inframe_deletion"
        if location == "Gene":
            cna_key = _copy_number_deletion_key(keys, location)
            if cna_key:
                return cna_key
        if "deletion" in keys:
            return "deletion"

    if has_insertion and variant_class in (None, "Sequence Variant") and location == "Exon":
        return _exon_shorthand_alteration(gene, exon, "insertion")

    return None


def _find_best_key(
    text: str,
    ontology: Dict[str, Dict[str, Any]],
    location: Optional[str] = None,
    variant_class: Optional[str] = None,
    gene: Optional[str] = None,
    exon: Optional[str] = None,
) -> Optional[str]:
    candidates = []
    if not text:
        return None

    for key, item in ontology.items():
        pattern = item.get("regex_pattern")
        if not pattern:
            continue

        if re.search(pattern, text, flags=re.IGNORECASE):
            candidates.append({
                "key": key,
                "priority": item.get("priority", 0),
                "variant_class": item.get("variant_class"),
                "minimal_location_type": item.get("minimal_location_type"),
            })

    if not candidates:
        return None

    indel_candidate = _resolve_indel_candidate(text, candidates, location, variant_class, gene, exon)
    if indel_candidate:
        return indel_candidate

    if variant_class:
        class_candidates = [
            candidate for candidate in candidates
            if candidate.get("variant_class") == variant_class
        ]
        if class_candidates:
            candidates = class_candidates

    if location:
        location_candidates = [
            candidate for candidate in candidates
            if candidate.get("minimal_location_type") == location
        ]
        if location_candidates:
            candidates = location_candidates

    candidates.sort(key=lambda x: x["priority"], reverse=True)
    if len(candidates) > 0:
        return candidates[0]["key"]
    return None


def get_normalized_mechanism(
    text: str,
    location: Optional[str] = None,
    variant_class: Optional[str] = None,
    gene: Optional[str] = None,
    exon: Optional[str] = None,
) -> Optional[str]:
    return _find_best_key(text, VARIANT_ALTERATION_ONTOLOGY, location, variant_class, gene, exon)

def get_normalized_function(text: str) -> Optional[str]:
    return _find_best_key(text, VARIANT_FUNCTION_ONTOLOGY)


if __name__ == "__main__":
    # examples = [
    #     "activating missense mutation",
    #     "MET exon 14 skipping",
    #     "EML4-ALK fusion oncogenic",
    #     "homozygous deletion",
    #     "splice donor mutation",
    #     "internal tandem duplication",
    #     "wild type",
    #     "ERBB2 overexpression",
    #     "likely oncogenic fusion",
    # ]
    #
    # for e in examples:
    #     print(f"Input: {e}")
    #     print("-" * 8)
    #     print(get_normalized_mechanism(e))
    #     print(get_normalized_function(e))
    #     print("-" * 80)

    res = get_normalized_mechanism("oncogenic mutation")
    print('ok')
