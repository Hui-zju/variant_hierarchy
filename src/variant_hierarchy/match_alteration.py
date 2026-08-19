from collections import deque
from typing import Any, Dict

from .normalize.ontology import VARIANT_ALTERATION_ONTOLOGY


def build_parents_map(ontology: Dict[str, Dict[str, Any]]) -> Dict[str, set[str]]:
    parents = {}
    for node, attrs in ontology.items():
        parents.setdefault(node, set())
        parent = attrs.get("parent_mechanism")
        if parent is not None:
            parents[node].add(parent)
            parents.setdefault(parent, set())
    return parents


def is_ancestor(ancestor: str, node: str, parents_map: dict[str, set[str]]) -> bool:
    if ancestor == node:
        return False

    q = deque([node])
    seen = set()
    while q:
        current = q.popleft()
        if current in seen:
            continue
        seen.add(current)
        for parent in parents_map.get(current, set()):
            if parent == ancestor:
                return True
            q.append(parent)
    return False


PARENTS_MAP = build_parents_map(VARIANT_ALTERATION_ONTOLOGY)


def alteration_match(alteration_a, alteration_b, parents_map=PARENTS_MAP):
    if len(alteration_a) > 0 and len(alteration_b) > 0:
        type_a = alteration_a["type"]
        type_b = alteration_b["type"]

        if alteration_a == alteration_b:
            return "equal"

        # Example: BRAF V600X contains BRAF V600E.
        if len(alteration_a) == 1 and type_a == type_b and "p_hgvs_posedit" in alteration_b:
            return "contains"

        if len(alteration_a) == 1 and type_a == type_b and "copyCount" in alteration_b:
            return "contains"

        if is_ancestor(type_a, type_b, parents_map):
            return "contains"

        return "unrelated"

    if len(alteration_a) == 0 and len(alteration_b) > 0:
        return "contains"
    if len(alteration_a) > 0 and len(alteration_b) == 0:
        return "unrelated"
    return "equal"


if __name__ == "__main__":
    print(is_ancestor("missense", "missense", PARENTS_MAP))
