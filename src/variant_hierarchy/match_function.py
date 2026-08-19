from collections import deque
from typing import Any, Dict

from .normalize.ontology import VARIANT_FUNCTION_ONTOLOGY


def build_function_parents_map(ontology: Dict[str, Dict[str, Any]]) -> Dict[str, set[str]]:
    parents = {}
    for node, attrs in ontology.items():
        functional_value = attrs.get("functional_value", node)
        parents.setdefault(functional_value, set())
        parent = attrs.get("parent_functional")
        if parent is not None:
            parents[functional_value].add(parent)
            parents.setdefault(parent, set())
    return parents


def is_ancestor(ancestor: str, node: str, parents_map: dict[str, set[str]]) -> bool:
    if ancestor == node:
        return True
    q = deque([node])
    seen = set()
    while q:
        x = q.popleft()
        if x in seen:
            continue
        seen.add(x)
        for p in parents_map.get(x, set()):
            if p == ancestor:
                return True
            q.append(p)
    return False


FUNCTION_PARENTS_MAP = build_function_parents_map(VARIANT_FUNCTION_ONTOLOGY)


def function_match(func_a, func_b, parents_map=FUNCTION_PARENTS_MAP):
    if func_a == func_b:
        return "equal"

    if len(func_a) == 0 and len(func_b) > 0:
        return "contains"
    if len(func_a) > 0 and len(func_b) == 0:
        return "unrelated"

    has_stricter_child_value = False
    for domain, value_a in func_a.items():
        if domain not in func_b:
            return "unrelated"

        value_b = func_b[domain]
        if value_a == value_b:
            continue

        if is_ancestor(value_a, value_b, parents_map):
            has_stricter_child_value = True
            continue

        return "unrelated"

    if has_stricter_child_value or len(func_a) < len(func_b):
        return "contains"
    return "equal"
