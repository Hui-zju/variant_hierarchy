import json
import re
from collections import defaultdict, deque
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Literal, Optional

from .normalize.gene_norm import gene_norm
from .normalize.variant_parser import parse_variant


DEFAULT_HIERARCHY_PATH = Path(__file__).resolve().parent / "data" / "hierarchy_non_redundant.json"


def _as_list(variants: str | Iterable[str]) -> list[str]:
    return [variants] if isinstance(variants, str) else list(variants)


class VariantHierarchy:
    def __init__(self, hierarchy_path: str | Path = DEFAULT_HIERARCHY_PATH):
        self.hierarchy_path = Path(hierarchy_path)
        with self.hierarchy_path.open("r", encoding="utf-8") as f:
            self.relations = json.load(f)

        self.children_by_parent: dict[str, list[str]] = defaultdict(list)
        self.parents_by_child: dict[str, list[str]] = defaultdict(list)
        self.nodes: set[str] = set()

        for relation in self.relations:
            parent = relation["parent"]
            child = relation["child"]
            self.children_by_parent[parent].append(child)
            self.parents_by_child[child].append(parent)
            self.nodes.add(parent)
            self.nodes.add(child)

    @staticmethod
    def _gene_from_name(name: str) -> str:
        return name.split(" ", 1)[0] if name else ""

    @staticmethod
    def _tree_node(name: str, children: Optional[list[dict]] = None) -> dict:
        node = {"name": name}
        if children is not None:
            node["children"] = children
        return node

    def normalize(self, variant: str, reference: Optional[str] = None) -> Optional[str]:
        if variant in self.nodes:
            return variant

        parsed = parse_variant(variant, reference=reference)
        if not parsed:
            return variant if variant in self.nodes else None

        display_name = parsed["display_name"]
        return display_name if display_name in self.nodes else display_name

    def normalize_gene(self, gene: str) -> Optional[str]:
        return gene_norm(gene)

    def find_gene_root(self, gene: str) -> Optional[str]:
        normalized_gene = self.normalize_gene(gene)
        if not normalized_gene:
            return None

        exact_root = f"{normalized_gene} mutation"
        if exact_root in self.nodes:
            return exact_root

        candidates = [
            node for node in self.nodes
            if self._gene_from_name(node).upper() == normalized_gene
            and node in self.children_by_parent
        ]
        if not candidates:
            return None
        return sorted(candidates, key=lambda node: (node != exact_root, len(node), node))[0]

    def parents(self, variant: str, reference: Optional[str] = None) -> list[str]:
        normalized = self.normalize(variant, reference=reference)
        if not normalized:
            return []
        return list(self.parents_by_child.get(normalized, []))

    def children(self, variant: str, reference: Optional[str] = None) -> list[str]:
        normalized = self.normalize(variant, reference=reference)
        if not normalized:
            return []
        return list(self.children_by_parent.get(normalized, []))

    def ancestors(self, variant: str, reference: Optional[str] = None) -> list[str]:
        normalized = self.normalize(variant, reference=reference)
        if not normalized:
            return []

        ancestors = []
        seen = set()
        queue = deque(self.parents_by_child.get(normalized, []))
        while queue:
            node = queue.popleft()
            if node in seen:
                continue
            seen.add(node)
            ancestors.append(node)
            queue.extend(self.parents_by_child.get(node, []))
        return ancestors

    def descendants(self, variant: str, reference: Optional[str] = None) -> list[str]:
        normalized = self.normalize(variant, reference=reference)
        if not normalized:
            return []

        descendants = []
        seen = set()
        queue = deque(self.children_by_parent.get(normalized, []))
        while queue:
            node = queue.popleft()
            if node in seen:
                continue
            seen.add(node)
            descendants.append(node)
            queue.extend(self.children_by_parent.get(node, []))
        return descendants

    def siblings(self, variant: str, reference: Optional[str] = None) -> list[str]:
        normalized = self.normalize(variant, reference=reference)
        if not normalized:
            return []

        siblings = []
        for parent in self.parents_by_child.get(normalized, []):
            siblings.extend(self.children_by_parent.get(parent, []))
        return [node for node in dict.fromkeys(siblings) if node != normalized]

    def relation(self, parent: str, child: str) -> bool:
        parent_name = self.normalize(parent)
        child_name = self.normalize(child)
        return bool(parent_name and child_name and child_name in self.children_by_parent.get(parent_name, []))

    def descendant_tree(self, root: str, reference: Optional[str] = None) -> Optional[dict]:
        root_name = self.normalize(root, reference=reference)
        if not root_name:
            root_name = self.find_gene_root(root)
        if not root_name:
            return None
        return self._build_descendant_tree(root_name, set())

    def _build_descendant_tree(self, root: str, seen: set[str]) -> dict:
        if root in seen:
            return self._tree_node(root, [])

        next_seen = set(seen)
        next_seen.add(root)
        children = [
            self._build_descendant_tree(child, next_seen)
            for child in self.children_by_parent.get(root, [])
        ]
        return self._tree_node(root, children)

    def gene_hierarchy(self, gene: str) -> Optional[dict]:
        root = self.find_gene_root(gene)
        if not root:
            return None

        normalized_gene = self.normalize_gene(gene)
        return {
            "type": "gene",
            "gene": normalized_gene,
            "root": root,
            "tree": self._build_descendant_tree(root, set()),
        }

    def ancestor_paths(self, variant: str, reference: Optional[str] = None) -> list[list[str]]:
        normalized = self.normalize(variant, reference=reference)
        if not normalized:
            return []

        paths = []
        stack = [(normalized, [normalized], {normalized})]
        while stack:
            current, path, seen = stack.pop()
            parents = self.parents_by_child.get(current, [])
            if not parents:
                paths.append(list(reversed(path)))
                continue

            for parent in parents:
                if parent not in seen:
                    stack.append((parent, path + [parent], seen | {parent}))

        return sorted(paths, key=lambda item: (len(item), item))

    def variant_hierarchy(self, variant: str, reference: Optional[str] = None) -> Optional[dict]:
        normalized = self.normalize(variant, reference=reference)
        if not normalized:
            return None

        paths = self.ancestor_paths(normalized)
        path_trees = [self._build_path_tree(path) for path in paths] or [
            self._build_descendant_tree(normalized, set())
        ]
        return {
            "type": "variant",
            "variant": normalized,
            "ancestor_paths": paths,
            "trees": path_trees,
            "descendant_tree": self._build_descendant_tree(normalized, set()),
        }

    def _build_path_tree(self, path: list[str]) -> dict:
        if not path:
            return self._tree_node("")

        root = self._tree_node(path[-1])
        if len(path) == 1:
            return self._build_descendant_tree(path[0], set())

        child = self._build_descendant_tree(path[-1], set())
        for node in reversed(path[:-1]):
            child = self._tree_node(node, [child])
        return child

    def hierarchy(
        self,
        query: str,
        reference: Optional[str] = None,
        mode: Literal["auto", "gene", "variant"] = "auto",
    ) -> Optional[dict]:
        if mode == "gene":
            return self.gene_hierarchy(query)
        if mode == "variant":
            return self.variant_hierarchy(query, reference=reference)
        if mode != "auto":
            raise ValueError("mode must be one of: 'auto', 'gene', 'variant'")

        if not re.search(r"\s", query.strip()) and self.find_gene_root(query):
            return self.gene_hierarchy(query)
        return self.variant_hierarchy(query, reference=reference)


@lru_cache(maxsize=1)
def get_default_hierarchy() -> VariantHierarchy:
    return VariantHierarchy()
