import csv
import json
from collections import defaultdict, deque
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
HIERARCHY_DIR = DATA_DIR / "hierarchy"
GRAPH_DIR = DATA_DIR / "graph"
DEFAULT_HIERARCHY_PATH = HIERARCHY_DIR / "hierarchy.json"
NON_REDUNDANT_HIERARCHY_PATH = HIERARCHY_DIR / "hierarchy_non_redundant.json"
NODE_PATH = GRAPH_DIR / "node.csv"
EDGE_PATH = GRAPH_DIR / "edge.csv"
NODE_WITH_INFO_PATH = GRAPH_DIR / "node_with_info.csv"
EDGE_WITH_INFO_PATH = GRAPH_DIR / "edge_with_info.csv"


class VariationHierarchicalRelation:
    """Query, reduce, print, and export variant hierarchy relations."""

    def __init__(self, relations=None, reduce=True):
        self.row_relations = relations if relations is not None else []
        self.variant_dict = self._collect_variant_info(self.row_relations)
        self.all_relations = self._deduplicate_relations(
            [self._normalize_relation(relation) for relation in self.row_relations]
        )
        self.h_relations = (
            self.get_non_redundant_relations(self.all_relations)
            if reduce
            else list(self.all_relations)
        )
        self._rebuild_indexes()

    @staticmethod
    def _collect_variant_info(relations):
        variant_dict = {}
        for relation in relations:
            for role in ("parent", "child"):
                variant = relation[role]
                if isinstance(variant, dict):
                    display_name = variant.get("display_name") or variant.get("displayName")
                    if display_name:
                        variant_dict[display_name] = variant
        return variant_dict

    def _normalize_relation(self, relation):
        parent_variant = relation["parent"]
        child_variant = relation["child"]
        parent = (
            parent_variant.get("display_name") or parent_variant.get("displayName")
            if isinstance(parent_variant, dict)
            else parent_variant
        )
        child = (
            child_variant.get("display_name") or child_variant.get("displayName")
            if isinstance(child_variant, dict)
            else child_variant
        )
        if not parent or not child:
            raise ValueError(f"Invalid hierarchy relation: {relation}")

        normalized = {
            "parent": parent,
            "child": child,
            "relation": relation.get("relation", "contain"),
        }

        if "source" in relation:
            normalized["source"] = relation["source"]
        else:
            normalized["source"] = {
                "parent": parent_variant.get("source", []) if isinstance(parent_variant, dict) else [],
                "child": child_variant.get("source", []) if isinstance(child_variant, dict) else [],
            }

        return normalized

    @staticmethod
    def _deduplicate_relations(relations):
        deduped = {}
        for relation in relations:
            key = (relation["parent"], relation["child"], relation["relation"])
            if key not in deduped:
                deduped[key] = relation
                continue

            old_sources = deduped[key].get("source", [])
            new_sources = relation.get("source", [])
            if isinstance(old_sources, list) and isinstance(new_sources, list):
                deduped[key]["source"] = sorted(set(old_sources) | set(new_sources))
        return list(deduped.values())

    def _rebuild_indexes(self):
        self.children_by_parent = defaultdict(list)
        self.parents_by_child = defaultdict(list)
        self.nodes = set()

        for relation in self.h_relations:
            parent = relation["parent"]
            child = relation["child"]
            self.children_by_parent[parent].append(child)
            self.parents_by_child[child].append(parent)
            self.nodes.add(parent)
            self.nodes.add(child)

    @staticmethod
    def _as_list(variations):
        return [variations] if isinstance(variations, str) else list(variations)

    def get_h_relations_dict(self):
        return {parent: list(children) for parent, children in self.children_by_parent.items()}

    def add_relation(self, parent, child, relation="contain", source=None):
        new_relation = {
            "parent": parent,
            "child": child,
            "relation": relation,
            "source": [] if source is None else source,
        }
        key = (parent, child, relation)
        if key not in {(item["parent"], item["child"], item["relation"]) for item in self.h_relations}:
            self.h_relations.append(new_relation)
            self.all_relations.append(new_relation)
            self._rebuild_indexes()

    def find_child(self, variations):
        children = []
        for variation in self._as_list(variations):
            children.extend(self.children_by_parent.get(variation, []))
        return children

    def find_parent(self, variations):
        parents = []
        for variation in self._as_list(variations):
            parents.extend(self.parents_by_child.get(variation, []))
        return parents

    def find_brother(self, variations):
        variation_list = self._as_list(variations)
        siblings = []
        for parent in self.find_parent(variation_list):
            siblings.extend(self.find_child(parent))
        return [node for node in dict.fromkeys(siblings) if node not in variation_list]

    def find_descendant(self, variations):
        descendants = []
        seen = set()
        queue = deque(self.find_child(variations))
        while queue:
            node = queue.popleft()
            if node in seen:
                continue
            seen.add(node)
            descendants.append(node)
            queue.extend(self.children_by_parent.get(node, []))
        return descendants

    def find_ancestor(self, variations):
        ancestors = []
        seen = set()
        queue = deque(self.find_parent(variations))
        while queue:
            node = queue.popleft()
            if node in seen:
                continue
            seen.add(node)
            ancestors.append(node)
            queue.extend(self.parents_by_child.get(node, []))
        return ancestors

    def find_descendant_(self, variations):
        return self.find_descendant(variations)

    def get_leaf_node(self):
        return self.nodes - set(self.children_by_parent)

    def find_leaf(self, variation):
        leaf_nodes = self.get_leaf_node()
        return [node for node in self.find_descendant(variation) if node in leaf_nodes]

    @staticmethod
    def _has_alternate_path(parent, child, children_by_parent):
        queue = deque(children_by_parent.get(parent, []))
        seen = set()
        while queue:
            node = queue.popleft()
            if node == child:
                return True
            if node in seen:
                continue
            seen.add(node)
            queue.extend(children_by_parent.get(node, []))
        return False

    def get_non_redundant_relations(self, relations=None):
        relations = self.h_relations if relations is None else relations
        non_redundant_relations = []

        for relation in relations:
            parent = relation["parent"]
            child = relation["child"]
            children_by_parent = defaultdict(list)
            for candidate in relations:
                if candidate is relation:
                    continue
                if candidate["parent"] == parent and candidate["child"] == child:
                    continue
                children_by_parent[candidate["parent"]].append(candidate["child"])

            if not self._has_alternate_path(parent, child, children_by_parent):
                non_redundant_relations.append(relation)

        return non_redundant_relations

    def save_relation_neo4j(
        self,
        node_file=NODE_PATH,
        edge_file=EDGE_PATH,
        with_info=False,
        node_with_info_file=NODE_WITH_INFO_PATH,
        edge_with_info_file=EDGE_WITH_INFO_PATH,
    ):
        nodes = sorted(self.nodes)
        node_ids = {node: f"e{index}" for index, node in enumerate(nodes)}

        Path(node_file).parent.mkdir(parents=True, exist_ok=True)
        with open(node_file, "w", newline="", encoding="utf-8") as fn:
            writer = csv.writer(fn)
            writer.writerow(["entity:ID", "name", ":LABEL"])
            for node in nodes:
                writer.writerow([node_ids[node], node, "Entity"])

        Path(edge_file).parent.mkdir(parents=True, exist_ok=True)
        with open(edge_file, "w", newline="", encoding="utf-8") as fe:
            writer = csv.writer(fe)
            writer.writerow([":START_ID", ":END_ID", ":TYPE"])
            for relation in self.h_relations:
                writer.writerow([
                    node_ids[relation["parent"]],
                    node_ids[relation["child"]],
                    relation["relation"],
                ])

        if not with_info:
            return

        Path(node_with_info_file).parent.mkdir(parents=True, exist_ok=True)
        with open(node_with_info_file, "w", newline="", encoding="utf-8") as fn:
            writer = csv.writer(fn)
            writer.writerow(["entity:ID", "name", ":LABEL", "variant_info"])
            for node in nodes:
                variant = self.variant_dict.get(node)
                writer.writerow([
                    node_ids[node],
                    node,
                    "Entity",
                    json.dumps(variant, ensure_ascii=False) if variant else "",
                ])

        Path(edge_with_info_file).parent.mkdir(parents=True, exist_ok=True)
        with open(edge_with_info_file, "w", newline="", encoding="utf-8") as fe:
            writer = csv.writer(fe)
            writer.writerow([
                ":START_ID",
                ":END_ID",
                ":TYPE",
                "START_name",
                "END_name",
                "relation",
                "source",
            ])
            for relation in self.h_relations:
                parent = relation["parent"]
                child = relation["child"]
                writer.writerow([
                    node_ids[parent],
                    node_ids[child],
                    relation["relation"],
                    parent,
                    child,
                    relation["relation"],
                    json.dumps(relation.get("source", []), ensure_ascii=False),
                ])

    def save_json(self, output_file=NON_REDUNDANT_HIERARCHY_PATH):
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(self.h_relations, f, ensure_ascii=False, indent=4)

    def print_hierarchy(self, parent_id=None, level=0):
        if parent_id is None:
            all_children = set(self.parents_by_child)
            root_nodes = sorted(node for node in self.nodes if node not in all_children)
            for root in root_nodes:
                self.print_hierarchy(root, level)
            return

        print("    " * level + parent_id)
        for child in self.children_by_parent.get(parent_id, []):
            self.print_hierarchy(child, level + 1)


def load_hierarchy(path=DEFAULT_HIERARCHY_PATH):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def get_neo4j(
    hierarchy_file=DEFAULT_HIERARCHY_PATH,
    node_file=NODE_PATH,
    edge_file=EDGE_PATH,
    non_redundant_file=NON_REDUNDANT_HIERARCHY_PATH,
    with_info=False,
    node_with_info_file=NODE_WITH_INFO_PATH,
    edge_with_info_file=EDGE_WITH_INFO_PATH,
):
    relation_graph = VariationHierarchicalRelation(load_hierarchy(hierarchy_file))
    relation_graph.save_json(non_redundant_file)
    relation_graph.save_relation_neo4j(
        node_file=node_file,
        edge_file=edge_file,
        with_info=with_info,
        node_with_info_file=node_with_info_file,
        edge_with_info_file=edge_with_info_file,
    )
    return relation_graph



if __name__ == "__main__":
    relation_dict = get_neo4j(with_info=True)
    # relation_dict.print_hierarchy("EGFR mutation", 3)
