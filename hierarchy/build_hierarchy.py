import json
import sys
from pathlib import Path

UTA_DSN = "postgresql://uta_admin:uta@localhost:5432/uta/uta_20241220"
SEQREPO_ROOT_DIR = "/usr/local/share/seqrepo/2024-12-20"
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from add_function_hierarchy import (
    get_function_hier,
    get_function_variants,
    get_variants_with_more_gene_variants,
)
from variant_hierarchy.contains import contains_variant
VARIANTS_PATH = PROJECT_ROOT / "data" / "variants" / "variants.json"
HIERARCHY_PATH = PROJECT_ROOT / "data" / "hierarchy" / "hierarchy.json"


def deduplicate_hierarchy(hierarchy):
    seen = set()
    result = []
    index = {}
    for edge in hierarchy:
        key = (
            edge["parent"]["display_name"],
            edge["child"]["display_name"],
            edge.get("relation", "contain"),
        )
        if key in seen:
            sources = index[key].setdefault("source", [])
            for source in edge.get("source", []):
                if source not in sources:
                    sources.append(source)
            continue
        seen.add(key)
        edge["source"] = list(edge.get("source", []))
        index[key] = edge
        result.append(edge)
    return result


def deduplicate_variants(variants):
    seen = set()
    result = []
    for variant in variants:
        display_name = variant.get("display_name", "")
        if not display_name or display_name in seen:
            continue
        seen.add(display_name)
        result.append(variant)
    return result


def get_hier(variants):
    # hgvs_count = 0
    hier_list = []
    for var1 in variants:
        if var1['variant_class'] != "HGVS Variant" and var1['variant_class'] != "Uncertain Variant" : # and var1['variant_class'] != "Function Variant"
            print(f"----------------{var1['display_name']} start!------------------------")
            for var2 in variants:
                # print(f"---{var2['displayName']}!---")
                if contains_variant(var1, var2):
                    hier_list.append({
                        "parent": var1,
                        "child": var2,
                        "relation": "contain",
                        "source": ["VH"],
                    })
                    print(f"{var1['display_name']}>{var2['display_name']}")
        # else:
        #     print(f"{var1['display_name']} is already a HGVS variant.")
            # hgvs_count += 1
            # if hgvs_count % 50 == 0:
            #     print(f"{var2['display_name']} is already a HGVS variant.")
    return hier_list

    # merged_dict = {}
    # for item in hier_list:
    #     attr = item['parent'] + '>' + item['child']
    #     merged_dict[attr] = item.copy()
    #
    # unique_hier_list = list(merged_dict.values())
    # unique_hier_list = [convert_set_to_list(item) for item in unique_hier_list]
    #
    #
    # # unique_list = [dict(t) for t in {frozenset(d.items()) for d in hier_list}]
    # with open("../data/hierarchy/hierarchy_unique.json", "w", encoding="utf-8") as f:
    #     json.dump(unique_hier_list, f, ensure_ascii=False, indent=4)



def main():
    with open(VARIANTS_PATH, "r", encoding="utf-8") as f:
        variants = json.load(f)
    variants = deduplicate_variants(
        variants
        + get_function_variants(variants)
        + get_variants_with_more_gene_variants(variants)
    )
    hier = deduplicate_hierarchy(get_hier(variants) + get_function_hier(variants))
    HIERARCHY_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(HIERARCHY_PATH, "w", encoding="utf-8") as f:
        json.dump(hier, f, ensure_ascii=False, indent=4)
    return hier


def test_contains_variant(variant_pairs):
    from variant_hierarchy.normalize.variant_parser import parse_variant

    parent =  variant_pairs["parent"]
    child =  variant_pairs["child"]
    parent_variant = parse_variant(parent)
    child_variant = parse_variant(child)
    result = contains_variant(parent_variant, child_variant)
    return result




if __name__ == '__main__':
    main()





