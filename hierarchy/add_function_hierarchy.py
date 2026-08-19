import csv
import sys
from dataclasses import asdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

FUNCTION_HIERARCHY_CONFIG = [
    {
        "source": "Oncokb",
        "file_path": PROJECT_ROOT / "data" / "sources" / "Oncokb" / "oncokb_variant_function_annotation.csv",
        "delimiter": ",",
        "gene": "gene",
        "child": "variant",
        "parent": "knownEffect",
    },
    {
        "source": "CGI",
        "file_path": PROJECT_ROOT / "data" / "sources" / "CGI" / "catalog_of_validated_oncogenic_mutations_latest" / "catalog_of_validated_oncogenic_mutations.tsv",
        "delimiter": "\t",
        "gene": "gene",
        "child": "protein",
        "parent": "oncogenic",
        "fixed_parent": True,
    },
    {
        "source": "PMC4232638",
        "file_path": PROJECT_ROOT / "data" / "sources" / "PMC4232638" / "13059_2014_484_MOESM2_ESM.csv",
        "delimiter": ",",
        "gene": "Gene",
        "child": "Amino acid change",
        "parent": "Type",
    },
]


SOP_FUNCTION_RULES = {
    "fusion": {
        "parent": "Likely Gain-of-function",
        "source": "OncoKB-SOP",
    },
    "truncating": {
        "parent": "Likely Loss-of-function",
        "source": "OncoKB-SOP",
    },
    "amplification": {
        "parent": "Gain-of-function",
        "source": "OncoKB-SOP",
    },
}


def get_gene(variant):
    return variant.get("location", {}).get("gene", "")


def make_edge(parent, child, source):
    return {
        "parent": parent,
        "child": child,
        "relation": "contain",
        "source": [source],
    }


def get_sop_function_hier(variants):
    from variant_hierarchy.normalize.variant_parser import parse_variant

    hier_list = []
    for child in variants:
        gene = get_gene(child)
        alteration_type = child.get("alteration", {}).get("type")
        rule = SOP_FUNCTION_RULES.get(alteration_type)
        if not gene or not rule:
            continue

        parent = parse_variant(rule["parent"], gene, variant_class=None)
        if parent:
            hier_list.append(make_edge(parent, child, rule["source"]))
    return hier_list


def get_function_variants(variants):
    result = []
    seen = set()
    for edge in get_sop_function_hier(variants):
        parent = edge["parent"]
        display_name = parent.get("display_name", "")
        if display_name and display_name not in seen:
            seen.add(display_name)
            result.append(parent)
    return result


def get_variants_with_more_gene_variants(variants):
    from variant_hierarchy.normalize.variant import build_gene_variant

    gene_list = sorted({
        variant.get("location", {}).get("gene")
        for variant in variants
        if variant.get("location", {}).get("gene")
        and variant.get("location", {}).get("gene") not in {"?", "v"}
    })
    existing_display_names = {
        variant.get("display_name")
        for variant in variants
        if variant.get("display_name")
    }

    result = []
    for gene in gene_list:
        new_variant = asdict(build_gene_variant(gene, "mutation"))
        new_variant["source"] = ["VH"]
        display_name = new_variant.get("display_name")
        if display_name not in existing_display_names:
            result.append(new_variant)
            existing_display_names.add(display_name)
    return result


def get_function_hier(variants=None):
    from variant_hierarchy.normalize.variant_parser import parse_variant

    hier_list = []
    for config in FUNCTION_HIERARCHY_CONFIG:
        with open(config["file_path"], "r", encoding="utf-8", newline="") as f:
            records = csv.DictReader(f, delimiter=config["delimiter"])
            for record in records:
                gene = str(record.get(config["gene"], "")).strip()
                child_text = str(record.get(config["child"], "")).strip()
                if config.get("fixed_parent"):
                    parent_text = config["parent"]
                else:
                    parent_text = str(record.get(config["parent"], "")).strip()

                parent = parse_variant(parent_text, gene, variant_class=None)
                child = parse_variant(child_text, gene, variant_class=None)
                if parent and child:
                    parent["source"] = [config["source"]]
                    child["source"] = [config["source"]]
                    hier_list.append(make_edge(parent, child, config["source"]))
    if variants:
        hier_list.extend(get_sop_function_hier(variants))
    return hier_list
