import re
import json
import sys
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from variant_hierarchy.normalize.variant_parser import parse_variant
from logic import LogicParser

from source_mappings.loader import load_yaml_mapping
COSMIC_MAPPING = load_yaml_mapping("cosmic_biomarker_aliases.yaml")


def filter_brackets_content(input_str):
    def condition(content):
        parser = LogicParser()
        parsed_result = parser.parse(content.strip())
        elements = parser.extract_elements(parsed_result)
        if len(elements) > 1 or ";" in content:
            return False
        else:
            return True

    def replace_func(match):
        content = match.group(2)
        chars = match.group(0)
        if condition(content):
            return ""
        else:
            return chars
    result = re.sub(r"(.)\(([^()]*?)\)", replace_func, input_str) #  r"(.)\((.*?)\)"
    return result


def norm_variant(biomarker):
    biomarker = COSMIC_MAPPING.get(biomarker.strip(), biomarker.strip())
    biomarker = filter_brackets_content(biomarker)

    # ----------------------------------------split variant------------------------------------------
    # and or
    parser = LogicParser()
    parsed_result = parser.parse(biomarker)
    elements = parser.extract_elements(parsed_result)
    if len(elements) > 1:
        results = []
        for element in elements:
            results.extend(norm_variant(element.strip()))
        return results

    # V600E/K
    match = re.match(r"([A-Za-z0-9_]+)_([A-Za-z])(\d+)([A-Za-z/]+)", biomarker)
    if match and "/" in biomarker:
        gene = match.group(1)
        ref_aa = match.group(2)
        position = match.group(3)
        alt_aa_set = match.group(4)
        alt_aa_list = alt_aa_set.split("/")  # ? "/" ??
        elements = [f"{gene}_{ref_aa}{position}{alt_aa}" for alt_aa in alt_aa_list]
        results = []
        for element in elements:
            results.extend(norm_variant(element.strip()))
        return results

    # EGFR_no_V600
    match = re.match(r"([a-z0-9_]+)_(no|not)_(.+)", biomarker, re.IGNORECASE)
    if match:
        gene = match.group(1)  # ???
        alteration = match.group(3)  # ??????
        biomarker = f"{gene}_{alteration}"


    # ----------------------------------------norm variant------------------------------------------

    parsed_variant = parse_variant(biomarker, variant_class=None)
    if parsed_variant:
        return [parsed_variant]
    print(f"Unable to parse variant from variantName (variantName={biomarker})")
    return []
    # raise Exception(f"Unable to parse variant from variantName (variantName={biomarker})")



def get_variant():
    file_path = "../data/sources/Cosmic/Actionability_AllData_v13_GRCh37.tsv"
    df = pd.read_csv(file_path, sep='\t')
    records = list(df['MUTATION_REMARK'].unique())

    variants = []
    for record in records:
        items = norm_variant(record)
        for item in items:
            variants.append({**item, **{"source": ["COSMIC"]},
                             **{"raw": {"alteration": record}}})

    unique_variants = list({json.dumps(v, sort_keys=True): v for v in variants}.values())
    # gene_variants = [variant for variant in unique_variants if variant['variantClass'] == "GeneVariant"]
    # unresolved_variants = [variant for variant in unique_variants if variant['variantClass'] == "UnresolvedVariant"]
    with open("cosmic_variants.json", "w", encoding="utf-8") as f:  # ../data/variants/
        json.dump(unique_variants, f, ensure_ascii=False, indent=4)
    return unique_variants


if __name__ == '__main__':
    variants = get_variant()

    # res = norm_variant("1p/19g")
    # variants = [item for variant in variants for item in variant]
    # dic = {'aa_change': '61', 'displayName': 'KRAS 61', 'reference1': 'KRAS'}
    # if dic in variants:
    #     print('error')

    # file_path = "/mnt/Data.Common/Chenh/PycharmProjects/InformationExtraction/dataset/graphkb_dataset/raw/Cosmic/Actionability_AllData_v13_GRCh37.tsv"
    # df = pd.read_csv(file_path, sep='\t')
    # unspecified_rows = df[df['MUTATION_REMARK'] == 'unspecified']
    # variants = list(df['MUTATION_REMARK'].unique())
    # variants = [norm_variant(item) for item in variants] # for variant in variants
    # variants = [item for variant in variants for item in variant]
    # test(variants)


    # typ_ = [item for item in variants if "(" not in item]
    # records = df.to_dict(orient='records')
    # res = norm_variant('BCR-ABL ?:? COSF1783 and IKZF1_unspecified and MS4A1_unspecified')
    # res1 = filter_brackets_content('t(11;14) (CCND1)')
    # res = norm_variant('ERBB2_unspecified and ((ESR1_unspecified or PGR_unspecified or (ESR1_unspecified and PGR_unspecified)) and (CCND1_unspecified or CCND1_unspecified or CCNE2_unspecified or CDK2_unspecified or CDK4_unspecified or CDK6_unspecified)') # EGFR_no_Exon19del and (EGFR_V600E/K or EGFR_exon20-21)
    # res = norm_variant('11q_deletion_(ETS1) and MS4A1_unspecified')  # EGFR_no_Exon19del and (EGFR_V600E/K or EGFR_exon20-21)

    # variants = [item for item in variants if 'Exon' in item]

    # variants = [split_variant(item) for item in variants]

    #
    # variant2 = [norm_variant(record['MUTATION_REMARK']) for record in records]
    
