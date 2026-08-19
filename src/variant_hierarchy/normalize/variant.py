from dataclasses import dataclass, field
from typing import Dict, Any, Optional
from .aa_position import get_hgvs_p_aa_range
from .ontology import VARIANT_ALTERATION_ONTOLOGY, VARIANT_FUNCTION_ONTOLOGY

@dataclass
class Variant:
    display_name: str
    variant_class: str  #HGVS Variant/Sequence Variant/Copy Number Variant/Structural Variant/Expression Variant/Epigenetic Variant/Function Variant
    alteration_coverage_of_location: Optional[str] = None  #Whole Location/Within Location/At Location/Between Locations/Unspecified
    location_type: Optional[str] = None  # Gene/Exon/Intron/Domain/Amino Acid/Nucleotide
    location: Optional[Dict[str, Any]] = field(default_factory=dict)
    alteration: Dict[str, Any] = field(default_factory=dict)
    function: Dict[str, Any] = field(default_factory=dict)


def get_display_name(
    gene: Optional[str] = None,
    location_label: Optional[str] = None,
    alteration_label: Optional[str] = None,
    function_label: Optional[str] = None
):
    parts = []

    if gene:
        parts.append(gene)

    if location_label:
        parts.append(location_label)

    if function_label:
        parts.append(function_label)

    if function_label and function_label.endswith("mutation") and alteration_label == "mutation":
        alteration_label = None

    if alteration_label:
        parts.append(alteration_label)

    return " ".join(parts)

def get_location_coverage(alteration_type, location_type):
    variant_class = VARIANT_ALTERATION_ONTOLOGY[alteration_type]["variant_class"] if alteration_type else "Unspecified"

    if location_type == "Amino Acid":
        coverage = "At Location"
    else:
        if alteration_type == "mutation":
            coverage = "Unspecified"
        elif alteration_type == "fusion":
            coverage = "Between Locations"

        elif variant_class == "Sequence Variant" and alteration_type != "deletion":
            coverage = "Within Location"
        else:
            coverage = "Whole Location"
    return coverage




def build_p_hgvs_variant(
    gene: str,
    p_hgvs_posedit: str,  # T790M E746_A750del  M695fs*26
    alteration_type: str
):
    display_name = get_display_name(
        gene=gene,
        location_label=p_hgvs_posedit
    )

    return Variant(
        display_name=display_name,
        variant_class="HGVS Variant",
        alteration_coverage_of_location="At Location",
        location_type="Amino Acid",
        location={
            "gene": gene,
            "p_hgvs_posedit": p_hgvs_posedit,
        },
        alteration={
            "type": alteration_type,
            "p_hgvs_posedit": p_hgvs_posedit,
        },
        function={}
    )

def build_amino_acid_variant(
    gene: str,
    aa_range: str,   #  "pos"  "pos1-pos2"
    alteration_type: str,
):
    type_definition = VARIANT_ALTERATION_ONTOLOGY[alteration_type] if alteration_type else None
    alteration_label = type_definition["canonical_label"] if alteration_type else None

    aa_pos = get_hgvs_p_aa_range(gene, aa_range)
    if aa_pos:
        if alteration_type == "missense":
            location_label = aa_pos + "X"
            alteration_label = None
        else:
            location_label = aa_pos
            alteration_label = alteration_label

        display_name = get_display_name(
            gene=gene,
            location_label=location_label,
            alteration_label=alteration_label
        )
    else:
        display_name = get_display_name(
            gene=gene,
            location_label=f"amino acid {aa_range}",
            alteration_label=alteration_label
        )

    return Variant(
        display_name=display_name,
        variant_class="Sequence Variant",
        alteration_coverage_of_location="At Location",
        location_type="Amino Acid",
        location={
            "gene": gene,
            "amino_acid": aa_range,
        },
        alteration={
            "type": alteration_type,
        },
        function={}
    )


def build_exon_variant(
    gene: str,
    exon_range: str, #  "pos"  "pos1-pos2"
    alteration_type: Optional[str] = None,
    function: Optional[str] = None
):
    type_definition = VARIANT_ALTERATION_ONTOLOGY[alteration_type] if alteration_type else None
    function_definition = VARIANT_FUNCTION_ONTOLOGY[function] if function else None

    alteration_label = type_definition["canonical_label"] if alteration_type else None
    function_label = function_definition["canonical_label"] if function else None
    display_name = get_display_name(
        gene=gene,
        location_label=f"exon {exon_range}",
        function_label=function_label,
        alteration_label=alteration_label
    )

    return Variant(
        display_name=display_name,
        variant_class=type_definition["variant_class"] if alteration_type else "Function Variant",
        alteration_coverage_of_location=get_location_coverage(alteration_type, location_type="Exon"),
        location_type="Exon",
        location={"gene": gene,"exon": exon_range,},
        alteration={"type": alteration_type} if alteration_type else {"type": "mutation"},
        function={function_definition["functional_domain"]: function_definition["functional_value"]} if function else {}
    )

def build_gene_variant(
    gene: str,
    alteration_type: Optional[str] = None,
    function: Optional[str] = None
):
    type_definition = VARIANT_ALTERATION_ONTOLOGY[alteration_type] if alteration_type else None
    function_definition = VARIANT_FUNCTION_ONTOLOGY[function] if function else None

    alteration_label = type_definition["canonical_label"] if alteration_type else None
    function_label = function_definition["canonical_label"] if function else None
    display_name = get_display_name(
        gene=gene,
        function_label=function_label,
        alteration_label=alteration_label
    )

    return Variant(
        display_name=display_name,
        variant_class=type_definition["variant_class"] if alteration_type else "Function Variant",
        alteration_coverage_of_location=get_location_coverage(alteration_type, "Gene"),
        location_type="Gene",
        location={"gene": gene},
        alteration={"type": alteration_type} if alteration_type else {"type": "mutation"},
        function={function_definition["functional_domain"]: function_definition["functional_value"]} if function else {}
    )


def build_copy_number_variant(
    gene: str,
    alteration_type: Optional[str] = None,
    copy_count:  Optional[int] = None,
):
    canonical_label = VARIANT_ALTERATION_ONTOLOGY[alteration_type] ["canonical_label"]
    alteration_label = canonical_label + " X" + str(copy_count) if copy_count else canonical_label
    display_name = get_display_name(
        gene=gene,
        alteration_label=alteration_label
    )

    alteration_dict = {"type": alteration_type}
    if copy_count:
        alteration_dict.update({"copyCount": copy_count})

    return Variant(
        display_name=display_name,
        variant_class="Copy Number Variant",
        alteration_coverage_of_location="Whole Location",
        location_type="Gene",
        location={ "gene": gene,},
        alteration=alteration_dict,
        function={}
    )

def build_fusion_variant(
    gene1: str,
    gene2: Optional[str] = None
):
    if gene2:
        location = {"gene": gene1, "gene2": gene2}
        display_name = f"{gene1}::{gene2} fusion"
    else:
        location = {"gene": gene1}
        display_name = f"{gene1} fusion"

    return Variant(
        display_name=display_name,
        variant_class="Structural Variant",
        alteration_coverage_of_location="Between Locations",
        location_type="Gene",
        location=location,
        alteration={"type": "fusion"},
        function= {}
    )


def build_unresolved_variant(
    variant: str,
    gene: Optional[str] = None,
):
    if gene:
        display_name = f"{gene} {variant}"
    else:
        display_name = variant

    return Variant(
        display_name=display_name,
        variant_class="Uncertain Variant",
        alteration_coverage_of_location=None,
        location_type=None,
        location={},
        alteration={},
        function= {}
    )
