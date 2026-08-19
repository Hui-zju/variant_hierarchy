# VariantHierarchy

VariantHierarchy is an open-source computational framework for semantic normalization and hierarchical organization of somatic categorical variants in precision oncology. Precision-oncology evidence is often described at different levels of abstraction: a knowledge base may annotate evidence for `EGFR exon 19 in-frame deletion` while a patient report contains `EGFR E746_A750del`, or evidence may be assigned to `BRAF V600X` while the patient variant is `BRAF V600E`. VariantHierarchy makes these cross-granularity relationships explicit, machine-readable, and reusable.

The framework represents each normalized variant concept across three semantic dimensions: genomic location, molecular alteration, and functional interpretation. It then infers directed containment relationships from broader concepts to more specific compatible concepts. The core implementation is maintained in `src/variant_hierarchy` and is distributed on PyPI as the `variant-hierarchy` package.

## Research Scope

VariantHierarchy focuses on somatic categorical variants from cancer precision oncology resources, including protein-level sequence variants, amino-acid site and region categories, exon-level categories, gene-level mutation classes, copy-number alterations, structural variants, expression and methylation states, and functional or oncogenicity assertions.

The current study integrates variant concepts from COSMIC, OncoKB, CIViC, and the Cancer Genome Interpreter (CGI). After normalization and deduplication, 12,472 unique source-derived variant concepts were retained. The final hierarchy contains semantically inferred containment relationships and knowledge-supported functional relationships, followed by redundancy removal to produce a non-redundant hierarchy graph.

## Framework Overview

VariantHierarchy is organized around four methodological components:

| Component | Role |
| --- | --- |
| Semantic representation | Represents each variant concept using genomic location, molecular alteration, and functional interpretation. |
| Controlled vocabularies | Standardize recurrent molecular alteration and functional interpretation terms and encode broader-narrower relationships. |
| Variant normalization | Converts heterogeneous source terms into structured normalized variant concepts with canonical display names. |
| Containment inference | Infers parent-child relationships between compatible concepts and constructs a directed hierarchy graph. |

## Semantic Representation

Each normalized variant concept is represented by a display name and three semantic dimensions:

| Dimension | Meaning | Examples |
| --- | --- | --- |
| Genomic location | The biological scope at which the variant is described. | gene, exon, amino-acid site, protein HGVS position/edit |
| Molecular alteration | The molecular event or alteration class. | mutation, missense, in-frame deletion, amplification, fusion, methylation |
| Functional interpretation | Biological consequence or oncogenicity assertion. | gain-of-function, loss-of-function, oncogenic mutation, inconclusive mutation |

The framework distinguishes concrete variants from categorical variants. In this project, a concrete variant refers to a protein-level sequence alteration that can be represented using HGVS protein nomenclature, such as `BRAF V600E`. Categorical variants are broader concepts that encompass multiple biologically related variants, such as `BRAF V600X`, `BRAF mutation`, or `BRAF oncogenic mutation`.

Variant concepts are further grouped into semantic categories including p.HGVS variants, categorical sequence variants, copy-number alteration categories, structural alteration categories, expression-state categories, methylation-state categories, and functional assertion categories.

## Variant Normalization

Variant normalization converts heterogeneous source terms into the semantic representation described above. The procedure includes source-specific preprocessing, gene normalization, variant parsing, and semantic assignment.

The normalization process covers:

- protein-level substitutions, deletions, duplications, delins, frameshift and nonsense variants;
- amino-acid site and region-level variant classes;
- exon-level mutation, deletion, insertion, splice, and exon-skipping categories;
- gene-level mutation and copy-number categories;
- structural variants including fusions, rearrangements, and translocations;
- expression and methylation states;
- functional and oncogenicity assertions.

Clinical and database shorthand is handled conservatively. For example, incomplete insertion shorthand such as `ALK T1151ins` is represented as an amino-acid position-level in-frame insertion rather than as a complete coordinate-resolvable protein HGVS edit.

Gene normalization uses a local alias table by default. If a VICC Gene Normalizer service is available at:

```text
http://localhost:8001
```

the runtime can use it as an enhanced gene-normalization backend and fall back to local aliases when needed.

## Hierarchy Construction

Containment inference determines whether one normalized variant concept semantically encompasses another. Candidate parent-child pairs are compared jointly across genomic location, molecular alteration, and functional interpretation. A parent contains a child when no dimension is incompatible and at least one dimension is broader.

Examples of containment relationships include:

- `BRAF V600X` contains `BRAF V600E`;
- `EGFR exon 19 in-frame deletion` contains `EGFR E746_A750del`;
- `TP53 mutation` contains `TP53 R175H`;
- `ALK fusion` contains compatible specific fusions such as `EML4::ALK fusion`;
- `PTEN loss-of-function` contains compatible loss-of-function mechanisms.

Location containment uses structured gene, exon, amino-acid, and p.HGVS representations. When precise positional comparison is required, local genomic infrastructure can be used to convert protein or region-level descriptions into comparable transcript-based intervals. This dynamic containment pathway may require:

- [VICC Gene Normalizer](https://github.com/cancervariants/gene-normalization) at `http://localhost:8001` for enhanced gene normalization;
- [UTA, the Universal Transcript Archive](https://github.com/biocommons/uta), for transcript alignment data;
- [SeqRepo](https://github.com/biocommons/biocommons.seqrepo) sequence data for local biological sequence retrieval.

A typical local UTA DSN is:

```text
postgresql://uta_admin:uta@localhost:5432/uta/uta_20241220
```

Precomputed hierarchy queries do not require local VICC, UTA, or SeqRepo services.

## Python Package

Install the runtime package from PyPI:

```bash
pip install variant-hierarchy
```

The package exposes precomputed hierarchy queries and lightweight normalization:

```python
from variant_hierarchy import (
    ancestors,
    children,
    descendants,
    gene_hierarchy,
    hierarchy,
    normalize_variant,
    parents,
    variant_hierarchy,
)

normalize_variant("BRAF V600")
parents("BRAF V600D")
children("BRAF V600X")
ancestors("EGFR E746_A750del")
descendants("BRAF V600")
```

The package supports both gene-level and variant-level hierarchy queries:

```python
gene_hierarchy("BRAF")
variant_hierarchy("BRAF V600D")

# Auto mode: gene input returns a gene hierarchy; variant input returns a variant hierarchy.
hierarchy("BRAF")
hierarchy("BRAF V600D")
```

Dynamic containment can be used when the required genomic infrastructure is available:

```python
from variant_hierarchy import contains_variant

contains_variant(
    "EGFR exon 19 in-frame deletion",
    "EGFR E746_A750del",
    uta_dsn="postgresql://uta_admin:uta@localhost:5432/uta/uta_20241220",
)
```

If the required genomic dependencies, UTA service, or SeqRepo data are unavailable, `contains_variant` raises `GenomicBackendUnavailable`. This does not affect precomputed hierarchy queries such as `parents()`, `children()`, `gene_hierarchy()`, or `variant_hierarchy()`.

## Hierarchy Visualization

An interactive hierarchy browser is available through GitHub Pages:

```text
https://hui-zju.github.io/variant_hierarchy/visualization/
```

The visualization supports gene hierarchy queries, variant node queries, ancestor paths, descendants, branch expansion and collapse, pan and zoom, and inspection of normalized variant metadata.



## Data Sources and Hierarchy Construction

The original third-party datasets are not redistributed in this repository. To reproduce the full analyses, download the following files from the original providers and place them under the specified local paths.

### Knowledge-base Sources

| Source | Download from | Local path                                                                               |
| --- | --- |------------------------------------------------------------------------------------------|
| OncoKB biomarker-drug associations | OncoKB | `data/sources/Oncokb/oncokb_biomarker_drug_associations.tsv`                             |
| CIViC Variant Summaries, 01-Nov-2024 | CIViC | `data/sources/CIVIC/01-Nov-2024-VariantSummaries.tsv`                                    |
| COSMIC Actionability v13, GRCh37 | COSMIC | `data/sources/Cosmic/Actionability_AllData_v13_GRCh37.tsv`                               |
| CGI biomarkers | Cancer Genome Interpreter | `data/sources/CGI/cancergenomeinterpreter.org_data_biomarkers_cgi_biomarkers_latest.tsv` |

### Functional-annotation Sources

| Source | Download from | Local path |
| --- | --- | --- |
| OncoKB variant function annotations | OncoKB | `data/sources/Oncokb/oncokb_variant_function_annotation.csv` |
| CGI validated oncogenic mutations | Cancer Genome Interpreter | `data/sources/CGI/catalog_of_validated_oncogenic_mutations_latest/catalog_of_validated_oncogenic_mutations.tsv` |
| Functionally validated cancer-related missense mutations | Supplementary data from PMID/PMC `PMC4232638` | `data/sources/PMC4232638/13059_2014_484_MOESM2_ESM.csv` |

### Construction Workflow

The full hierarchy is built from the source files in three main steps:

1. Run `normalize/variant_norm.py` to normalize heterogeneous source variant terms and generate standardized variant records under `data/variants/`.
2. Run `hierarchy/build_hierarchy.py` to infer containment and function-supported relationships and write the full hierarchy graph to `data/hierarchy/hierarchy.json`.
3. Run `hierarchy/hierarchy_class.py` to remove redundant transitive relations and export graph files, including `data/hierarchy/hierarchy_non_redundant.json` and Neo4j-ready node/edge tables under `data/graph/`.


## Evaluation

VariantHierarchy is evaluated at two levels.

First, curated challenge sets are used to evaluate variant normalization and containment inference. The challenge sets are stored under `data/challenge_sets/`: `normalization_challenge_set.tsv` samples heterogeneous variant terms across predefined categories, and `hierarchy_relation_challenge_set.tsv` samples candidate parent-child pairs across positive and negative relation types. Expert-reviewed labels are used as the reference standard. The stratum definitions and file schema are documented in `data/challenge_sets/STRATA.md`.

Second, application-level utility is evaluated by comparing exact matching with hierarchy-aware matching for precision oncology evidence retrieval. Patient-variant records from GENIE/MSKCC and MSK-IMPACT cohorts are normalized and matched against evidence-bearing variant concepts from COSMIC, OncoKB, CIViC, and CGI. Hierarchy-aware matching additionally allows patient variants to retrieve evidence assigned to compatible broader ancestor concepts.

### Patient Cohorts

| Cohort | Download from | Required local files |
| --- | --- | --- |
| GENIE/MSKCC OncoKB-annotated cohort | `https://github.com/oncokb/oncokb-datahub/tree/main/PUBLICATION/2023/CANCER_DISCOVERY` | `data/patient/genie_mskcc_processed/genie_mskcc_samples_with_2017_oncokb_annotation.txt`<br>`data/patient/genie_mskcc_processed/genie_mskcc_samples_with_2022_oncokb_annotation.txt` |
| MSK-IMPACT 2017 | `https://github.com/cBioPortal/datahub/tree/master/public/msk_impact_2017` | `data/patient/msk_impact_2017/data_clinical_patient.txt`<br>`data/patient/msk_impact_2017/data_clinical_sample.txt`<br>`data/patient/msk_impact_2017/data_cna.txt`<br>`data/patient/msk_impact_2017/data_sv.txt`<br>`data/patient/msk_impact_2017/data_mutations_api_full.json` |

For MSK-IMPACT 2017, mutation records are obtained from the cBioPortal API molecular profile `msk_impact_2017_mutations` and cached locally as `data_mutations_api_full.json`.
