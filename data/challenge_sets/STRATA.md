# Variant Challenge Set Strata

These TSV files are compact candidate challenge sets. Each row keeps only an input and the expected output derived from the current normalized variants and hierarchy logic.

## TSV columns

### normalization_challenge_set.tsv

- `id`: record id.
- `stratum_id`: difficulty type code; definitions are listed below.
- `input_text`: original source text used as normalization input.
- `expected_variant_name`: expert-approved normalized variant name.

`normalization_challenge_set.tsv` follows the Variant schema order. `normalization_challenge_set_old.tsv` preserves the previous reviewed ordering for comparison.

### unsupported_terms_taxonomy.tsv

- `id`: unsupported-term record id.
- `raw_term`: original term that the current system did not normalize.
- `reference`: source-provided gene, feature, or biomarker context.
- `category_id`: unsupported category code.
- `category_name`: unsupported category label.
- `recommendation`: whether to exclude, describe as a limitation, or consider for future parser extension.

### hierarchy_relation_challenge_set.tsv

- `id`: record id.
- `stratum_id`: difficulty type code; definitions are listed below.
- `parent_variant`: candidate parent variant name.
- `child_variant`: candidate child variant name.
- `expected_contains`: binary label for whether `parent_variant` contains `child_variant`; `1` means contains and `0` means does not contain.

`hierarchy_relation_challenge_set.tsv` is the manually reviewed hierarchy challenge set ordered according to the Variant schema and the practical evidence-retrieval direction.

## Normalization strata

The normalization strata follow the Variant schema in `normalize/variant.py`. Each stratum contains 10 reviewed samples. The top-level design follows `variant_class`; location, alteration, and function define subtype boundaries.

- `N01` HGVS Variant: concrete protein-level variants
- `N02` Sequence Variant: amino-acid site
- `N03` Sequence Variant: gene-level sequence mechanism
- `N04` Sequence Variant: region-level sequence category
- `N05` Sequence Variant: splice-related category
- `N06` Structural Variant
- `N07` Copy Number Variant
- `N08` Expression Variant and Epigenetic Variant
- `N09` Function Variant
- `N10` Broad classified gene-level variant

## Hierarchy relation strata

The hierarchy strata follow parent-child `variant_class` combinations in the Variant schema. Each stratum contains 5 reviewed samples. H01-H11 are positive containment cases (`expected_contains = 1`), and H12-H21 are negative controls (`expected_contains = 0`). The ordering emphasizes practical evidence retrieval from narrow patient-like variants to broader categorical evidence variants.

- `H01` Sequence -> HGVS: gene-level category to concrete protein variant
- `H02` Sequence -> HGVS: amino-acid site category to concrete protein variant
- `H03` Sequence -> HGVS: region-level category to concrete protein variant
- `H04` Sequence -> Sequence: gene-level category to amino-acid site category
- `H05` Sequence -> Sequence: gene-level category to sequence mechanism category
- `H06` Sequence -> Sequence: gene-level category to region-level category
- `H07` Sequence -> Sequence: region-level category to narrower region-level category
- `H08` Structural -> Structural: fusion abstraction to partner-specific fusion
- `H09` Gene-level anchor -> same-gene classified variant
- `H10` Function hierarchy: ontology or knowledge-supported compatible category
- `H11` Alteration ontology: broad mechanism to narrower mechanism
- `H12` Negative control: different gene
- `H13` Negative control: amino-acid site mismatch
- `H14` Negative control: region mismatch
- `H15` Negative control: gene-level root cross-gene
- `H16` Negative control: incompatible mechanism
- `H17` Negative control: same-site HGVS sibling
- `H18` Negative control: reversed hierarchy
- `H19` Negative control: fusion partner mismatch
- `H20` Negative control: incompatible function categories
- `H21` Negative control: function-mechanism incompatibility

## Unsupported-term taxonomy

This file is not part of the positive normalization gold standard. It is an exclusion/error-analysis table for terms that are outside the current challenge-set scope or need a future ontology/parser decision.

- `U01` Non-variant biomarker or phenotype.
- `U02` Pharmacogenomic/SNP/star-allele identifier.
- `U03` Nucleotide-level HGVS or intronic cDNA notation.
- `U04` Cytoband or chromosomal alteration.
- `U05` Protein state, localization, isoform, or expression assay term.
- `U06` Ambiguous loss/truncation ontology term.
- `U07` Potentially standardizable cancer variant form not yet parsed.
- `U08` Complex fusion/transcript/cytogenetic expression.
- `U09` Invalid, obsolete, or non-human gene/feature symbol.
- `U10` Other unsupported or ambiguous expression.

## Output files

- `normalization_challenge_set.tsv`
- `normalization_challenge_set_old.tsv`
- `hierarchy_relation_challenge_set.tsv`
- `unsupported_terms_taxonomy.tsv`
