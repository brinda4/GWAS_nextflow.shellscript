# GWAS Analysis Pipeline Methodology

## 1. Population Demarcation and Cohort Mapping
In the NHGRI-EBI GWAS Catalog, genetic associations are reported at the level of individual studies (`STUDY ACCESSION`). Individual study participants are grouped into ancestry categories defined by EMBL-EBI ontology standards:
- **Discovery Stage (`STAGE == 'initial'`)**: To ensure genetic allele frequencies and effect sizes correspond to the population under study, ancestry categories are derived from the discovery stage. If a study lacks stage metadata, all stages are evaluated.
- **Ancestry Harmonization**:
  - **European Cohorts**: Defined by membership in `European`.
  - **South Asian Cohorts**: Defined by membership in `South Asian`.
  - **Multi-Ancestry / Shared Cohorts**: Studies incorporating both European and South Asian cohorts (e.g., global trans-ethnic meta-analyses).

## 2. Allele Frequency Parsing and Quality Thresholds
- **Target Interval**: Low-frequency variants are formally defined as alleles with population frequency between 0.5% and 2.0% inclusive:
  \[
  0.005 \le \text{RAF} \le 0.020
  \]
- **Missing Value Handling**:
  - `NR` (Not Reported): 51.7% of associations in the GWAS Catalog lack an explicit risk allele frequency in the summary table. These records are cataloged in the missingness audit and excluded from allele-frequency distribution statistics to prevent bias.
  - Text Annotations: Cases where frequencies are annotated with population labels (e.g. `0.015 (European)`) are parsed with regular expressions to retrieve the numeric scalar.

## 3. Gene-Disease-Population Mapping
- **Gene Extraction**:
  - Ensembl mapped genes (`MAPPED_GENE`) are prioritized for functional objectivity. If absent, author-reported genes (`REPORTED GENE(S)`) are used.
  - Intergenic annotations (`GENE1 - GENE2`) and transcript clusters are tokenized into discrete gene symbols while preserving internal punctuation for recognized biological entities (e.g., `HLA-A`).
- **Population Classification**:
  Each gene is categorized into one of three ancestral regimes:
  1. **European Only**: Discovered exclusively in European GWAS cohorts.
  2. **South Asian Only**: Discovered exclusively in South Asian GWAS cohorts.
  3. **Both Populations (Shared)**: Identified in both European and South Asian studies, demonstrating ancestral cross-population generalizability.
