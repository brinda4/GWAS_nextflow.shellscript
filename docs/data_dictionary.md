# GWAS Catalog Data Dictionary

This document details the complete 38-column schema of the NHGRI-EBI GWAS Catalog v1.0.2 association dataset and the 12-column schema of the Ancestry dataset.

---

## 1. GWAS Catalog Associations Dataset (`gwas_catalog_v1.0.2-associations`)

| Index | Column Header | Data Type | Domain | Description & Curation Rules |
| :---: | :--- | :--- | :--- | :--- |
| **0** | `DATE ADDED TO CATALOG` | Date (`YYYY-MM-DD`) | Study Metadata | Date the study was officially curated and released in the catalog. |
| **1** | `PUBMEDID` | Integer | Study Metadata | PubMed identification number of the published article. |
| **2** | `FIRST AUTHOR` | String | Study Metadata | Last name and initials of the primary author (e.g. `Malik R`). |
| **3** | `DATE` | Date (`YYYY-MM-DD`) | Study Metadata | Online publication date (epub date if available). |
| **4** | `JOURNAL` | String | Study Metadata | Abbreviated journal name (e.g. `Nat Genet`, `Lancet`). |
| **5** | `LINK` | URL | Study Metadata | Europe PMC / PubMed URL linking to publication. |
| **6** | `STUDY` | String | Study Metadata | Full title of the scientific manuscript. |
| **7** | `DISEASE/TRAIT` | String | Disease / Trait | Author-reported clinical condition, disease, or continuous trait examined. |
| **8** | `INITIAL SAMPLE SIZE` | String | Population | Cohort breakdown for Stage 1 (discovery), detailing sample sizes and ancestral groups. |
| **9** | `REPLICATION SAMPLE SIZE` | String | Population | Cohort breakdown for Stage 2 (replication) stages. Marked as `'NA'` if unapplicable. |
| **10** | `REGION` | String | Genomic Position | Cytogenetic band/region associated with variant (e.g. `4q25`, `6p21.3`). |
| **11** | `CHR_ID` | String | Genomic Position | Chromosome identifier (`1`-`22`, `X`, `Y`, `MT`). |
| **12** | `CHR_POS` | Integer | Genomic Position | Base-pair coordinate on reference genome build GRCh38.p14. |
| **13** | `REPORTED GENE(S)` | String | Gene | Gene symbol(s) reported as candidate causal genes by the original authors. |
| **14** | `MAPPED_GENE` | String | Gene | Gene(s) objectively mapped by Ensembl. Intergenic variants formatted as `GENE1 - GENE2`. |
| **15** | `UPSTREAM_GENE_ID` | Integer | Gene | Entrez Gene ID for nearest gene upstream of variant locus. |
| **16** | `DOWNSTREAM_GENE_ID` | Integer | Gene | Entrez Gene ID for nearest gene downstream of variant locus. |
| **17** | `SNP_GENE_IDS` | String | Gene | Entrez Gene ID(s) if variant falls within gene boundary (intronic, exonic, UTR). |
| **18** | `UPSTREAM_GENE_DISTANCE` | Integer | Genomic Position | Linear distance in base pairs from variant to nearest upstream gene. |
| **19** | `DOWNSTREAM_GENE_DISTANCE` | Integer | Genomic Position | Linear distance in base pairs from variant to nearest downstream gene. |
| **20** | `STRONGEST SNP-RISK ALLELE` | String | Variant / Mutation | Strongest associated SNP plus the reported risk/effect allele (e.g. `rs13143308-A`). |
| **21** | `SNPS` | String | Variant / Mutation | dbSNP rsID identifier(s). Can contain multiple variants for haplotypes or interactions. |
| **22** | `MERGED` | Binary (`0`/`1`) | Variant / Mutation | Denotes whether the original rsID has been merged into a subsequent record (`1` = yes). |
| **23** | `SNP_ID_CURRENT` | Integer | Variant / Mutation | Current active numeric rs identifier in dbSNP. |
| **24** | `CONTEXT` | String | Variant / Mutation | Predicted functional consequence (e.g. `intron_variant`, `missense_variant`, `intergenic_variant`). |
| **25** | `INTERGENIC` | Binary (`0`/`1`) | Genomic Position | Binary flag indicating if variant is intergenic (`1`) or intragenic (`0`). |
| **26** | `RISK ALLELE FREQUENCY` | Float / String | Allele Frequency | Reported risk/effect allele frequency (RAF) in controls. `'NR'` if not reported. |
| **27** | `P-VALUE` | Float / Scientific | Association Metric | Published p-value for the strongest association (rounded to 1 significant digit). |
| **28** | `PVALUE_MLOG` | Float | Association Metric | \(-\log_{10}(	ext{P-value})\). |
| **29** | `P-VALUE (TEXT)` | String | Association Metric | Contextual information for p-value (e.g. subgroup, smoking status, males only). |
| **30** | `OR or BETA` | Float | Association Metric | Reported odds ratio (OR) or beta effect coefficient. |
| **31** | `95% CI (TEXT)` | String | Association Metric | Published or estimated 95% confidence interval for effect size. |
| **32** | `PLATFORM [SNPS PASSING QC]` | String | Study Metadata | Genotyping array manufacturer and number of imputed/analyzed variants. |
| **33** | `CNV` | Binary (`Y`/`N`) | Study Metadata | Indicates whether association pertains to copy-number variation (`Y`/`N`). |
| **34** | `MAPPED_TRAIT` | String | Disease / Trait | Standardized ontology label from the Experimental Factor Ontology (EFO). |
| **35** | `MAPPED_TRAIT_URI` | URI | Disease / Trait | Permanent URI linking to the EFO term (e.g. `http://www.ebi.ac.uk/efo/EFO_0000305`). |
| **36** | `STUDY ACCESSION` | String | Study Metadata | Unique accession code for the GWAS Catalog study (e.g. `GCST005842`). |
| **37** | `GENOTYPING TECHNOLOGY` | String | Study Metadata | Technology employed (e.g. `Genome-wide genotyping array`, `Targeted array`). |

---

## 2. GWAS Catalog Ancestry Dataset (`gwas_catalog-ancestry`)

| Header | Description |
| :--- | :--- |
| `STUDY ACCESSION` | GWAS Catalog accession identifier matching the association records. |
| `PUBMEDID` | PubMed ID of the publication. |
| `FIRST AUTHOR` | First author name. |
| `DATE` | Publication date. |
| `INITIAL SAMPLE DESCRIPTION` | Free-text description of discovery sample composition. |
| `REPLICATION SAMPLE DESCRIPTION` | Free-text description of replication sample composition. |
| `STAGE` | Analysis phase: `initial` (discovery) or `replication`. |
| `NUMBER OF INDIVDUALS` | Sample size for the specific ancestry cohort. |
| `BROAD ANCESTRAL CATEGORY` | Standardized category label (e.g. `European`, `South Asian`, `East Asian`). |
| `COUNTRY OF ORIGIN` | Self-reported or geographic origin of individuals. |
| `COUNTRY OF RECRUITMENT` | Geographic location where cohort was recruited. |
| `ADDITIONAL ANCESTRY DESCRIPTION` | Detailed author notes or specific sub-population descriptors. |
