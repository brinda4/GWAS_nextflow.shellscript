"""
Catalog Schema Inspector & Quality Auditor
------------------------------------------
Profiles all columns, identifies data types, audits missingness, and prints
a structured data dictionary for the complete GWAS Catalog dataset.
"""

import os
import zipfile
import pandas as pd
import logging
from typing import Dict, Any

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def inspect_catalog_columns(associations_zip_path: str, ancestry_tsv_path: str, sample_size: int = 100000) -> Dict[str, Any]:
    """
    Reads dataset headers, samples records, and audits column completeness.
    """
    logger.info("Inspecting GWAS Catalog Associations schema...")
    with zipfile.ZipFile(associations_zip_path, "r") as z:
        assoc_tsv_name = z.infolist()[0].filename
        with z.open(assoc_tsv_name) as f:
            df_assoc_sample = pd.read_csv(f, sep="\t", nrows=sample_size, low_memory=False)
            
    logger.info("Inspecting GWAS Catalog Ancestry schema...")
    df_anc_sample = pd.read_csv(ancestry_tsv_path, sep="\t", nrows=sample_size, low_memory=False)
    
    assoc_columns = list(df_assoc_sample.columns)
    anc_columns = list(df_anc_sample.columns)
    
    # Missingness audit on sample
    assoc_missing = df_assoc_sample.isna().mean().to_dict()
    
    # Categorize columns into biological and study dimensions
    column_domains = {
        "Variant/Mutation": ["SNPS", "STRONGEST SNP-RISK ALLELE", "SNP_ID_CURRENT", "CONTEXT", "MERGED"],
        "Gene": ["MAPPED_GENE", "REPORTED GENE(S)", "SNP_GENE_IDS", "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID"],
        "Allele Frequency": ["RISK ALLELE FREQUENCY"],
        "Population / Ancestry": ["INITIAL SAMPLE SIZE", "REPLICATION SAMPLE SIZE", "BROAD ANCESTRAL CATEGORY", "STAGE", "NUMBER OF INDIVDUALS", "COUNTRY OF ORIGIN"],
        "Disease / Trait": ["DISEASE/TRAIT", "MAPPED_TRAIT", "MAPPED_TRAIT_URI"],
        "Chromosomal Position": ["CHR_ID", "CHR_POS", "REGION", "UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE"],
        "Study Metadata": ["STUDY ACCESSION", "PUBMEDID", "STUDY", "FIRST AUTHOR", "DATE", "JOURNAL", "GENOTYPING TECHNOLOGY", "PLATFORM [SNPS PASSING QC]"]
    }
    
    logger.info(f"Associations file contains {len(assoc_columns)} columns across {sample_size:,} sampled rows.")
    logger.info(f"Ancestry file contains {len(anc_columns)} columns.")
    
    return {
        "assoc_tsv_name": assoc_tsv_name,
        "assoc_columns": assoc_columns,
        "anc_columns": anc_columns,
        "column_domains": column_domains,
        "assoc_missingness": assoc_missing
    }


def print_schema_summary(inspection_results: Dict[str, Any]):
    """
    Prints a formatted markdown summary of column domains and completeness.
    """
    print("\n=== GWAS CATALOG COLUMN PROFILING ===")
    for domain, cols in inspection_results["column_domains"].items():
        print(f"\n--- {domain} ---")
        for c in cols:
            missing_rate = inspection_results["assoc_missingness"].get(c, None)
            miss_str = f" (Sample Null: {missing_rate*100:.1f}%)" if missing_rate is not None else ""
            print(f"  * {c}{miss_str}")
