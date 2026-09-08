"""
Gene-Disease-Population Landscape Mapper (Analysis 2)
-----------------------------------------------------
Maps unique genes to diseases and complex traits across European, South Asian,
and multi-ancestry cohorts. Distinguishes European-only, South Asian-only,
and shared disease-associated genes.
"""

import os
import zipfile
import pandas as pd
import numpy as np
import logging
from typing import Dict, Any, List, Set
from gwas_pipeline.filter_variants import build_ancestry_map, tokenize_genes

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def map_genes_to_diseases_and_populations(associations_zip_path: str, ancestry_tsv_path: str, output_dir: str, chunk_size: int = 100000, study_anc_map: Dict[str, Set[str]] = None) -> Dict[str, Any]:
    """
    Streams associations, extracts gene-disease associations in EUR and SAS cohorts,
    and categorizes genes by ancestral exclusivity.
    """
    os.makedirs(output_dir, exist_ok=True)
    if study_anc_map is None:
        study_anc_map = build_ancestry_map(ancestry_tsv_path)
    
    logger.info("Streaming full catalog for Gene-Disease-Population mapping...")
    gene_disease_pop_tuples = []
    
    with zipfile.ZipFile(associations_zip_path, "r") as z:
        assoc_tsv_name = z.infolist()[0].filename
        with z.open(assoc_tsv_name) as f:
            for chunk in pd.read_csv(f, sep="\t", chunksize=chunk_size, low_memory=False,
                                     usecols=["STUDY ACCESSION", "MAPPED_GENE", "REPORTED GENE(S)", "DISEASE/TRAIT"]):
                
                # Check population
                pops = chunk["STUDY ACCESSION"].map(study_anc_map)
                is_eur = pops.apply(lambda s: isinstance(s, set) and "European" in s)
                is_sas = pops.apply(lambda s: isinstance(s, set) and "South Asian" in s)
                
                valid_mask = (is_eur | is_sas) & chunk["DISEASE/TRAIT"].notna()
                if valid_mask.any():
                    sub = chunk[valid_mask]
                    eur_arr = is_eur[valid_mask].values
                    sas_arr = is_sas[valid_mask].values
                    gene_col = sub["MAPPED_GENE"].fillna(sub["REPORTED GENE(S)"]).values
                    trait_col = sub["DISEASE/TRAIT"].values
                    
                    for g_str, trait, h_eur, h_sas in zip(gene_col, trait_col, eur_arr, sas_arr):
                        g_list = tokenize_genes(g_str)
                        if g_list and trait:
                            t_clean = str(trait).strip()
                            for g in g_list:
                                gene_disease_pop_tuples.append((g, t_clean, bool(h_eur), bool(h_sas)))

    logger.info(f"Extracted {len(gene_disease_pop_tuples):,} gene-disease association occurrences.")
    df_gdp = pd.DataFrame(gene_disease_pop_tuples, columns=["GENE", "DISEASE_TRAIT", "HAS_EUR", "HAS_SAS"])
    
    # 1. Triplet aggregated table (Gene, Disease/Trait, Population Category, Association Count)
    logger.info("Aggregating gene-disease-population triplets...")
    df_triplets = df_gdp.groupby(["GENE", "DISEASE_TRAIT"]).agg({
        "HAS_EUR": "any",
        "HAS_SAS": "any",
        "GENE": "count"
    }).rename(columns={"GENE": "ASSOCIATION_COUNT"}).reset_index()
    
    df_triplets["POPULATION_CONTEXT"] = np.where(
        df_triplets["HAS_EUR"] & df_triplets["HAS_SAS"], "Both EUR & SAS",
        np.where(df_triplets["HAS_EUR"], "European Only", "South Asian Only")
    )
    
    triplet_out_path = os.path.join(output_dir, "gene_disease_population_mapping.tsv.gz")
    logger.info(f"Saving {len(df_triplets):,} triplets to {triplet_out_path}...")
    df_triplets.to_csv(triplet_out_path, sep="\t", index=False, compression="gzip")
    
    # 2. Gene-level classification across all diseases
    logger.info("Classifying genes by overall ancestral involvement...")
    df_gene_class = df_gdp.groupby("GENE").agg({
        "HAS_EUR": "any",
        "HAS_SAS": "any",
        "DISEASE_TRAIT": ["nunique", lambda s: "; ".join(pd.Series(s).value_counts().head(3).index)]
    }).reset_index()
    
    df_gene_class.columns = ["GENE", "HAS_EUR", "HAS_SAS", "TOTAL_DISEASES", "TOP_ASSOCIATED_TRAITS"]
    df_gene_class["POPULATION_CLASSIFICATION"] = np.where(
        df_gene_class["HAS_EUR"] & df_gene_class["HAS_SAS"], "Both European & South Asian",
        np.where(df_gene_class["HAS_EUR"], "European Only", "South Asian Only")
    )
    
    gene_class_out_path = os.path.join(output_dir, "population_gene_classification.tsv")
    df_gene_class.to_csv(gene_class_out_path, sep="\t", index=False)
    
    total_genes = len(df_gene_class)
    eur_only = (df_gene_class["POPULATION_CLASSIFICATION"] == "European Only").sum()
    sas_only = (df_gene_class["POPULATION_CLASSIFICATION"] == "South Asian Only").sum()
    both_pop = (df_gene_class["POPULATION_CLASSIFICATION"] == "Both European & South Asian").sum()
    total_diseases = df_gdp["DISEASE_TRAIT"].nunique()
    
    stats = {
        "total_unique_diseases": int(total_diseases),
        "total_unique_genes": int(total_genes),
        "eur_only_genes": int(eur_only),
        "sas_only_genes": int(sas_only),
        "both_pop_genes": int(both_pop)
    }
    
    logger.info(f"Gene-disease mapping complete: {total_genes:,} genes ({eur_only:,} EUR-only, {sas_only:,} SAS-only, {both_pop:,} shared).")
    return stats
