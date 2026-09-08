"""
Master CLI Runner for GWAS Analysis Pipeline
============================================
Executes downloading, column inspection, low-frequency variant filtering (Analysis 1),
gene-disease-population mapping (Analysis 2), and publication visualization.
"""

import os
import sys
import argparse
import logging
import json

from gwas_pipeline.downloader import download_gwas_catalog_data
from gwas_pipeline.inspector import inspect_catalog_columns, print_schema_summary
from gwas_pipeline.filter_variants import filter_low_frequency_records, build_ancestry_map
from gwas_pipeline.gene_disease_mapper import map_genes_to_diseases_and_populations
from gwas_pipeline.visualizer import (
    plot_allele_frequency_distribution,
    plot_gene_disease_population_landscape,
    plot_population_gene_overlap
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("GWAS-Runner")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DATA_DIR = os.path.join(BASE_DIR, "data", "raw")
PROCESSED_DATA_DIR = os.path.join(BASE_DIR, "data", "processed")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")


def main():
    parser = argparse.ArgumentParser(description="Execute complete NHGRI-EBI GWAS Catalog analysis pipeline.")
    parser.add_argument("--download", action="store_true", help="Download raw GWAS Catalog files if missing.")
    parser.add_argument("--inspect", action="store_true", help="Inspect and document catalog columns and missingness.")
    parser.add_argument("--filter", action="store_true", help="Execute Analysis 1 (RAF 0.5%-2.0% in EUR & SAS).")
    parser.add_argument("--map-genes", action="store_true", help="Execute Analysis 2 (Gene-Disease-Population mapping).")
    parser.add_argument("--plot", action="store_true", help="Generate publication-quality figures.")
    parser.add_argument("--all", action="store_true", help="Execute all pipeline stages sequentially.")
    
    args = parser.parse_args()
    
    # Default to --all if no specific stage requested
    if not (args.download or args.inspect or args.filter or args.map_genes or args.plot or args.all):
        args.all = True
        
    assoc_zip = os.path.join(RAW_DATA_DIR, "gwas_catalog_associations.zip")
    ancestry_tsv = os.path.join(RAW_DATA_DIR, "gwas_catalog_ancestry.tsv")
    
    if args.download or args.all:
        logger.info("=== STEP 1: DOWNLOADING DATA ===")
        download_gwas_catalog_data(RAW_DATA_DIR)
        
    if args.inspect or args.all:
        logger.info("=== STEP 2: INSPECTING CATALOG COLUMNS & SCHEMA ===")
        results = inspect_catalog_columns(assoc_zip, ancestry_tsv)
        print_schema_summary(results)
        
    study_anc_map = None
    if args.filter or args.map_genes or args.all:
        study_anc_map = build_ancestry_map(ancestry_tsv)

    if args.filter or args.all:
        logger.info("=== STEP 3: ANALYSIS 1 - LOW-FREQUENCY FILTERING (0.5% - 2.0% RAF) ===")
        a1_stats = filter_low_frequency_records(assoc_zip, ancestry_tsv, PROCESSED_DATA_DIR, study_anc_map=study_anc_map)
        logger.info(f"Analysis 1 Results:\n{json.dumps(a1_stats, indent=2, default=int)}")
        
    if args.map_genes or args.all:
        logger.info("=== STEP 4: ANALYSIS 2 - GENE-DISEASE-POPULATION MAPPING ===")
        a2_stats = map_genes_to_diseases_and_populations(assoc_zip, ancestry_tsv, PROCESSED_DATA_DIR, study_anc_map=study_anc_map)
        logger.info(f"Analysis 2 Results:\n{json.dumps(a2_stats, indent=2, default=int)}")
        
    if args.plot or args.all:
        logger.info("=== STEP 5: GENERATING PUBLICATION FIGURES ===")
        filtered_path = os.path.join(PROCESSED_DATA_DIR, "low_freq_variants_eur_sas.tsv.gz")
        triplets_path = os.path.join(PROCESSED_DATA_DIR, "gene_disease_population_mapping.tsv.gz")
        gene_class_path = os.path.join(PROCESSED_DATA_DIR, "population_gene_classification.tsv")
        
        plot_allele_frequency_distribution(filtered_path, FIGURES_DIR)
        plot_gene_disease_population_landscape(triplets_path, FIGURES_DIR)
        plot_population_gene_overlap(gene_class_path, FIGURES_DIR)
        
    logger.info("=== GWAS PIPELINE COMPLETE ===")


if __name__ == "__main__":
    main()
