"""
Low-Frequency Variant & Gene Filtering Module (Analysis 1)
---------------------------------------------------------
Filters the complete GWAS Catalog for variants with Risk Allele Frequency (RAF)
between 0.5% and 2.0% inclusive (0.005 <= RAF <= 0.020) in European and South Asian cohorts.
"""

import os
import re
import zipfile
import gzip
import pandas as pd
import numpy as np
import logging
from typing import Dict, Set, Tuple, List, Any

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def build_ancestry_map(ancestry_tsv_path: str) -> Dict[str, Set[str]]:
    """
    Builds a lookup mapping study accession to broad ancestral categories.
    Prioritizes discovery/initial stage cohorts.
    """
    logger.info(f"Loading ancestry metadata from {ancestry_tsv_path}...")
    df_anc = pd.read_csv(ancestry_tsv_path, sep="\t", low_memory=False)
    
    study_anc_map = {}
    for acc, grp in df_anc.groupby("STUDY ACCESSION"):
        initial_grp = grp[grp["STAGE"] == "initial"]
        cats = set()
        for cat_str in initial_grp["BROAD ANCESTRAL CATEGORY"].dropna():
            for part in str(cat_str).split(","):
                p = part.strip()
                if p:
                    cats.add(p)
        if not cats:  # fallback to all stages if initial is unannotated
            for cat_str in grp["BROAD ANCESTRAL CATEGORY"].dropna():
                for part in str(cat_str).split(","):
                    p = part.strip()
                    if p:
                        cats.add(p)
        study_anc_map[acc] = cats
        
    logger.info(f"Successfully mapped ancestry for {len(study_anc_map):,} studies.")
    return study_anc_map


def parse_numeric_raf(val: Any) -> float:
    """
    Extracts a numeric float from the raw RISK ALLELE FREQUENCY field.
    Handles floats, scientific notation, and ignores 'NR' or text annotations.
    """
    if pd.isna(val):
        return np.nan
    s = str(val).strip()
    if s.upper() == "NR" or not s:
        return np.nan
    match = re.search(r"[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[eE][-+]?\d+)?", s)
    if match:
        try:
            return float(match.group(0))
        except ValueError:
            return np.nan
    return np.nan


def tokenize_genes(gene_str: Any) -> List[str]:
    """
    Tokenizes mapped or reported gene strings into discrete, clean gene symbols.
    Handles 'GENE1 - GENE2' (intergenic) and comma-separated transcripts.
    """
    if pd.isna(gene_str) or not str(gene_str).strip() or str(gene_str).strip().upper() == "NR":
        return []
    s = str(gene_str).strip()
    tokens = re.split(r"\s+-\s+|[,;]\s*", s)
    genes = [t.strip() for t in tokens if t.strip() and t.strip().upper() != "NR"]
    return genes


def tokenize_variants(snp_str: Any) -> List[str]:
    """
    Tokenizes variant identifiers (rsIDs or genomic variants).
    """
    if pd.isna(snp_str) or not str(snp_str).strip() or str(snp_str).strip().upper() == "NR":
        return []
    s = str(snp_str).strip()
    tokens = re.split(r"[,;x\s]+", s)
    variants = [t.strip() for t in tokens if t.strip().startswith("rs") or t.strip().isdigit()]
    if not variants and s and s.upper() != "NR":
        variants = [s]
    return variants


def filter_low_frequency_records(associations_zip_path: str, ancestry_tsv_path: str, output_dir: str, chunk_size: int = 100000, study_anc_map: Dict[str, Set[str]] = None) -> Dict[str, Any]:
    """
    Streams the complete associations table, filters for RAF in [0.005, 0.020] for EUR and SAS,
    and writes out processed result tables.
    """
    os.makedirs(output_dir, exist_ok=True)
    if study_anc_map is None:
        study_anc_map = build_ancestry_map(ancestry_tsv_path)
    
    logger.info("Streaming full associations dataset and filtering low-frequency variants...")
    total_records = 0
    raf_numeric_count = 0
    raf_nr_count = 0
    raf_na_count = 0
    
    filtered_chunks = []
    
    with zipfile.ZipFile(associations_zip_path, "r") as z:
        assoc_tsv_name = z.infolist()[0].filename
        with z.open(assoc_tsv_name) as f:
            for chunk in pd.read_csv(f, sep="\t", chunksize=chunk_size, low_memory=False):
                total_records += len(chunk)
                
                # Map populations
                pops = chunk["STUDY ACCESSION"].map(study_anc_map)
                is_eur = pops.apply(lambda s: isinstance(s, set) and "European" in s)
                is_sas = pops.apply(lambda s: isinstance(s, set) and "South Asian" in s)
                
                # Audit RAF
                raw_raf = chunk["RISK ALLELE FREQUENCY"]
                raf_na_count += raw_raf.isna().sum()
                raf_nr_count += (raw_raf.astype(str).str.strip().str.upper() == "NR").sum()
                
                num_raf = raw_raf.apply(parse_numeric_raf)
                raf_numeric_count += num_raf.notna().sum()
                
                # Target filter: RAF in [0.005, 0.020] AND (EUR or SAS)
                in_raf_range = num_raf.notna() & (num_raf >= 0.005) & (num_raf <= 0.020)
                in_target_pop = is_eur | is_sas
                
                match_mask = in_raf_range & in_target_pop
                if match_mask.any():
                    sub = chunk[match_mask].copy()
                    sub["PARSED_RAF"] = num_raf[match_mask]
                    sub["IS_EUR"] = is_eur[match_mask]
                    sub["IS_SAS"] = is_sas[match_mask]
                    sub["POPULATION_GROUP"] = np.where(
                        sub["IS_EUR"] & sub["IS_SAS"], "Both EUR & SAS",
                        np.where(sub["IS_EUR"], "European Only", "South Asian Only")
                    )
                    filtered_chunks.append(sub[[
                        "STUDY ACCESSION", "PUBMEDID", "DISEASE/TRAIT", "MAPPED_TRAIT",
                        "CHR_ID", "CHR_POS", "SNPS", "STRONGEST SNP-RISK ALLELE",
                        "MAPPED_GENE", "REPORTED GENE(S)", "PARSED_RAF",
                        "IS_EUR", "IS_SAS", "POPULATION_GROUP"
                    ]])
                    
                if total_records % 500000 == 0:
                    logger.info(f"Processed {total_records:,} records...")

    logger.info(f"Finished processing all {total_records:,} catalog records.")
    df_filtered = pd.concat(filtered_chunks, ignore_index=True)
    
    # Save processed filtered dataset (compressed TSV)
    out_filtered_path = os.path.join(output_dir, "low_freq_variants_eur_sas.tsv.gz")
    logger.info(f"Writing {len(df_filtered):,} filtered records to {out_filtered_path}...")
    df_filtered.to_csv(out_filtered_path, sep="\t", index=False, compression="gzip")
    
    # Extract unique low-frequency variants and genes
    unique_variants = set()
    unique_genes = set()
    variant_gene_records = []
    
    for _, row in df_filtered.iterrows():
        snps = tokenize_variants(row["SNPS"])
        genes = tokenize_genes(row["MAPPED_GENE"]) if pd.notna(row["MAPPED_GENE"]) else tokenize_genes(row["REPORTED GENE(S)"])
        pop_grp = row["POPULATION_GROUP"]
        raf = row["PARSED_RAF"]
        trait = row["DISEASE/TRAIT"]
        
        for s in snps:
            unique_variants.add(s)
        for g in genes:
            unique_genes.add(g)
            
        for s in (snps if snps else [None]):
            for g in (genes if genes else [None]):
                variant_gene_records.append({
                    "VARIANT": s,
                    "GENE": g,
                    "POPULATION_GROUP": pop_grp,
                    "RAF": raf,
                    "DISEASE_TRAIT": trait
                })
                
    df_var_gene = pd.DataFrame(variant_gene_records).drop_duplicates()
    out_summary_path = os.path.join(output_dir, "low_freq_genes_and_variants_summary.tsv")
    df_var_gene.to_csv(out_summary_path, sep="\t", index=False)
    
    stats = {
        "total_catalog_records": int(total_records),
        "raf_numeric_count": int(raf_numeric_count),
        "raf_nr_count": int(raf_nr_count),
        "raf_na_count": int(raf_na_count),
        "filtered_records": int(len(df_filtered)),
        "eur_records": int(df_filtered["IS_EUR"].sum()),
        "sas_records": int(df_filtered["IS_SAS"].sum()),
        "both_records": int((df_filtered["IS_EUR"] & df_filtered["IS_SAS"]).sum()),
        "unique_variants": int(len(unique_variants)),
        "unique_genes": int(len(unique_genes))
    }
    
    logger.info(f"Low-frequency filtering complete: {stats['filtered_records']:,} records, {stats['unique_variants']:,} variants, {stats['unique_genes']:,} genes.")
    return stats
