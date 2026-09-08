"""
Data Downloader Module
----------------------
Handles automated retrieval of official NHGRI-EBI GWAS Catalog data releases.
"""

import os
import time
import urllib.request
import zipfile
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Official NHGRI-EBI download endpoints
ASSOCIATIONS_URL = "https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false"
ANCESTRY_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv"


def download_file(url: str, output_path: str, chunk_size: int = 1024 * 1024) -> str:
    """
    Downloads a remote file with progress tracking.
    """
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        logger.info(f"Target file already exists: {output_path} ({os.path.getsize(output_path)/(1024*1024):.1f} MB). Skipping download.")
        return output_path

    logger.info(f"Connecting to {url}...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (GWAS-Pipeline/1.0)"})
    t0 = time.time()
    
    with urllib.request.urlopen(req) as resp, open(output_path, "wb") as out_f:
        total_header = resp.headers.get("Content-Length")
        total_size = int(total_header) if total_header else None
        downloaded = 0
        
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            out_f.write(chunk)
            downloaded += len(chunk)
            if total_size:
                pct = (downloaded / total_size) * 100
                logger.info(f"Downloaded {downloaded/(1024*1024):.1f} MB / {total_size/(1024*1024):.1f} MB ({pct:.1f}%)")
            else:
                if downloaded % (10 * 1024 * 1024) < chunk_size:
                    logger.info(f"Downloaded {downloaded/(1024*1024):.1f} MB...")
                    
    logger.info(f"Successfully downloaded {output_path} ({downloaded/(1024*1024):.1f} MB) in {time.time()-t0:.1f}s.")
    return output_path


def download_gwas_catalog_data(raw_data_dir: str) -> dict:
    """
    Ensures both the associations zip and ancestry tsv datasets are downloaded.
    """
    os.makedirs(raw_data_dir, exist_ok=True)
    
    assoc_zip = os.path.join(raw_data_dir, "gwas_catalog_associations.zip")
    ancestry_tsv = os.path.join(raw_data_dir, "gwas_catalog_ancestry.tsv")
    
    download_file(ASSOCIATIONS_URL, assoc_zip)
    download_file(ANCESTRY_URL, ancestry_tsv)
    
    return {
        "associations_zip": assoc_zip,
        "ancestry_tsv": ancestry_tsv
    }
