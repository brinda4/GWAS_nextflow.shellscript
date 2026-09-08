# Raw Data Directory

Raw datasets downloaded from the NHGRI-EBI GWAS Catalog are stored here.
To download the raw files automatically, run:

```bash
python run_pipeline.py --download
# or
./run_pipeline.sh
```

### Source URLs:
- Associations: https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false
- Ancestry: https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv
