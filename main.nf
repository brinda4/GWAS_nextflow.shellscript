#!/usr/bin/env nextflow
/*
========================================================================================
   GWAS Analysis Pipeline - Nextflow DSL2 Workflow
========================================================================================
   Authors: GWAS Research Team
   Source Data: NHGRI-EBI GWAS Catalog
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def helpMessage() {
    log.info"""
    =================================================================================
    GWAS Analysis Nextflow Pipeline v${workflow.manifest.version}
    =================================================================================
    Usage:
      nextflow run main.nf [options]

    Input/Output Options:
      --raw_data_dir       Directory for raw downloaded catalog data [default: data/raw]
      --associations_zip   Path to gwas_catalog_associations.zip [default: data/raw/gwas_catalog_associations.zip]
      --ancestry_tsv       Path to gwas_catalog_ancestry.tsv [default: data/raw/gwas_catalog_ancestry.tsv]
      --outdir             Output directory for results [default: results]

    Filtering Parameters:
      --min_raf            Minimum Risk Allele Frequency (inclusive) [default: 0.005]
      --max_raf            Maximum Risk Allele Frequency (inclusive) [default: 0.020]
      --chunk_size         Streaming chunk size [default: 100000]

    Execution Profiles:
      -profile standard    Run locally
      -profile conda       Use Conda environment (environment.yml)
      -profile docker      Use Docker container
      -profile slurm       Execute on SLURM cluster
    =================================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

log.info """
=================================================================================
   GWAS Analysis Nextflow Pipeline v${workflow.manifest.version}
=================================================================================
Associations Zip : ${params.associations_zip}
Ancestry TSV     : ${params.ancestry_tsv}
RAF Window       : [${params.min_raf}, ${params.max_raf}] (0.5% - 2.0%)
Output Directory : ${params.outdir}
=================================================================================
"""

/*
----------------------------------------------------------------------------------------
   PROCESS 1: DOWNLOAD DATA
----------------------------------------------------------------------------------------
*/
process DOWNLOAD_DATA {
    tag "Download GWAS Catalog Data"
    publishDir "${params.raw_data_dir}", mode: 'copy'

    output:
    path "gwas_catalog_associations.zip", emit: assoc_zip
    path "gwas_catalog_ancestry.tsv",     emit: ancestry_tsv

    script:
    """
    python -c "from gwas_pipeline.downloader import download_file; download_file('${params.associations_url}', 'gwas_catalog_associations.zip'); download_file('${params.ancestry_url}', 'gwas_catalog_ancestry.tsv')"
    """
}

/*
----------------------------------------------------------------------------------------
   PROCESS 2: INSPECT SCHEMA & QUALITY AUDIT
----------------------------------------------------------------------------------------
*/
process INSPECT_CATALOG {
    tag "Inspect Schema & Missingness"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path assoc_zip
    path ancestry_tsv

    output:
    path "catalog_inspection_report.txt", emit: inspection_report

    script:
    """
    python -c "
from gwas_pipeline.inspector import inspect_catalog_columns
import json
results = inspect_catalog_columns('${assoc_zip}', '${ancestry_tsv}')
with open('catalog_inspection_report.txt', 'w') as f:
    f.write('=== GWAS CATALOG SCHEMA INSPECTION ===\\n')
    f.write(f'Associations File: {results[\"assoc_tsv_name\"]}\\n')
    f.write(f'Number of Columns: {len(results[\"assoc_columns\"])}\\n\\n')
    f.write('Columns and Missingness:\\n')
    for k, v in results[\"assoc_missingness\"].items():
        f.write(f'  {k}: {v*100:.2f}% null\\n')
"
    """
}

/*
----------------------------------------------------------------------------------------
   PROCESS 3: FILTER LOW-FREQUENCY VARIANTS (0.5% - 2.0% RAF)
----------------------------------------------------------------------------------------
*/
process FILTER_LOW_FREQUENCY {
    tag "Analysis 1: Low-Frequency Variants"
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path assoc_zip
    path ancestry_tsv

    output:
    path "low_freq_variants_eur_sas.tsv.gz",        emit: filtered_tsv
    path "low_freq_genes_and_variants_summary.tsv", emit: var_gene_summary
    path "analysis1_stats.json",                    emit: a1_stats

    script:
    """
    python -c "
from gwas_pipeline.filter_variants import filter_low_frequency_records
import json
stats = filter_low_frequency_records('${assoc_zip}', '${ancestry_tsv}', '.', chunk_size=${params.chunk_size})
with open('analysis1_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)
"
    """
}

/*
----------------------------------------------------------------------------------------
   PROCESS 4: MAP GENES, DISEASES, AND POPULATIONS
----------------------------------------------------------------------------------------
*/
process MAP_GENE_DISEASE {
    tag "Analysis 2: Gene-Disease-Population Mapping"
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path assoc_zip
    path ancestry_tsv

    output:
    path "gene_disease_population_mapping.tsv.gz", emit: triplets_tsv
    path "population_gene_classification.tsv",     emit: gene_class_tsv
    path "analysis2_stats.json",                   emit: a2_stats

    script:
    """
    python -c "
from gwas_pipeline.gene_disease_mapper import map_genes_to_diseases_and_populations
import json
stats = map_genes_to_diseases_and_populations('${assoc_zip}', '${ancestry_tsv}', '.', chunk_size=${params.chunk_size})
with open('analysis2_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)
"
    """
}

/*
----------------------------------------------------------------------------------------
   PROCESS 5: GENERATE PUBLICATION FIGURES
----------------------------------------------------------------------------------------
*/
process GENERATE_VISUALIZATIONS {
    tag "Generate Publication Plots"
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path filtered_tsv
    path triplets_tsv
    path gene_class_tsv

    output:
    path "*.png", emit: png_plots
    path "*.pdf", emit: pdf_plots

    script:
    """
    python -c "
from gwas_pipeline.visualizer import (
    plot_allele_frequency_distribution,
    plot_gene_disease_population_landscape,
    plot_population_gene_overlap
)
plot_allele_frequency_distribution('${filtered_tsv}', '.')
plot_gene_disease_population_landscape('${triplets_tsv}', '.')
plot_population_gene_overlap('${gene_class_tsv}', '.')
"
    """
}

/*
----------------------------------------------------------------------------------------
   PROCESS 6: SYNTHESIZE SUMMARY REPORTS
----------------------------------------------------------------------------------------
*/
process GENERATE_REPORT {
    tag "Generate Summary Reports"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path a1_stats
    path a2_stats

    output:
    path "gwas_executive_summary.csv",  emit: summary_csv
    path "final_workflow_summary.txt", emit: summary_txt

    script:
    """
    python -c "
import json
import pandas as pd

with open('${a1_stats}') as f:
    s1 = json.load(f)
with open('${a2_stats}') as f:
    s2 = json.load(f)

summary_rows = [
    {'Metric': 'Total Catalog Association Records', 'Value': s1['total_catalog_records'], 'Category': 'Overview'},
    {'Metric': 'Numeric RAF Records', 'Value': s1['raf_numeric_count'], 'Category': 'Data Quality'},
    {'Metric': 'Unreported RAF (NR)', 'Value': s1['raf_nr_count'], 'Category': 'Data Quality'},
    {'Metric': 'Low-Freq Filtered Records (0.5%-2.0% RAF)', 'Value': s1['filtered_records'], 'Category': 'Analysis 1'},
    {'Metric': 'Low-Freq Unique Variants (rsIDs)', 'Value': s1['unique_variants'], 'Category': 'Analysis 1'},
    {'Metric': 'Low-Freq Unique Genes', 'Value': s1['unique_genes'], 'Category': 'Analysis 1'},
    {'Metric': 'Total Unique Diseases/Traits', 'Value': s2['total_unique_diseases'], 'Category': 'Analysis 2'},
    {'Metric': 'Total Disease-Associated Genes', 'Value': s2['total_unique_genes'], 'Category': 'Analysis 2'},
    {'Metric': 'European-Only Genes', 'Value': s2['eur_only_genes'], 'Category': 'Analysis 2'},
    {'Metric': 'South Asian-Only Genes', 'Value': s2['sas_only_genes'], 'Category': 'Analysis 2'},
    {'Metric': 'Shared Genes (Both Populations)', 'Value': s2['both_pop_genes'], 'Category': 'Analysis 2'}
]
df = pd.DataFrame(summary_rows)
df.to_csv('gwas_executive_summary.csv', index=False)

with open('final_workflow_summary.txt', 'w') as out:
    out.write('=== GWAS NEXTFLOW PIPELINE EXECUTION SUMMARY ===\\n')
    for r in summary_rows:
        out.write(f\"{r['Category']} | {r['Metric']}: {r['Value']:,}\\n\")
"
    """
}

/*
----------------------------------------------------------------------------------------
   WORKFLOW EXECUTION ENTRYPOINT
----------------------------------------------------------------------------------------
*/
workflow {
    // Check if raw data files exist, otherwise trigger download
    if (file(params.associations_zip).exists() && file(params.ancestry_tsv).exists()) {
        ch_assoc_zip = Channel.fromPath(params.associations_zip)
        ch_ancestry   = Channel.fromPath(params.ancestry_tsv)
    } else {
        DOWNLOAD_DATA()
        ch_assoc_zip = DOWNLOAD_DATA.out.assoc_zip
        ch_ancestry   = DOWNLOAD_DATA.out.ancestry_tsv
    }

    // Step 1: Column inspection
    INSPECT_CATALOG(ch_assoc_zip, ch_ancestry)

    // Step 2: Analysis 1 (Low-frequency variant filtering)
    FILTER_LOW_FREQUENCY(ch_assoc_zip, ch_ancestry)

    // Step 3: Analysis 2 (Gene-disease-population mapping)
    MAP_GENE_DISEASE(ch_assoc_zip, ch_ancestry)

    // Step 4: Publication figures
    GENERATE_VISUALIZATIONS(
        FILTER_LOW_FREQUENCY.out.filtered_tsv,
        MAP_GENE_DISEASE.out.triplets_tsv,
        MAP_GENE_DISEASE.out.gene_class_tsv
    )

    // Step 5: Summary report
    GENERATE_REPORT(
        FILTER_LOW_FREQUENCY.out.a1_stats,
        MAP_GENE_DISEASE.out.a2_stats
    )
}

workflow.onComplete {
    log.info """
    =================================================================================
    Pipeline completed at: ${workflow.complete}
    Execution duration   : ${workflow.duration}
    Execution status     : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Output directory     : ${params.outdir}
    =================================================================================
    """
}
