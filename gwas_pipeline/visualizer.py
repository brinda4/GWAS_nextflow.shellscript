"""
Publication-Quality Visualizer
------------------------------
Generates high-resolution (300 DPI) publication-grade figures in PNG and vector PDF:
1. Allele-frequency distribution and comparison between European and South Asian cohorts.
2. Gene-disease-population landscape (stacked bar plot of top traits & ancestral contribution).
3. Population gene overlap summary.
"""

import os
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Configure publication typography and theme
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
plt.rcParams["axes.edgecolor"] = "#2b2b2b"
plt.rcParams["axes.linewidth"] = 1.0


def plot_allele_frequency_distribution(filtered_tsv_gz_path: str, output_dir: str):
    """
    Figure 1: Allele-frequency distribution comparison between European and South Asian populations.
    Panel A: Kernel Density Estimates (KDE) and Histograms.
    Panel B: Boxplot & Violin Plot comparing European vs South Asian RAF.
    """
    logger.info("Generating Figure 1: Allele Frequency Distribution...")
    df = pd.read_csv(filtered_tsv_gz_path, sep="\t", compression="gzip")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    # Palette
    palette = {
        "European Only": "#2b5c8f",
        "South Asian Only": "#d95f02",
        "Both EUR & SAS": "#7570b3"
    }
    
    # Panel A: Distribution plot (KDE + Histogram)
    ax0 = axes[0]
    sns.histplot(
        data=df[df["POPULATION_GROUP"].isin(["European Only", "South Asian Only", "Both EUR & SAS"])],
        x="PARSED_RAF",
        hue="POPULATION_GROUP",
        palette=palette,
        bins=30,
        kde=True,
        stat="density",
        common_norm=False,
        element="step",
        alpha=0.35,
        ax=ax0
    )
    ax0.set_title("A. Allele Frequency Density (0.5% - 2.0% RAF)", fontsize=13, fontweight="bold", pad=12)
    ax0.set_xlabel("Risk Allele Frequency (RAF)", fontsize=11, fontweight="bold")
    ax0.set_ylabel("Density", fontsize=11, fontweight="bold")
    ax0.set_xlim(0.0045, 0.0205)
    ax0.grid(True, linestyle="--", alpha=0.5)
    
    # Annotate summary stats
    eur_med = df[df["IS_EUR"]]["PARSED_RAF"].median()
    sas_med = df[df["IS_SAS"]]["PARSED_RAF"].median()
    ax0.axvline(eur_med, color="#2b5c8f", linestyle=":", linewidth=1.5, label=f"EUR Median: {eur_med:.4f}")
    ax0.axvline(sas_med, color="#d95f02", linestyle=":", linewidth=1.5, label=f"SAS Median: {sas_med:.4f}")
    ax0.legend(title="Population Cohort", frameon=True, facecolor="white", edgecolor="#cccccc")
    
    # Panel B: Box / Violin plot
    ax1 = axes[1]
    order = ["European Only", "Both EUR & SAS", "South Asian Only"]
    sns.boxplot(
        data=df,
        x="POPULATION_GROUP",
        y="PARSED_RAF",
        order=order,
        palette=palette,
        width=0.45,
        boxprops=dict(alpha=0.8),
        ax=ax1
    )
    # Overlay swarm/strip summary
    sns.stripplot(
        data=df.sample(min(len(df), 2000), random_state=42),
        x="POPULATION_GROUP",
        y="PARSED_RAF",
        order=order,
        palette=palette,
        size=2,
        jitter=0.2,
        alpha=0.4,
        ax=ax1
    )
    ax1.set_title("B. Low-Frequency Variation by Cohort Type", fontsize=13, fontweight="bold", pad=12)
    ax1.set_xlabel("Population Group", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Risk Allele Frequency (RAF)", fontsize=11, fontweight="bold")
    ax1.set_ylim(0.0045, 0.0205)
    ax1.grid(True, linestyle="--", alpha=0.5, axis="y")
    
    # Sample counts on x-axis
    counts = df["POPULATION_GROUP"].value_counts()
    labels = [f"{grp}\n(n={counts.get(grp, 0):,})" for grp in order]
    ax1.set_xticklabels(labels, fontsize=10)
    
    plt.tight_layout()
    
    png_path = os.path.join(output_dir, "fig1_allele_frequency_distribution.png")
    pdf_path = os.path.join(output_dir, "fig1_allele_frequency_distribution.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    logger.info(f"Figure 1 saved to {png_path} and {pdf_path}.")


def plot_gene_disease_population_landscape(triplets_tsv_gz_path: str, output_dir: str, top_n_traits: int = 20):
    """
    Figure 2: Relationship between Genes, Diseases/Traits, and Populations.
    Panel A: Stacked bar plot of gene discovery across top traits by population.
    Panel B: Heatmap matrix of ancestral gene proportion across high-impact disease domains.
    """
    logger.info("Generating Figure 2: Gene-Disease-Population Landscape...")
    df = pd.read_csv(triplets_tsv_gz_path, sep="\t", compression="gzip")
    
    # Top traits by total distinct genes
    trait_gene_counts = df.groupby("DISEASE_TRAIT")["GENE"].nunique().sort_values(ascending=False)
    top_traits = trait_gene_counts.head(top_n_traits).index.tolist()
    
    # Filter for top traits
    df_top = df[df["DISEASE_TRAIT"].isin(top_traits)]
    
    # Cross-tabulate: Trait x Population Context (distinct genes)
    trait_pop_genes = df_top.groupby(["DISEASE_TRAIT", "POPULATION_CONTEXT"])["GENE"].nunique().unstack(fill_value=0)
    
    # Ensure all columns present
    for col in ["European Only", "Both EUR & SAS", "South Asian Only"]:
        if col not in trait_pop_genes.columns:
            trait_pop_genes[col] = 0
            
    trait_pop_genes = trait_pop_genes[["European Only", "Both EUR & SAS", "South Asian Only"]]
    trait_pop_genes = trait_pop_genes.loc[top_traits]  # preserve sorted order
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), dpi=300, gridspec_kw={"width_ratios": [1.4, 1.0]})
    
    # Panel A: Horizontal Stacked Bar Chart
    ax0 = axes[0]
    palette = ["#2b5c8f", "#7570b3", "#d95f02"]
    
    # Plot bars
    trait_pop_genes.plot(
        kind="barh",
        stacked=True,
        color=palette,
        ax=ax0,
        width=0.75,
        edgecolor="#222222",
        linewidth=0.5
    )
    ax0.invert_yaxis()  # top trait at the top
    ax0.set_title(f"A. Gene Discovery Breakdown for Top {top_n_traits} Complex Traits", fontsize=13, fontweight="bold", pad=12)
    ax0.set_xlabel("Number of Associated Genes", fontsize=11, fontweight="bold")
    ax0.set_ylabel("Disease or Trait", fontsize=11, fontweight="bold")
    ax0.legend(title="Population Context", frameon=True, facecolor="white", edgecolor="#cccccc")
    ax0.grid(True, linestyle="--", alpha=0.5, axis="x")
    
    # Panel B: Normalized Proportions Heatmap
    ax1 = axes[1]
    trait_pop_pct = trait_pop_genes.div(trait_pop_genes.sum(axis=1), axis=0) * 100
    
    sns.heatmap(
        trait_pop_pct,
        cmap="Blues",
        annot=True,
        fmt=".1f",
        cbar_kws={"label": "Proportion of Discovered Genes (%)"},
        linewidths=0.5,
        linecolor="#dddddd",
        ax=ax1
    )
    ax1.set_title("B. Population Discovery Proportions (%)", fontsize=13, fontweight="bold", pad=12)
    ax1.set_xlabel("Population Group", fontsize=11, fontweight="bold")
    ax1.set_ylabel("")
    ax1.set_yticklabels([])  # aligned with panel A
    
    plt.tight_layout()
    
    png_path = os.path.join(output_dir, "fig2_gene_disease_population_landscape.png")
    pdf_path = os.path.join(output_dir, "fig2_gene_disease_population_landscape.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    logger.info(f"Figure 2 saved to {png_path} and {pdf_path}.")


def plot_population_gene_overlap(gene_class_tsv_path: str, output_dir: str):
    """
    Figure 3: Global Ancestral Specificity and Sharing of 30,054 Disease-Linked Genes.
    """
    logger.info("Generating Figure 3: Population Gene Overlap...")
    df = pd.read_csv(gene_class_tsv_path, sep="\t")
    
    class_counts = df["POPULATION_CLASSIFICATION"].value_counts()
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    colors = ["#7570b3", "#2b5c8f", "#d95f02"]
    
    wedges, texts, autotexts = ax.pie(
        class_counts.values,
        labels=[f"{k}\n(n={v:,})" for k, v in class_counts.items()],
        autopct="%1.1f%%",
        startangle=140,
        colors=colors,
        explode=(0.04, 0.04, 0.12),
        wedgeprops=dict(width=0.6, edgecolor="#ffffff", linewidth=2)
    )
    
    plt.setp(autotexts, size=11, weight="bold", color="white")
    plt.setp(texts, size=10, weight="bold")
    ax.set_title(f"Global Ancestral Sharing of Disease-Linked Genes\n(Total = {len(df):,} Unique Genes)", fontsize=13, fontweight="bold", pad=14)
    
    png_path = os.path.join(output_dir, "fig3_population_gene_overlap.png")
    pdf_path = os.path.join(output_dir, "fig3_population_gene_overlap.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    logger.info(f"Figure 3 saved to {png_path} and {pdf_path}.")
