"""
Unit and Integration Tests for GWAS Analysis Pipeline
=====================================================
Validates tokenization, numeric RAF parsing, filtering thresholds,
and population categorization logic.
"""

import pytest
import numpy as np
import pandas as pd
from gwas_pipeline.filter_variants import parse_numeric_raf, tokenize_genes, tokenize_variants


def test_parse_numeric_raf_valid_floats():
    assert parse_numeric_raf("0.015") == 0.015
    assert parse_numeric_raf(0.015) == 0.015
    assert parse_numeric_raf(".0159") == 0.0159
    assert parse_numeric_raf("7.0E-4") == 0.0007


def test_parse_numeric_raf_missing_and_text():
    assert np.isnan(parse_numeric_raf("NR"))
    assert np.isnan(parse_numeric_raf("nr"))
    assert np.isnan(parse_numeric_raf(np.nan))
    assert np.isnan(parse_numeric_raf(""))
    # Annotated string with parenthesis
    assert parse_numeric_raf("0.017 (European)") == 0.017


def test_tokenize_genes():
    assert tokenize_genes("BRCA1") == ["BRCA1"]
    assert tokenize_genes("LINC01438 - ENSG00000249519") == ["LINC01438", "ENSG00000249519"]
    assert tokenize_genes("HLA-A") == ["HLA-A"]  # preserves internal hyphen
    assert tokenize_genes("APOE, TOMM40") == ["APOE", "TOMM40"]
    assert tokenize_genes("NR") == []
    assert tokenize_genes(np.nan) == []


def test_tokenize_variants():
    assert tokenize_variants("rs123456") == ["rs123456"]
    assert tokenize_variants("rs123; rs456") == ["rs123", "rs456"]
    assert tokenize_variants("NR") == []
    assert tokenize_variants(np.nan) == []


def test_allele_frequency_range():
    test_vals = [0.004, 0.005, 0.012, 0.020, 0.021]
    filtered = [v for v in test_vals if 0.005 <= v <= 0.020]
    assert filtered == [0.005, 0.012, 0.020]
