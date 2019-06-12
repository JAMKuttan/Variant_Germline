import pytest
import pandas as pd
from io import StringIO
import os

test_output_path = os.path.dirname(os.path.abspath(__file__)) + \
		'/../output/'

@pytest.mark.output
def test_output_VCF():
    assert os.path.exists(os.path.join(test_output_path, 'GM12878/GM12878.annot.vcf.gz'))
    annot_vcf = test_output_path + 'GM12878/GM12878.annot.vcf.gz'
    annot_output = open(annot_vcf, encoding="utf-8").readlines()
    assert 'FILTER=<ID=PASS,Description="All filters passed">' in annot_output[1]

@pytest.mark.output
def test_output_Coverage_Histogram():
    assert os.path.exists(os.path.join(test_output_path, 'GM12878.coverage_histogram.png'))
    assert os.path.getsize(os.path.join(test_output_path, 'GM12878.coverage_histogram.png')) > 0

@pytest.mark.output
def test_output_CDF():
    assert os.path.exists(os.path.join(test_output_path, 'coverage_cdf.png'))
    assert os.path.getsize(os.path.join(test_output_path, 'coverage_cdf.png')) > 0

@pytest.mark.output
def test_output_QC():
    assert os.path.exists(os.path.join(test_output_path, 'sequence.stats.txt'))
    assert os.path.getsize(os.path.join(test_output_path, 'sequence.stats.txt')) > 0
