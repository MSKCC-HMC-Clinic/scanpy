import pytest

from pathlib import Path
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import os

# pytest.importorskip("execute_r_script")


def test_arguments(
    rscript_path=None
):
    arguments = ['_data/test.csv', 'a,b,c']

    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='test_arguments.R', arguments=arguments)

def test_fgsea(
    rscript_path=None
):

    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='example_fgsea.R')
    # run time was 30.9 sec
    # sce.tl.remove_temp()
