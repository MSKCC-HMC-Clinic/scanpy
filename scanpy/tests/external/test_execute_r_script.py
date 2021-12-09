import pytest

import scanpy as sc
import pandas as pd
import scanpy.external as sce

pytest.importorskip("execute_r_script")

def test_arguments(
    rscript_path=None
):
    import scanpy.external as sce
    arguments = ['_data/test.csv', 'a,b,c']

    sce.tl.execute_r_script(rscript_path, r_filename='./_scripts/test_arguments.R', arguments=arguments)
    #TODO: create tests for different sizes of .csv's.../anndata objects