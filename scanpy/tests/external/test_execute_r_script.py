import pytest

from pathlib import Path
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import os
import filecmp
import shutil
from scanpy import settings

# pytest.importorskip("execute_r_script")


def test_arguments(
    rscript_path=None,
    verbosity=False
):
    """Test providing arguments to python subprocess command and cache creation and removal"""

    # clear cache
    cachedir = settings.cachedir
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)

    arguments = ['test.txt', 'contents']

    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='example_arguments.R', arguments=arguments, verbosity=verbosity)

    # define paths for test and master files for comparison
    test_file_path = os.path.join(cachedir,'test.txt')
    master_file_path = './_data/test.txt'

    # test that the file created by arguments provided to execute_r_script
    assert os.path.exists(test_file_path)
    assert filecmp.cmp(test_file_path, master_file_path, shallow=False)

    # test the remove_cache_file and clear_cache functions
    sce.tl.remove_cache_file('test.txt')
    assert not os.path.exists(test_file_path)
    sce.tl.clear_cache()
    assert not os.path.exists(cachedir)


def test_fgsea(
    rscript_path=None,
    verbosity=False
):
    """Test calling fgsea on default arguments and reproducabililty.
        Previously tested against a set output (from running it once), but this output
        differs between MAC OS and Windows machines.
    """

    # clear cache
    cachedir = settings.cachedir
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)
    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='example_fgsea.R', arguments=['first_fgseaRes.csv'], verbosity=verbosity)

    sce.tl.execute_r_script(rscript_path, r_filename='example_fgsea.R', arguments=['second_fgseaRes.csv'], verbosity=verbosity)

    # define paths for first and second files for comparison
    first_test_file_path = os.path.join(cachedir,'first_fgseaRes.csv')
    second_test_file_path = os.path.join(cachedir,'second_fgseaRes.csv')

    # test that the file created by arguments provided to execute_r_script
    assert os.path.exists(first_test_file_path)
    assert os.path.exists(second_test_file_path)
    # assert filecmp.cmp(test_file_path, master_file_path, shallow=False)

    first_fgsea_sample_df = pd.read_csv(first_test_file_path)
    second_fgsea_sample_df = pd.read_csv(second_test_file_path)
    assert first_fgsea_sample_df.equals(second_fgsea_sample_df)

    # test the remove_cache_dir function to delete all files under directory
    sce.tl.clear_cache()
    assert not os.path.exists(cachedir)
