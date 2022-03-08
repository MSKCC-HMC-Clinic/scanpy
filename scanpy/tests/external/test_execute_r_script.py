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

    cachedir = settings.cachedir
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)

    arguments = ['test.txt', 'contents']
    
    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='example_arguments.R', arguments=arguments, verbosity=verbosity)

    # define paths for test and master files for comparison
    test_file_path = os.path.join(cachedir,'test.txt')
    print(test_file_path)
    master_file_path = './_data/test.txt'

    # test that the file created by arguments provided to execute_r_script
    assert os.path.exists(test_file_path)
    assert filecmp.cmp(test_file_path, master_file_path, shallow=False)

    # test the remove_cache_file and clear_cache functions
    # sce.tl.remove_cache_file('test.txt')
    assert not os.path.exists(test_file_path)
    # sce.tl.clear_cache()
    assert not os.path.exists(cachedir)


def test_fgsea(
    rscript_path=None,
    verbosity=False
):
    cachedir = settings.cachedir
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)
    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='example_fgsea.R', verbosity=verbosity)

    # define paths for test and master files for comparison
    test_file_path = os.path.join(cachedir,'fgseaRes.csv')
    master_file_path = './_data/fgseaRes.csv'

    # test that the file created by arguments provided to execute_r_script
    assert os.path.exists(test_file_path)
    assert filecmp.cmp(test_file_path, master_file_path, shallow=False)

    # test the remove_cache_dir function to delete all files under directory
    sce.tl.clear_cache()
    assert not os.path.exists(cachedir)
