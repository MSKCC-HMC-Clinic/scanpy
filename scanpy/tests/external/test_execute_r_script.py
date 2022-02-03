import pytest

from pathlib import Path
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import os
import filecmp

# pytest.importorskip("execute_r_script")

HERE = os.getcwd()
TEMP_DIR_PATH = HERE.replace('tests\\external','external\\tl\\_tmp')


def test_arguments(
    rscript_path=None
):
    arguments = ['test.txt', 'contents']

    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='example_arguments.R', arguments=arguments)

    # define paths for test and master files for comparison
    test_file_path = Path(TEMP_DIR_PATH,'test.txt')
    master_file_path = './_data/test.txt'

    # test that the file created by arguments provided to execute_r_script
    assert os.path.exists(test_file_path)
    assert filecmp.cmp(test_file_path, master_file_path, shallow=False)

    # test the remove_temp_file and remove_temp_dir functions
    sce.tl.remove_temp_file('test.txt')
    assert not os.path.exists(test_file_path)
    sce.tl.remove_temp_dir()
    assert not os.path.exists(TEMP_DIR_PATH)


def test_fgsea(
    rscript_path=None
):
    # .R script is located in /external/tl/_scripts
    sce.tl.execute_r_script(rscript_path, r_filename='example_fgsea.R')

    # define paths for test and master files for comparison
    test_file_path = Path(TEMP_DIR_PATH,'fgseaRes.csv')
    master_file_path = './_data/fgseaRes.csv'

    # test that the file created by arguments provided to execute_r_script
    assert os.path.exists(test_file_path)
    assert filecmp.cmp(test_file_path, master_file_path, shallow=False)

    # test the remove_temp_dir function to delete all files under directory
    sce.tl.remove_temp_dir()
    assert not os.path.exists(TEMP_DIR_PATH)
