from typing import Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import shutil

from ... import settings


HERE = Path(__file__).parent


def execute_r_script(
    rscript_path: str = None,
    r_filename: str = None,
    arguments: Optional[list] = None,
    verbosity: Optional[bool] = False,
):
    """\

    Helper function that runs given R file in terminal using python subprocess module with the Rscript command

    Parameters
    ----------
    rscript_path
        Path to local directory of Rscript.exe
    r_filename
        Filename of R script to run
    arguments
        optional list of arguments to add to command, such as input filenames specific
        to the r_filename.R script
    verbosity
        Option to write subprocess output to terminal. Default false, which sets
        stdout and stderr to NULL

    Returns
    -------
    Boolean of whether script ran successfully

    Notes
    -------
    This function requires the local path to the Rscript command to run the
    provided script using the subprocess module

    Any files created by running the r_filename.R script will be saved in a
    scanpy._settings.ScanpyConfig.cachedir directory (default './cache/')
    Any external function that calls execute_r_script will be
    responsible for removing this cache directory and its files

    The r_filename.R script provided as a parameter must be contained in the
    scanpy/external/tl/_scripts directory

    Example
    -------
    >>> is_success = scanpy.external.tl.execute_r_script('PATH_TO_Rscript.exe', 'test.R')

    """

    if rscript_path is None:
        raise ValueError('Please provide the local path to your Rscript.exe executable')
    if r_filename is None:
        raise ValueError('Please provide the filename of your R script')

    try:

        # create cache directory for temporary files
        cachedir = settings.cachedir
        if not os.path.exists(cachedir):
            os.mkdir(cachedir)

        rfilepath = os.path.join(HERE,'_scripts', r_filename)

        command = [rscript_path, rfilepath, '--vanilla', cachedir]

        # add optional arguments that are specific to the r_filename.R script
        if arguments is not None:
            command += arguments

        if not verbosity:
            exit_code = subprocess.call(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        else:
            exit_code = subprocess.call(command)

        if exit_code == 0:
            return True
        else:
            return False

    except os.error:
        return False


def clear_cache():

    """\

    Helper function that recursively removes the scanpy._settings.ScanpyConfig.cachedir directory and all files in it

    Notes
    -------
    This function would be called in a python script after processing any files generated
    from the execute_r_script call

    Example
    -------
    # assume running 'test.R' generates a 'test.csv' file in the settings.cachedir directory
    >>> is_success = scanpy.external.tl.execute_r_script('PATH_TO_Rscript.exe', 'test.R')
    # process 'test.csv'
    >>> test_csv = pd.read_csv('test.csv')
    # remove settings.cachedir directory and its contents
    >>> clear_cache()
    """
    temp_dir = settings.cachedir
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


def remove_cache_file(
    filename: str = None
):

    """\

    Helper function that removes a specified file from the scanpy._settings.ScanpyConfig.cachedir directory

    Notes
    -------
    This function would be called in a python script after processing any files generated
    from the execute_r_script call

    Example
    -------
    # assume running 'test.R' generates a 'test.csv' and 'test.txt' file in the settings.cachedir directory
    >>> is_success = scanpy.external.tl.execute_r_script('PATH_TO_Rscript.exe', 'test.R')
    # process 'test.csv'
    >>> test_csv = pd.read_csv('test.csv')
    # remove just the 'test.csv' file
    >>> remove_cache_file('test.csv')
    """
    temp_file = os.path.join(settings.cachedir, filename)
    if os.path.exists(temp_file):
        os.remove(temp_file)
