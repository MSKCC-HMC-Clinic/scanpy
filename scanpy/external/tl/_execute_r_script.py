from typing import Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess, os

HERE = Path(__file__).parent 

def execute_r_script (
    rscript_path: str,
    arguments: Optional['list'] = None,
    r_filename: str = 'test.R',
):
    """\
    
    Parameters
    ----------
    rscript_path
        Path to local directory of Rscript.exe
    r_filename
        Filename of R script to run
    arguments
        optional list of arguments to add to command
   
    Returns
    -------
    Boolean of whether script ran successfully
    
    Notes
    -------
	This function requires the local path to the Rscript command to run the 
    provided script using the subprocess module

    This function will prompt the user in terminal for the Rscript command path

    Example
    -------
    >>> is_success = scanpy.external.tl.execute_r_script('test.R')
 
    """
 
    path2rscript = Path(HERE, r_filename)
    # command = "/Library/Frameworks/R.framework/Versions/3.6/Resources/Rscript"
    # command = rscript_path
    arg = '--vanilla'

    # add optional arguments for the R script
    if arguments is not None:
        arguments += [arg]
        arg = " ".join(arguments)
    subprocess.call ([rscript_path, arg, path2rscript])