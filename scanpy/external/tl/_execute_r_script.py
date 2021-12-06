from typing import Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess, os

HERE = Path(__file__).parent 

def testR (
    r_filename: str = 'test.R'
):
    """\
    
    Parameters
    ----------
    r_filename
        Filename of R script to run
   
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
    command = "/Library/Frameworks/R.framework/Versions/3.6/Resources/Rscript"
    arg = '--vanilla'
    subprocess.call ([command, arg, path2rscript])