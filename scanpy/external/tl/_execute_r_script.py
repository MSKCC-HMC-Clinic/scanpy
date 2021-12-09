from typing import Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess, os

import os

# idea, just be able to write 'Rscript' as the rscript_path in subprocess call...
# os.environ['R_HOME'] = 'C:/Program Files/R/R-4.1.2/bin' # or '/Library/Frameworks/R.framework/Versions/3.6/Resources/'
# os.environ['R_USER'] = 'C:/Program Files/R/R-4.1.2/bin' # or '/Library/Frameworks/R.framework/Versions/3.6/Resources/'

HERE = Path(__file__).parent 

def execute_r_script (
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
        optional list of arguments to add to command
    verbosity
        Option to write subprocess output to terminal. Default false, which sets
        stdout and stderr to NULL
    cache @ TODO
        Option to save intermediate R files. Default false.
   
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

    if rscript_path is None:
        raise ValueError('Please provide the local path to your Rscript.exe executable')
    if r_filename is None:
        raise ValueError('Please provide the filename of your R script')

    try:
        path2rscript = Path(HERE, r_filename)

        command = [rscript_path, path2rscript, '--vanilla']

        # add optional arguments for the R script
        if arguments is not None:
            command += arguments
        
        print(command)

        if verbosity == False: 
            exit_code = subprocess.call(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        else:
            exit_code = subprocess.call(command)

        if exit_code == 0:
            return True
        else:
            return False

    except:
        return False
