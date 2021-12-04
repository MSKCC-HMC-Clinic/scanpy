from typing import Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess, os

HERE = Path(__file__).parent 

def testR (rscript="test.R"):

    # current_dir = Path(os.getcwd())
    # rdir = current_dir / rfile
    # scanpy uses here notation, but getting __file__ undefined errors
    # HERE = Path(__file__).parent 
    print(HERE)
    rdir = Path(HERE, rscript)

    # rcode = subprocess("./test.R")
    subprocess.call (["C:/Program Files/R/R-4.1.2/bin/Rscript", "--vanilla", rdir])
    # subprocess.call (["C:/Program Files/R/R-4.1.2/bin/Rscript", "--vanilla", "C:/Users/qjuli/Downloads/Harvey Mudd/Clinic/scanpy/scanpy/external/tl/test.R"])



