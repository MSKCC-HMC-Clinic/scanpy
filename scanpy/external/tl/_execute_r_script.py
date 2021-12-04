from typing import Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess, os

HERE = Path(__file__).parent 

def testR (rscript="test.R"):

    path2rscript = Path(HERE, rscript)
    command = "C:/Program Files/R/R-4.1.2/bin/Rscript"
    arg = '--vanilla'
    subprocess.call ([command, arg, path2rscript])



