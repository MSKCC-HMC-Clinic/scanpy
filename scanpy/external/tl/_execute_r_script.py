from typing import Optional, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess, os

HERE = Path(__file__).parent 

def testR (rscript="test.R"):
    path2rscript = Path(HERE, rscript)
    command = "/Library/Frameworks/R.framework/Versions/3.6/Resources/Rscript"
    arg = '--vanilla'
    subprocess.call ([command, arg, path2rscript])