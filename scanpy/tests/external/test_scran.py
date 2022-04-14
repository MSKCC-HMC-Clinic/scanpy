import pytest
from typing import Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scanpy.external as sce
import os
import shutil
import pandas.api.types as ptypes
from anndata import AnnData

# pytest.importorskip("scran")

def test_scran(rscript_path):
    pbmc = sc.datasets.pbmc3k()
    sce.tl.scran(pbmc, rscript_path, verbosity=True)

