# import pytest

# import scanpy as sc
# import scanpy.external as sce
# import pandas as pd
# from matplotlib.testing.compare import compare_images
# from pathlib import Path

# pytest.importorskip("cluster_covariance")

# HERE: Path = Path(__file__).parent
# # ROOT = HERE / '_images'
# # FIGS = HERE / 'figures'

# def test_joyplot():
#     """
#     Test that joyplt works
#     """
#     adata = sc.datasets.pbmc68k_reduced()
#     #     adata.to_df()
#     df = sce.tl.cluster_covariance(adata)
#     # df
#     sce.tl.joyplot(df = df, grouping = "leiden", title = "Gene expression per cluster", 
#                             color = ['skyblue','lightpink','red'], x_range = (-4,4), output_file = "joyplottest.png",x_label="Expression range")


#     assert compare_images('joyplottest.png',str(HERE)+'/master_joyplot.png', tol = 10)

# def main():
#     test_joyplot()
#     print("Test passed")

# if __name__=="__main__":
#     main()   