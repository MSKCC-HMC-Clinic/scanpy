from typing import Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def gsea_barplot(
    df: pd.DataFrame, 
    score: Optional['str'] = 'pval',
    score_cutoff: Optional[float] = None,
    top_term: Optional[float] = 10,
    figsize: Optional[Tuple[float, float]] = (6.5,6),
    title: Optional['str'] = None,
    x_label: Optional['str'] = None,
    y_label: Optional['str'] = None,
    out_dir: Optional['str'] = 'gseabarplot.png',
    color = '#1f77b4', # matplotlib default: https://matplotlib.org/stable/users/dflt_style_changes.html#id1
):
  
    """\
    Horizontal barplot of GSEA results of genes by enrichment terms given a dataframe

    Allows visualizations of either pval, nes, and fdr from outputs using either gseapy or fgsea.
    Dataframe likely generated from scanpy.external.tl.gsea() function.
    
    Parameters
    ----------
    df
        Dataframe where each row is an enrichment term with GSEA scores as columns
    score
        Chosen score to display.
        Must be exact name of column as it appears in dataframe.
        gseapy supports: 'es', 'nes', 'pval', 'fdr'
        fgsea supports: 'es', 'nes', 'pval', 'padj', 'log2error'
    score_cutoff
        Cutoff for displayed values.
        Example usage: for pval < 0.05
    top_term
        Cutoff for displayed number of enrichment terms (bars)
    figsize
        width, height
    title
        Legend title. Appears on top of the color bar. Use '\n' to add line breaks.
    x_label
        User defined x_label
    y_label
        User defined y_label
    out_file
        Name of output plot image (include ".png" extension).

    Returns
    -------
    Returns horizontal barplot with a bar for each of the top enrichment terms
    by the chosen score
    
    Notes
    -------
    Is designed to take in output df from scanpy.external.tl.gsea() function.
    The score supported is dependent on the type of gsea previously run, 
    either 'gseapy' or 'fgsea'

    Example
    -------
    >>> gsea_df = scanpy.external.tl.gsea(...)
    >>> gsea_barplot(df = gsea_df, score = 'pval')

    """
    try:
      dd = df.sort_values(by=score)
      fig = plt.figure(figsize=figsize)
      dd = df.sort_values(by=score).iloc[-top_term:,:]

      if score_cutoff:
        dd = dd[dd[score] > score_cutoff]  

      ax = fig.add_subplot(111)

      # user defined colors for positive and negative values
      if (dd[score] > 0).any() and (dd[score] < 0).any():
        if len(color) == 2:
          color_tmp = ['']*len(dd[score])
          for i, s in enumerate(dd[score]):
            if s >= 0:
              color_tmp[i] = color[0]
            else:
              color_tmp[i] = color[1]
          color = color_tmp

      bar = dd.plot.barh( 
          y=score,
          alpha=0.75, 
          color = color,
          fontsize=16, 
          ax=ax)
      
      xlabel = x_label if x_label else score
      ylabel = y_label if y_label else "Term"
      bar.set_xlabel(xlabel, fontsize=16, fontweight='bold')
      bar.set_ylabel(ylabel, fontsize=16, fontweight='bold')

      bar.legend_.remove()
      bar.set_title(title, fontsize=24, fontweight='bold')

      plt.savefig(out_dir)     

    except KeyError:
      print('ERROR: key \'' + str(score) + '\' is not in provided dataframe')