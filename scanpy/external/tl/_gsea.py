from typing import Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy.external as sce
import os
import shutil
import pandas.api.types as ptypes
from ... import settings

HERE = Path(__file__).parent


def gsea_barplot(
    df: pd.DataFrame,
    score: Optional['str'] = 'nes',
    score_cutoff: Optional[float] = None,
    top_term: Optional[float] = 10,
    figsize: Optional[Tuple[float, float]] = (6.5,6),
    title: Optional['str'] = None,
    x_label: Optional['str'] = None,
    y_label: Optional['str'] = None,
    out_dir: Optional['str'] = 'gseabarplot.png',
    color: Optional[Union['str',list]] = '#1f77b4',
):

    """\
    Horizontal barplot of GSEA results of genes by enrichment terms given a dataframe
    Allows visualizations of either pval, nes, and fdr from outputs using either gseapy or fgsea.
    Dataframe generated from scanpy.external.tl.gsea() function.

    Parameters
    ----------
    df
        Dataframe where each row is an enrichment term with GSEA scores as columns
    score
        Chosen score to display. Default is 'nes'
        Must be exact name of column as it appears in dataframe.
        gseapy supports: 'es', 'nes', 'pval', 'fdr'
        fgsea supports: 'es', 'nes', 'pval', 'padj', 'log2err'
    score_cutoff
        Cutoff for displayed values.
        Example usage: when displaying p-values (using score='pval'), the user
        can define a cutoff for significance, such as pval < 0.05
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
    color
        Either string or list of colors when mapping positive and negative colors
        Default color is matplotlib default: https://matplotlib.org/stable/users/dflt_style_changes.html#id1

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
        fig = plt.figure(figsize=figsize)
        df = df.sort_values(by=score).iloc[-top_term:,:]

        if score_cutoff:
            df = df[df[score] > score_cutoff]

        ax = fig.add_subplot(111)

        # user defined colors for positive and negative values
        colors = color
        if (df[score] > 0).any() and (df[score] < 0).any():
            if len(color) == 2:
                color_tmp = [''] * len(df[score])
            for i, s in enumerate(df[score]):
                if s >= 0:
                    color_tmp[i] = color[0]
                else:
                    color_tmp[i] = color[1]
            colors = color_tmp

        bar = df.plot.barh(
            y=score,
            alpha=0.75,
            color=colors,
            fontsize=16,
            ax=ax)

        xlabel = x_label if x_label else score
        ylabel = y_label if y_label else "Term"
        bar.set_xlabel(xlabel, fontsize=16, fontweight='bold')
        bar.set_ylabel(ylabel, fontsize=16, fontweight='bold')

        bar.legend_.remove()
        bar.set_title(title, fontsize=16, fontweight='bold')

        plt.savefig(out_dir)

        return ax

    except KeyError:
        print('ERROR: key \'' + str(score) + '\' is not in provided dataframe')


def gsea(
    input_gene_ranking_file: str = None,
    hallmark_gene_sets_file: str = None,
    type: str = 'gseapy',
    out_dir: Optional['str'] = None,
    rscript_path: Optional['str'] = None,
    verbosity: Optional[bool] = False,
    cache: Optional[bool] = False
):

    """\

    Python interface for 2 preranked GSEA functions: gseapy (Fang, 2020) and fgsea (Korotkevich et al., 2019).
    Users provide a set of genes to be considered, a set of enrichment terms
    (hallmark gene sets that are a defined library or set of pathways),
    and select either gseapy (default) or fgsea.
    The set of genes may be optionally pre-ranked with user-defined weights.
    The output is a table of results for each enrichment term.

    Parameters
    ----------
    input_gene_ranking_file
        Name of .rnk or .csv file of two columns. First column is gene names (string type) and
        second column (numberic type) is their respective rankings by weight
    hallmark_gene_sets_file
        Name of .gmt or .csv file
    type
        Type of GSEA: “gseapy” (default) or “fgsea”
    out_dir
        Optional name of output file, if a string is provided it will be saved as .csv
    rscript_path
        Required for fgsea
    verbosity
        Passed in to execute_r_script. Default false
    cache
        for fgsea only: to keep any intermediate files generated by fgsea.R
        Default false

    Returns
    -------
    Pandas dataframe of GSEA results with the following columns:
        {es: enrichment score,
        nes: normalized enrichment score,
        pval: p-value,
        padj: BH adjusted p-value (fgsea specific),
        log2err: expected error p-value logarithm standard deviation (fgsea specific),
        fdr: false discovery rate, (gseapy specific)
        geneset_size: size of gene set, (gseapy specific)
        matched_size: genes matched to the data,
        genes: gene names from the data set, (gseapy specific)
        ledge_genes: leading edge genes
        }

    Notes
    -------
    For hallmark_gene_sets, user must either provide a .gmt file or .csv of genes

    Output csv has columns: ['es', 'nes', 'pval', 'padj', 'log2err', 'fdr', 'geneset_size', 'matched_size', 'genes', 'ledge_genes']]

    Resulting table has columns of both gseapy and fgsea results, but may be NaN
    based on the type of GSEA chosen

    gseapy specific columns: 'fdr' (false discovery rate), 'geneset_size', 'genes'
    fgsea specific columns: 'padj' (adjusted p-value), 'log2err'

    For both fgsea and gseapy, the bounds for the minimum and maximum number
    of genes that are also in the gene set are as follows
    minsize = 15
    maxsize = 500

    Example
    -------
    >>> gsea_df = scanpy.external.tl.gsea('rank.rnk', 'genes.gmt', 'gseapy')
    >>> fgsea_df = scanpy.external.tl.gsea(input_gene_ranking_file='./_data/gsea_data.gsea_data.rnk', hallmark_gene_sets_file='./_data/h.all.v7.4.symbols.gmt', type='fgsea', out_dir='_data/test_fgsea_df.csv', rscript_path=rscript_path)


    """

    if type not in ['gseapy', 'fgsea']:
        raise ValueError(
            'Please select either `fgsea` or `gseapy` (default).'
        )

    if Path(input_gene_ranking_file).suffix not in ['.rnk', '.csv']:
        raise ValueError(
                        'Only .rnk and .csv are supported for input gene ranking file'
                    )

    if Path(hallmark_gene_sets_file).suffix not in ['.gmt', '.csv']:
        raise ValueError(
                        'Only .gmt and .csv are supported for input hallmark gene sets file'
                    )

    # create cachedir if applicable
    cachedir = settings.cachedir
    if not os.path.exists(cachedir):
        os.mkdir(cachedir)

    if type == 'gseapy':

        try:
            import gseapy as gp
        except ImportError:
            raise ImportError('Please install GSEApy: `pip install gseapy`.')

        # gseapy has a logging bug, where it is never closed and makes it difficult
        # to remove files created by running gsea, currently putting it into the cache
        # but can't clear the cache during runtime, which does not match `cache` parameter above
        gsea_cache = os.path.join(cachedir, 'gseapy')
        if not os.path.exists(gsea_cache):
            os.mkdir(gsea_cache)

        # read in preranked gene list, .rnk or .csv accepted
        if Path(input_gene_ranking_file).suffix == '.rnk':
            # .rnk files are tab delimited
            rnk = pd.read_csv(input_gene_ranking_file, header=None, sep="\t")
        else:
            # assert that the .csv has two columns, one of type str for gene names and another of type float for rankings
            rnk = pd.read_csv(input_gene_ranking_file, header=None)

            try:
                assert ptypes.is_string_dtype(rnk.iloc[:,0])
                assert ptypes.is_numeric_dtype(rnk.iloc[:,1])
            except AssertionError:
                raise AssertionError('.csv format for input_gene_ranking_file is string type for gene list (first column) and numeric type for gene ranking values (second column).')

        # run gseapy, returns a prerank object
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=hallmark_gene_sets_file,
            outdir=str(gsea_cache),
            no_plot=True,
            seed=0
        )

        gseapy_df = pre_res.res2d

        # add fgsea columns to gseapy so outputs of both functions are equal
        gseapy_df['padj'] = np.NaN
        gseapy_df['log2err'] = np.NaN

        # change order of columns
        gseapy_df = gseapy_df[['es', 'nes', 'pval', 'padj', 'log2err', 'fdr', 'geneset_size', 'matched_size', 'genes', 'ledge_genes']]

        if out_dir is not None:
            gseapy_df.to_csv(out_dir, index=True)

        # gseapy_temp_dir = 'GSEA_Prerank'
        # if os.path.exists(gseapy_temp_dir):
        #     shutil.rmtree(gseapy_temp_dir)

        return gseapy_df

    else:

        args = [input_gene_ranking_file, hallmark_gene_sets_file]

        sce.tl.execute_r_script(rscript_path, 'fgsea.R', args, verbosity=verbosity)

        # read in the file generated by fgsea.R
        fgseaRes_path = os.path.join(cachedir, "fgseaRes.csv")
        fgsea_df = pd.read_csv(fgseaRes_path)

        if not cache:
            sce.tl.clear_cache()

        # rename columns for shared output format
        fgsea_df = fgsea_df.rename(columns={'pathway':'Term', 'ES':'es', 'NES':'nes', 'leadingEdge': 'ledge_genes','size':'matched_size'})
        fgsea_df = fgsea_df.set_index('Term')

        # gsea uses ' ' while gseapy uses ';' as delimiter for ledge_genes, match gseapy format
        fgsea_df = fgsea_df.replace(to_replace=r' ', value=';', regex=True)

        # add gseapy columns to fgsea so outputs of both functions are equal
        fgsea_df['fdr'] = np.NaN
        fgsea_df['geneset_size'] = np.NaN
        fgsea_df['genes'] = np.NaN

        fgsea_df = fgsea_df[['es', 'nes', 'pval', 'padj', 'log2err', 'fdr', 'geneset_size', 'matched_size', 'genes', 'ledge_genes']]

        if out_dir is not None:
            fgsea_df.to_csv(out_dir, index=True)

        return fgsea_df

