from typing import Optional, Tuple, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy.external as sce


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
        fgsea supports: 'es', 'nes', 'pval', 'padj', 'log2error'
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
        dd = df.sort_values(by=score)
        fig = plt.figure(figsize=figsize)
        dd = df.sort_values(by=score).iloc[-top_term:,:]

        if score_cutoff:
            dd = dd[dd[score] > score_cutoff]

        ax = fig.add_subplot(111)

        # user defined colors for positive and negative values
        colors = color
        if (dd[score] > 0).any() and (dd[score] < 0).any():
            if len(color) == 2:
                color_tmp = [''] * len(dd[score])
            for i, s in enumerate(dd[score]):
                if s >= 0:
                    color_tmp[i] = color[0]
                else:
                    color_tmp[i] = color[1]
            colors = color_tmp

        bar = dd.plot.barh(
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
    input_gene_ranking_file: str,
    hallmark_gene_sets_file: Optional['str'] = None,
    hallmark_gene_sets_list: Optional[list] = None,
    type: str = 'gseapy',
    out_dir: Optional['str'] = None,
    rscript_path: Optional['str'] = None,
    verbosity: Optional[bool] = False,

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
        Name of .rnk or .csv file of two columns. First column is gene names and
        second column is their respective rankings by weight
    hallmark_gene_sets_file
        Name of .gmt file
    hallmark_gene_sets_list
        List of hallmark gene set names
    type
        Type of GSEA: “gseapy” (default) or “fgsea”
    out_dir
        Optional name of output file, if a string is provided it will be saved as .csv
    rscript_path
        Required for fgsea
    verbosity
        Passed in to execute_r_script. Default false

    Returns
    -------
    Pandas dataframe of GSEA results

    Notes
    -------
    For hallmark_gene_sets, user must either provide a .gmt file or list of genes

    Resulting table has columns of both gseapy and fgsea results, but may be NaN
    based on the type of GSEA chosen

    gseapy specific columns: 'fdr' (false discovery rate)
    fgsea specific columns: 'padj' (adjusted p-value), 'log2error'

    For both fgsea and gseapy, the bounds for the minimum and maximum number
    of genes that are also in the gene set are as follows
    minsize = 15
    maxsize = 500

    Example
    -------
    >>> gsea_df = scanpy.external.tl.gsea('rank.rnk', 'genes.gmt', 'gseapy')

    """
    if hallmark_gene_sets_file is not None and hallmark_gene_sets_list is not None:
        raise ValueError(
            'You can\'t specify both a hallmark_gene_sets_file, hallmark_gene_sets_list. ' 'Please select only one.'
        )

    if type == 'gseapy':

        import gseapy as gp

        # read in preranked gene list, .rnk or .csv accepted
        rnk = pd.read_csv(input_gene_ranking_file, header=None, sep="\t")

        # both .gmt or list of pathways accepted for hallmark_gene_sets
        hallmark_gene_sets = hallmark_gene_sets_list if hallmark_gene_sets_list else hallmark_gene_sets_file

        # run gseapy, returns a prerank object
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=hallmark_gene_sets,
            no_plot=True,
            seed=0
        )

        gseapy_df = pre_res.res2d

        if out_dir is not None:
            gseapy_df.to_csv(out_dir, index=True)

        return gseapy_df

    elif type == 'fgsea':
        # run fgsea
        print("fgsea")

        # create the arguments
        args = [input_gene_ranking_file]

        if hallmark_gene_sets_file is not None:
            args += ['--file', hallmark_gene_sets_file]
        else:
            # modify the list to be read in as r input
            # r_list = ",".join(hallmark_gene_sets_list)
            args += ['--list', hallmark_gene_sets_list]

        # execute R script
        # handles missing rscript_path error
        sce.tl.execute_r_script(rscript_path, './_scripts/fgsea.R', args, verbosity=verbosity)

        # TODO
        # read in intermediate file, process, and delete intermediate file

    else:
        raise ValueError(
            'Please select either "fgsea" or "gseapy" (default).'
        )
