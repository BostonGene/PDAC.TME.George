import warnings

import numpy as np
import pandas as pd


class GeneSet(object):
    def __init__(self, name, descr, genes):
        self.name = name
        self.descr = descr
        self.genes = set(genes)
        self.genes_ordered = list(genes)

    def __str__(self):
        return '{}\t{}\t{}'.format(self.name, self.descr, '\t'.join(self.genes))
    
def read_gene_sets(gmt_file):
    """
    Return dict - {geneset_name : GeneSet object}

    :param gmt_file: str, path to .gmt file
    :return: dict
    """
    gene_sets = {}
    with open(gmt_file) as handle:
        for line in handle:
            items = line.strip().split('\t')
            name = items[0].strip()
            description = items[1].strip()
            genes = set([gene.strip() for gene in items[2:]])
            gene_sets[name] = GeneSet(name, description, genes)

    return gene_sets


def ssgsea_score(ranks, genes):
    """
    In-house implementation of ssGSEA score
    
    :param ranks: pd.DataFrame, DataFrame with calculated ranks of genes with genes in index and samples in columns 
    :param genes: list, list of genes (a.k.a geneset) which ssgsea scores should be calculated for
    :return: pd.Series, calculated ssGSEA scores for used geneset
    """

    common_genes = list(set(genes).intersection(set(ranks.index)))
    if not len(common_genes):
        return pd.Series([0.0]*len(ranks.columns), index=ranks.columns)
    sranks = ranks.loc[common_genes]
    return (sranks**1.25).sum()/(sranks**0.25).sum() - (len(ranks.index) - len(common_genes) + 1) / 2


def ssgsea_formula(data, gene_sets, rank_method='max'):
    """
    Return DataFrame with ssgsea scores
    Only overlapping genes will be analyzed

    :param data: pd.DataFrame, DataFrame with samples in columns and variables in rows
    :param gene_sets: dict, keys - processes, values - bioreactor.gsea.GeneSet
    :param rank_method: str, 'min' or 'max'. Ranking method for pandas
    :return: pd.DataFrame, ssgsea scores, index - genesets, columns - patients
    """
    ranks = data.rank(method=rank_method, na_option='bottom')
    return pd.DataFrame(
        {
            gs_name: ssgsea_score(ranks, gene_sets[gs_name].genes)
            for gs_name in list(gene_sets.keys())
        }
    ).T


def median_scale(data, clip=None):
    c_data = (data - data.median()) / data.mad()
    if clip is not None:
        return c_data.clip(-clip, clip)
    return c_data


def read_dataset(file, sep='\t', header=0, index_col=0, comment=None):
    return pd.read_csv(file, sep=sep, header=header, index_col=index_col,
                       na_values=['Na', 'NA', 'NAN'], comment=comment)


def item_series(item, indexed=None):
    """
    Creates a series filled with item with indexes from indexed (if Series-like) or numerical indexes (size=indexed)
    :param item: value for filling
    :param indexed:
    :return: pd.Series
    """
    if indexed is not None:
        if hasattr(indexed, 'index'):
            return pd.Series([item]*len(indexed), index=indexed.index)
        elif type(indexed) is int and indexed > 0:
            return pd.Series([item]*indexed, index=np.arange(indexed))
    return pd.Series()

def to_common_samples(df_list=()):
    """
    Accepts a list of dataframes. Returns all dataframes with only intersecting indexes
    :param df_list: list of pd.DataFrame
    :return: pd.DataFrame
    """
    cs = set(df_list[0].index)
    for i in range(1, len(df_list)):
        cs = cs.intersection(df_list[i].index)

    if len(cs) < 1:
        warnings.warn('No common samples!')
    return [df_list[i].loc[list(cs)] for i in range(len(df_list))]


def sort_by_terms_order(data, t_order: list, vector = None):
    """
    Sort "data" into blocks with values from "t_order". If "vector" is provided, sort each block by corresponding values in "vector"
    :param data: pd.Series
    :param t_order: list, values for blocks to sort "data" into
    :param vector: pd.Series, same index as data, which values to sort each block by
    :return: np.array, 1 dimensional
    """

    x = []
    for term in t_order:
        indices = data[data == term].index

        if len(indices):
            if vector is not None:
                x.append(vector.reindex(indices).dropna().sort_values().index)
            else:
                x.append(indices)

    return np.concatenate(x)
