#! /usr/bin/env python
"""
SNAKEMAKE - RNA seq merge counts
Version: 1.0
"""

import os
from sys import argv
from sys import stderr
from scipy.stats import rankdata
from bioinfokit.analys import norm
#from pygmnormalize.utils import percentile

import numpy as np
import pandas as pd

def total_count_normalization(matrix):
    """
    Total count normalization
    
    Parameters
    ----------
    matrix : array_like
        Matrix to normalize.
        
    Returns
    -------
    array_like
        Normalized matrix.
    """
    return matrix / matrix.sum(axis=0)

def percentile_normalization(matrix, p, saving_memory=False):
    """
    Percentile normalization
    
    Parameters
    ----------
    matrix : array_like
        Matrix to normalize.
    p : float in range of [0,100]
        Percentile to compute, which must be between 0 and 100 inclusive.
    saving_memory : bool
        Parameter for activation of RAM saving mode. This may take longer.
        
    Returns
    -------
    array_like
        Normalized matrix.
    """
    return matrix / percentile(matrix, p)

def quartile_normalization(matrix, q, saving_memory=False):
    """
    Quartile normalization

    Parameters
    ----------
    matrix : array_like
        Matrix to normalize.
    q : string from {"lower", "median", "upper"} or quartile number (1, 2 or 3)
        The names of quartiles to compute in accordance:
        "lower" = 1,
        "median" = 2,
        "upper" = 3.
    saving_memory : bool
        Parameter for activation of RAM saving mode. This may take longer.

    Returns
    -------
    array_like
        Normalized matrix.
    """
    d = {"upper": 75, "lower": 25, "median": 50, 3: 75, 1: 25, 2: 50}
    assert q in d, 'Unexpected quartile for normalization: "' + str(q) + '"'
    return percentile_normalization(matrix, d[q])

def percentile(matrix, p, saving_memory=False):
    """
    Estimation of percentile for each column without zero-rows.
    
    Parameters
    ----------
    matrix : array_like
        Matrix to calculate percentile.
    p : float in range of [0,100]
        Percentile to compute, must be between 0 and 100 inclusive.
    saving_memory : bool
        Parameter for activation of RAM saving mode. This may take longer.
        
    Returns
    -------
    array_like
        Calculated percentile for each column.
    """
    if saving_memory:
        if not isinstance(matrix, np.ndarray):
            matrix = np.array(matrix)
        mask = [np.any(r > 0) for r in matrix] 
        return np.array([np.percentile(c[mask], p) for c in matrix.T])
    
    return np.percentile(matrix[np.any(matrix > 0, axis=1)], p, axis=0)

def normalize(norm_type, column):
    f = open("merged.txt", "r")
    counts = pd.read_csv(f , sep="\t", index_col=list(range(column)))
    f.close()
    f = open("merged.txt", "r")
    gene_length = pd.read_csv(f , sep="\t")
    f.close()
    f = open("normalized.txt", "w")
    
    if norm_type == "Median_of_ratios":
        norm_counts = np.log2(normalize_counts(counts) + 1)
        norm_type = 'rnaseq_median_of_ratios'
    elif norm_type == "TMM":
        norm_counts = tmm_normalization(counts)
        norm_type = 'rnaseq_tmm'
    elif norm_type == "CPM":
        cpm_norm = (counts * 1e6) / counts.sum()
        norm_counts = cpm_norm
        norm_type = 'rnaseq_cpm'
    elif norm_type == "TPM":
        a = counts.divide(gene_length['Length'].iloc[0], axis='index') * 1e3
        norm_counts = (a * 1e6) / a.sum()
        norm_type = 'rnaseq_tpm'
    elif norm_type == "RPKM":
        norm_counts = (counts.divide(gene_length['Length'].index, axis='index') * 1e9) / counts.sum()
        norm_type = 'rnaseq_fpkm'
    elif norm_type == "Total_count":
        norm_counts = total_count_normalization(counts)
        norm_type = 'rnaseq_total'
    else :
        norm_counts = quartile_normalization(counts, "upper")

    temp_list_counts = []
    temp_list_norm = []
    samples = list(norm_counts.columns)
    
    for i in range(len(samples)):
        temp_series_norm = norm_counts.iloc[:,i]
        temp_series_counts =counts.iloc[:,i]

        series_norm = temp_series_norm.to_frame()
        series_counts = temp_series_counts.to_frame()

        series_norm.insert(0,'sample_id',samples[i])
        series_counts.insert(0,'sample_id',samples[i])

        series_norm.dropna('columns')
        series_counts.dropna('columns')

        temp_list_norm.append(series_norm)
        temp_list_counts.append(series_counts)
    
    norm_counts = pd.concat(temp_list_norm, ignore_index=False)
    norm_counts[norm_type] = norm_counts[norm_counts.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    counts = pd.concat(temp_list_counts, ignore_index=False)
    norm_counts['rnaseq_count'] = counts[counts.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    order = norm_counts.columns.tolist()
    norm_counts[['sample_id', 'rnaseq_count', norm_type]].to_csv(f , sep="\t", index=True)
    #norm_counts.to_csv(f , sep="\t", index=True)
    f.close()


def normalize_all():
    f = open("merged.txt", "r")
    counts = pd.read_csv(f , sep="\t", index_col=list(range(6)))
    f.close()
    f = open("merged.txt", "r")
    gene_length = pd.read_csv(f , sep="\t")
    f.close()
    f = open("normalized.txt", "w")

    #if norm_type == "Median_of_ratios":
    mor_counts = np.log2(normalize_counts(counts) + 1)
    #elif norm_type == "TMM":
    tmm_counts = tmm_normalization(counts)
    #elif norm_type == "CPM":
    cpm_counts = (counts * 1e6) / counts.sum()
    #cpm_counts = cpm_norm
    #elif norm_type == "TPM":
    a = counts.divide(gene_length['Length'].iloc[0], axis='index') * 1e3
    tpm_counts = (a * 1e6) / a.sum()
    #elif norm_type == "RPKM":
    rpkm_counts = (counts.divide(gene_length['Length'].index, axis='index') * 1e9) / counts.sum()
    
    temp_mor_list = []
    temp_tmm_list = []
    temp_cpm_list = []
    temp_tpm_list = []
    temp_rpkm_list = []

    samples = list(mor_counts.columns)
    for i in range(len(samples)):
        temp_mor_series = mor_counts.iloc[:,i]
        temp_tmm_series = tmm_counts.iloc[:,i]
        temp_cpm_series = cpm_counts.iloc[:,i]
        temp_tpm_series = tpm_counts.iloc[:,i]
        temp_rpkm_series = rpkm_counts.iloc[:,i]

        mor_series = temp_mor_series.to_frame()
        tmm_series = temp_tmm_series.to_frame()
        cpm_series = temp_cpm_series.to_frame()
        tpm_series = temp_tpm_series.to_frame()
        rpkm_series = temp_rpkm_series.to_frame()

        mor_series.insert(0,'sample_id',samples[i])
        tmm_series.insert(0,'sample_id',samples[i])
        cpm_series.insert(0,'sample_id',samples[i])
        tpm_series.insert(0,'sample_id',samples[i])
        rpkm_series.insert(0,'sample_id',samples[i])

        mor_series.dropna('columns')
        tmm_series.dropna('columns')
        cpm_series.dropna('columns')
        tpm_series.dropna('columns')
        rpkm_series.dropna('columns')
        
        temp_mor_list.append(mor_series)
        temp_tmm_list.append(tmm_series)
        temp_cpm_list.append(cpm_series)
        temp_tpm_list.append(tpm_series)
        temp_rpkm_list.append(rpkm_series)

    mor_counts = pd.concat(temp_mor_list, ignore_index=False)
    mor_counts['median_of_ratios'] = mor_counts[mor_counts.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    tmm_counts = pd.concat(temp_tmm_list, ignore_index=False)
    mor_counts['tmm'] = tmm_counts[tmm_counts.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    cpm_counts = pd.concat(temp_cpm_list, ignore_index=False)
    mor_counts['cpm'] = cpm_counts[cpm_counts.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    tpm_counts = pd.concat(temp_tpm_list, ignore_index=False)
    mor_counts['tpm'] = tpm_counts[tpm_counts.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    rpkm_counts = pd.concat(temp_rpkm_list, ignore_index=False)
    mor_counts['rpkm'] = rpkm_counts[rpkm_counts.columns[1:]].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1)
    
    mor_counts[['sample_id', 'median_of_ratios', 'tmm', 'cpm', 'tpm', 'rpkm']].to_csv(f , sep="\t", index=True)
    #norm_counts.to_csv(f , sep="\t", index=True)
    f.close()


def tmm_normalization(matrix, index_ref=None, trim_fold_change=0.3, trim_abs_expr=0.05, saving_memory=False):
    """
    Trimmed mean of M-values normalization

    Parameters
    ----------
    matrix : array_like
        Matrix to normalize.
    index_ref:
        Index of reference column.
    trim_fold_change:
        Percent of trimmed for folder change.
    trim_abs_expr:
        Percent of trimmed for absolute expression.
    saving_memory : bool
        Parameter for activation of RAM saving mode. This may take longer.

    Returns
    -------
    array_like
        Normalized matrix.
    """
    matrix_np = np.array(matrix)                      # better speed of calculating
    np.seterr(divide='ignore', invalid='ignore')      # for divide on zeros in log2

    # Calculation log2(tmm_factor)
    def log2_tmm(index_vec):
        # select the necessary vectors
        curr_vec = matrix_np[:, index_vec]
        ref_vec = matrix_np[:, index_ref]

        # total number molecules in cells
        total_curr_vec = np.sum(curr_vec)
        total_ref_vec = np.sum(ref_vec)

        # select significant genes
        check_inf = (~np.isinf(matr_a[:, index_vec])) & (~np.isinf(matr_m[:, index_vec]))
        ranks = rankdata(matr_a[:, index_vec][check_inf], method='ordinal')
        bool_a = (ranks > len(ranks) * trim_abs_expr) & (ranks < len(ranks) * (1 - trim_abs_expr))
        ranks = rankdata(matr_m[:, index_vec][check_inf], method='ordinal')
        bool_m = (ranks > len(ranks) * trim_fold_change) & (ranks < len(ranks) * (1 - trim_fold_change))
        curr_vec = curr_vec[check_inf]
        ref_vec = ref_vec[check_inf]
        bool_curr_vec = curr_vec > 0
        bool_ref = ref_vec > 0
        bool_result = bool_curr_vec & bool_ref & bool_a & bool_m

        # calculation of required values
        w_vec = 1 / ((total_curr_vec - curr_vec[bool_result]) / (total_curr_vec * curr_vec[bool_result]) +
                     (total_ref_vec - ref_vec[bool_result]) / (total_ref_vec * ref_vec[bool_result]))
        m_vec = np.log2(curr_vec[bool_result] / total_curr_vec) - np.log2(ref_vec[bool_result] / total_ref_vec)

        # calculation log2(tmm_factor)
        w_sum = np.sum(w_vec)
        if np.isclose(w_sum, 0) or np.isinf(w_sum):
            print("Unexpected sum of weights for vector {}: '{}'".format(index_vec, w_sum), file=stderr)
            return 0

        return np.sum(w_vec * m_vec) / w_sum

    # find index of reference column
    f75 = percentile(matrix_np, 75, saving_memory)
    if index_ref is None:
        index_ref = np.argmin(abs(f75 - np.mean(f75)))
    elif not isinstance(index_ref, int) and isinstance(matrix, pd.DataFrame):
        index_ref = np.where(matrix.columns.values == (index_ref))[0][0]

    # find matrix A and M described expression levels of genes
    matr_norm = matrix_np / np.sum(matrix_np, axis=0)
    matr_a = np.log2(matr_norm * matr_norm[:, index_ref].reshape(matr_norm.shape[0], 1)) / 2
    matr_m = np.log2(matr_norm / matr_norm[:, index_ref].reshape(matr_norm.shape[0], 1))

    # calculation tmm_factor and normalization of input data
    tmm_factor = 2 ** np.array([log2_tmm(i) for i in range(matrix_np.shape[1])])
    return matrix / tmm_factor


def normalize_counts(counts):
    """Normalizes expression counts using DESeq's median-of-ratios approach."""

    with np.errstate(divide="ignore"):
        size_factors = estimate_size_factors(counts)
        return counts / size_factors


def estimate_size_factors(counts):
    """Calculate size factors for DESeq's median-of-ratios normalization."""

    def _estimate_size_factors_col(counts, log_geo_means):
        log_counts = np.log(counts)
        mask = np.isfinite(log_geo_means) & (counts > 0)
        return np.exp(np.median((log_counts - log_geo_means)[mask]))

    log_geo_means = np.mean(np.log(counts), axis=1)
    size_factors = np.apply_along_axis(
        _estimate_size_factors_col, axis=0,
        arr=counts, log_geo_means=log_geo_means)

    return size_factors


def merge(input_file, columns):
    frames = (pd.read_csv(fp, sep="\t", index_col=list(range(columns-1)))
        for fp in input_file)
    merged = pd.concat(frames, axis=1)
    
    # Extract sample names.
    f = open("merged.txt", "w")
    merged.to_csv(f , sep="\t", index=True)
    f.close()


def main():
    norm_type=argv[1]
    files = argv[2:]
    df = pd.read_csv(files[0], sep="\t")
    cols=len(df.axes[1])
    del df
    merge(files, cols)
    if norm_type == 'ALL':
        normalize_all()
    else:
        normalize(norm_type, cols)

if __name__ == '__main__':
    main()
