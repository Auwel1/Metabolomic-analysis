# -*- coding: utf-8 -*-

"""
Author : COFFE Aurelien
bunch of functions to work on metabolomic data
"""

import re
import pandas as pd
import numpy as np
import scipy.stats as sc
import _utils as ut


# separate amino acid from other compounds

def data_reader(datapath):
    xl = re.compile(r".*xl.*$")

    if xl.match(datapath):
        data_frame = pd.read_excel(datapath, header=0)
    else:
        data_frame = pd.read_csv(datapath, header=0, sep=None, index_col=0, engine='python')
    return data_frame


def amino_extract(data):
    aa_table = data.loc[data['metabolite'].isin(ut.aa)]
    noaa_table = data.loc[not data['metabolite'].isin(ut.aa)]
    return aa_table, noaa_table


def mean_extract(df):
    """
    extract mean columns into a new dataframe
    """
    pat = re.compile(r".*_mean")
    subset = df[list(df.columns[df.columns.str.match(pat)])]
    return subset


def sd_extract(df):
    """
    extract mean columns into a new dataframe
    """
    pat = re.compile(r".*_sd")
    subset = df[list(df.columns[df.columns.str.match(pat)])]
    return subset


def lineage_extract(df, r):
    """
      create a list with the cellular lineage names
      """
    history = []
    for i in df.columns:
        m = re.match(r, i)
        if m is not None:
            lineage = m.group(2)
            if lineage not in history:
                history.append(lineage)
    return history


# extract a subset from a sample
def sample_extract(df, sample):
    p = re.compile(sample)
    df_sample = df[df['sample'].str.match(p)]
    return df_sample


# extract a subset from precise metabolite
def metabolite_extract(df, metabo):
    metabo_sub = df[df['metabolite'] == metabo]
    return metabo_sub


# sum for the isotopologue
def iso_sum(df, col):
    sums = []
    for i in set(df['sample']):
        tmp = df[df['sample'] == i]
        summ = tmp[col].aggregate(func=sum)
        sums.append(summ)

    return sums


def cell_list(cell_file):
    return pd.read_csv(cell_file, header=None, index_col=0)


# create the isotopologue mean fractions dataframe
def iso_lists(df, mean_list, sd_list):

    for z in set(df["isotopologue"]):  # for each isotopologue
        iso = df[df["isotopologue"] == z]  # create a new subset again with the same isotopologue from each replicate
        # measure the mean and standard deviation between replicate and add them to a list
        mean_list.append(np.mean(list(iso["isotopologue_fraction"])))
        sd_list.append(np.std(list(iso["isotopologue_fraction"]), ddof=1))
    return mean_list, sd_list

def iso_means(mean_list, sd_list, result, samp):
    # add the lists as new columns to the dataframe
    if len(mean_list) != 0:
        result[samp + '_mean'] = mean_list
        result[samp + '_sd'] = sd_list
    return result


# create the mean C13 enrichment dataframe
def enrichment_lists(df, mean_list, sd_list):
    enr = list(set(df['mean_enrichment']))  # extract the mean enrichment

    # creating the columns
    mean_list.append(np.mean(enr) * 100)  # measure the mean between replicate (transformed in a percentage)
    if len(enr) == 1:
        sd_list.append(0)
    else:
        sd_list.append(np.std(enr, ddof=1) * 100)  # measure the standard deviation between replicate (percentage)
    return mean_list, sd_list


def enrichment_means(mean_list, sd_list, enr_means, samp):
    # add the columns to the dataframe
    if len(mean_list) != 0:
        enr_means[samp + '_mean'] = mean_list
        enr_means[samp + '_sd'] = sd_list
    return enr_means


# create the dataframe with corrected global area from each metabolite and each cell type

def area_lists(df, mean_list, sd_list):
    sums = iso_sum(df, 'corrected_area')  # measure the sum of areas of a metabolite

    # create the means and standard deviation columns (measured between each replicate)
    mean_list.append(sc.gmean(sums))
    sd_list.append(np.std(sums, ddof=1))
    return mean_list, sd_list

def area_means(mean_list, sd_list, area_mean, samp, cellnumber=None):
    # add columns to the dataframe
    if len(mean_list) != 0:
        if cellnumber is not None:
            area_mean[samp + '_mean'] = list(map((lambda x: x / float(cellnumber.loc[samp])), mean_list))
            area_mean[samp + '_sd'] = list(map((lambda x: x / float(cellnumber.loc[samp])), sd_list))
        else:
            area_mean[samp + '_mean'] = mean_list
            area_mean[samp + '_sd'] = sd_list

    return area_mean

