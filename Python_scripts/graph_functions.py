# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import re
import data_pipe_functions as f
import _utils as ut
import seaborn as sns


def histo_print(df, pdf):
    nrows = 4
    flag = 0
    fig = plt.figure()
    for m in df.index:
        flag += 1
        if flag == 1: fig = plt.figure(figsize=(8, 15))
        if flag <= 4:
            ax = fig.add_subplot(nrows, 1, flag)
            ax.bar(df.columns, list(df.loc[m]))
            plt.title(m)
            plt.ylim([0, 100])
            plt.gcf().subplots_adjust(hspace=0.5)
        else:
            flag = 0
            pdf.savefig(fig)
            plt.close()



def hist_enrichment_single(data_frame):
    """
    Create PDF docs by sample with histograms about C13 enrichment for each metabolite
    :param data_frame:
    """
    r = ut.sample_graph_pattern
    sub = f.mean_extract(data_frame)
    history = f.lineage_extract(sub, r)

    for l in history:
        pdf = matplotlib.backends.backend_pdf.PdfPages(l + "_enrichment_hist.pdf")
        if len(l) < 2:
            lin_sub = sub[list(sub.columns[sub.columns.str.match(r'.*_' + l + '_.*')])]
            histo_print(lin_sub, pdf)

        else:
            lin_sub = sub[list(sub.columns[sub.columns.str.match(r'.*' + l + '.*')])]
            histo_print(lin_sub, pdf)
        pdf.close()


def time_samp(df):
    """
    extract time point through sample names
    """
    samples = sorted(df.columns)
    pat = ut.sample_graph_pattern
    times = []
    for i in samples:
        mat = re.match(pat, i)
        if mat is not None and mat.group(1) not in times:
            times.append(mat.group(1))
    return times


def pattern_extracting(df_mean, df_std, patt, means_list=None, sd_list=None):
    """
    This function extract subsets from tables which correspond to particular time points or samples

    :param patt:
    :param df_mean: main mean table
    :param df_std: main standard deviation table
    :param means_list: list of mean subsets
    :param sd_list: list of standard deviation subset
    :return:
    """
    time = re.compile(r".*\d.*")
    if re.match(time, patt):
        t_sub = df_mean[list(df_mean.columns[df_mean.columns.str.match(r'.*' + patt + '.*')])]
        s_sub = df_std[list(df_std.columns[df_std.columns.str.match(r'.*' + patt + '.*')])]
        means_list.append(t_sub)
        sd_list.append(s_sub)
        return means_list, sd_list
    else:
        t_sub = df_mean[list(df_mean.columns[df_mean.columns.str.match(r'.*_' + patt + '_.*')])]
        s_sub = df_std[list(df_std.columns[df_std.columns.str.match(r'.*_' + patt + '_.*')])]
        return t_sub, s_sub


def multi_metabo_enr(metabolite, flag, fig, mean_list, std_list, pdf, colors, pt):
    """
    sub_function to build enrichment histograms, work on the metabolite loop
    :param metabolite:
    :param flag: count the number of figures on a page, 4 figure per page
    :param fig: one page
    :param mean_list: list of mean subsets
    :param std_list: list of standard deviation subsets
    :param pdf: pdf file to save figures
    :return:
    """
    if flag <= 4:
        ax = fig.add_subplot(4, 1, flag)
        for i in range(0, len(mean_list)):
            ax.bar(mean_list[i].columns, mean_list[i].loc[metabolite], width=0.5, color=colors)
            ax.errorbar(mean_list[i].columns, mean_list[i].loc[metabolite],
                        std_list[i].loc[metabolite], fmt='none', colors='black')

        plt.title(metabolite)
        plt.ylim([0, 100])
        plt.gcf().subplots_adjust(hspace=0.7, top=0.95, bottom=0.1)
        plt.xticks(rotation=45)
    else:
        flag = 0
        pdf.savefig(fig)
        plt.close()
        pt -= 1
    pt += 1
    return flag, pt


def multi_samp_enr_hist(df, output, name):
    """
    Create PDF file with histograms presenting each metabolite C13 enrichment for each metabolite and each sample
    through time
    :param name: file name
    :param output: output directory path
    :param df: dataframe path
    :return:
    """
    df = f.data_reader(df)
    sub = f.mean_extract(df)
    sd = f.sd_extract(df)
    times = time_samp(sub)
    tables = []
    sds = []
    colors = ['red', 'orange', 'yellow', 'green']
    sns.set()
    pdf = matplotlib.backends.backend_pdf.PdfPages(output + name + ut.graph_types[0])
    for t in times:
        tables, sds = pattern_extracting(sub, sd, t, tables, sds)
    nrows = 4
    flag = 0
    fig = plt.figure(figsize=(8, 16))
    pt = 0
    while pt < len(sub.index):
        m = sub.index[pt]
        flag += 1
        if flag == 1:
            fig = plt.figure(figsize=(8, 16))
        flag, pt = multi_metabo_enr(m, flag, fig, tables, sds, pdf, colors, pt)
    pdf.close()


def metabo_CID(metabolite, flag, figure, mean_subset, std_subset, pdf, pt, fig):
    """
    sub function creating the different figures included in the CID graphs pdf

    :param metabolite: current metabolite
    :param flag: number of graph per figure
    :param figure: one page
    :param mean_subset: the subset with the mean values
    :param std_subset: subset with standard deviation values
    :param pdf: output pdf document
    :param pt: iterator applied to the while loop
    :return:
    """
    if flag <= 4:
        figure.add_subplot(2, 2, flag)
        m_sub = mean_subset[mean_subset.index.str.match(metabolite)]
        ms_sub = std_subset[std_subset.index.str.match(metabolite)]
        n = len(m_sub.index)
        bottom = 3 * [0]
        mvar = []

        for i in range(0, n, 1):
            mvar.append('M' + str(i))
            plt.bar(m_sub.columns, m_sub.iloc[i], bottom=bottom, width=0.5, yerr=ms_sub.iloc[i])
            bottom = bottom + m_sub.iloc[i]
        plt.ylim([0, 1.2])
        plt.title(metabolite)
        plt.legend(mvar, loc='upper center', bbox_to_anchor=(0.5, 1.4), ncol=5)
        plt.ylabel("C13 fractions")

    else:
        flag = 0
        pdf.savefig(figure)
        plt.close()
        pt -= 1
    pt += 1
    return flag, pt, fig


def CID_hist(df, output, name):
    """
    Create a pdf file including isotopologue fraction graphs

    :param name: file name
    :param output: output directory path
    :param df: dataframe path
    :return:
    """
    data_frame = f.data_reader(df)
    r = ut.sample_graph_pattern
    sub, sd = f.mean_extract(data_frame), f.sd_extract(data_frame)
    lineage = f.lineage_extract(sub, r)
    metabolites = sorted(set(sub.index))
    nrows = 2
    for l in lineage:
        flag = 0
        pdf = matplotlib.backends.backend_pdf.PdfPages(output + name + '_' + l + '_' + ut.graph_types[2])
        t_sub, s_sub = pattern_extracting(sub, sd, l)
        sns.set()
        sns.set_palette('Paired', n_colors=12, )
        fig = plt.figure(nrows, figsize=(12, 10))
        pt = 0
        while pt < len(metabolites):
            m = metabolites[pt]
            flag += 1
            if flag == 1:
                fig = plt.figure(nrows, figsize=(12, 10))
                plt.gcf().subplots_adjust(hspace=0.8, bottom=0.1, top=0.8)
            flag, pt, fig = metabo_CID(m, flag, fig, t_sub, s_sub, pdf, pt, fig)
        plt.close()
        pdf.close()


def area_metabo(metabolite, flag, fig, mean_subset, std_subset, pdf, pt):
    """
    Create histograms from the total area looping on metabolite list

    :param metabolite: current metabolite in the loop
    :param flag: the graph place in the page
    :param fig: one page
    :param mean_subset: the subset with the mean values
    :param std_subset: subset with standard deviation values
    :param pdf: pdf doc
    :param pt: iteration number
    :return:
    """
    flag += 1
    if flag == 1:
        fig = plt.figure(2, figsize=(12, 10))
        sns.set_palette('Dark2')
        plt.gcf().subplots_adjust(hspace=0.8, bottom=0.1, top=0.8)
    if flag <= 4:
        fig.add_subplot(2, 2, flag)
        m_sub = mean_subset[mean_subset.index.str.match(metabolite)]
        ms_sub = std_subset[std_subset.index.str.match(metabolite)]
        n = len(m_sub.index)

        for i in range(0, n, 1):
            plt.bar(m_sub.columns, m_sub.iloc[i], width=0.5, yerr=ms_sub.iloc[i])
            lin = np.poly1d(np.polyfit(x=[0, 1, 2], y=m_sub.iloc[i], deg=1))
            plt.plot(m_sub.columns, lin([0, 1, 2]), '--k', linestyle='--')
            plt.text(1.005, max(m_sub.iloc[i]), lin)
        plt.title(metabolite)
        plt.ylabel("Mean normalized area")

    else:
        flag = 0
        pdf.savefig()
        plt.close()
        pt -= 1
    pt += 1
    return flag, pt, fig


def area_quanti_hist(df, output, name):
    """
    Create a pdf document with total area from each metabolite at each point of time. One pdf per cell sample

    :param name: file name
    :param output: output directory path
    :param df: dataframe path
    :return:
    """
    data_frame = f.data_reader(df)
    r = ut.sample_graph_pattern
    sub, sd = f.mean_extract(data_frame), f.sd_extract(data_frame)
    lineage = f.lineage_extract(sub, r)
    metabolites = sorted(set(sub.index))
    sns.set()
    for l in lineage:
        flag = 0
        pt = 0
        t_sub, s_sub = pattern_extracting(sub, sd, l)
        pdf = matplotlib.backends.backend_pdf.PdfPages(output + name + '_' + l + ut.graph_types[1])
        fig = plt.figure(2, figsize=(12, 10))
        while pt < len(metabolites):
            m = metabolites[pt]
            flag, pt, fig = area_metabo(m, flag, fig, t_sub, s_sub, pdf, pt)
        pdf.close()
        plt.close()


def coeff_metabo(m, pt, flag, mean_subset, lineage, fig, pdf):
    """
    Create graphs from leading coefficients  looping on each metabolite

    :param m: current metabolite
    :param pt: iterator
    :param flag: graph place on the page
    :param mean_subset: dataframe with mean values
    :param lineage: cell sample
    :param fig: the page
    :param pdf: the output document
    :return:
    """
    flag += 1
    if flag == 1:
        fig = plt.figure(figsize=(8, 16))
    if flag <= 4:
        ax = fig.add_subplot(4, 1, flag)
        m_sub = mean_subset[mean_subset.index.str.match(m)]
        coeffs = []
        plt.title(m)
        for l in lineage:
            l_sub = m_sub[list(m_sub.columns[mean_subset.columns.str.match(r'.*_' + l + '_.*')])]
            lin = np.poly1d(np.polyfit(x=[0, 1, 2], y=l_sub.iloc[0], deg=1))
            coeffs.append(lin.coeffs[0])
        plt.bar(lineage, coeffs, width=0.4, color=['red', 'orange', 'yellow', 'green'])
        plt.ylabel('time evolving coefficient')
    else:
        flag = 0
        pdf.savefig(fig)
        plt.close()
        pt -= 1
    pt += 1
    return flag, pt, fig


def coeff_evolving(df, output, name):
    """
    Create PDF document with leading coefficients histograms from total area graphs

    :param name:
    :param output:
    :param df:
    :return:
    """
    data_frame = f.data_reader(df)
    r = ut.sample_graph_pattern
    sub, sd = f.mean_extract(data_frame), f.sd_extract(data_frame)
    lineage = f.lineage_extract(sub, r)
    metabolites = sorted(set(sub.index))
    flag = 0
    fig = plt.figure(figsize=(8, 16))
    sns.set()
    pt = 0
    pdf = matplotlib.backends.backend_pdf.PdfPages(output + name + ut.graph_types[3])
    while pt < len(metabolites):
        m = metabolites[pt]
        flag, pt, fig = coeff_metabo(m, pt, flag, sub, lineage, fig, pdf)
    pdf.close()
