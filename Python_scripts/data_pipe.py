# -*- coding: utf-8 -*-

import data_pipe_functions as f
import pandas as pd
import re
import _utils as ut
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-df', '--Data_frame', type=str, help='\nDataset already corrected by IsoCor', required=True)
parser.add_argument('-o', '--output', type=str, help='\nOutput directory path', default='./data_frame_results/')
parser.add_argument('-n', '--name', type=str, help='\nName your result table', default='result')
parser.add_argument('-a', '--Areas', help='\nOutput the result table of the averaged sum of Areas ',
                    action='store_true')
parser.add_argument('-i', '--Isotopologue', help='\nOutput the result table of isotopologue fractions ',
                    action='store_true')
parser.add_argument('-e', '--Enrichment', help='\nOutput the result table of isotope enrichment', action='store_true')
parser.add_argument('-c', '--Cell_list', help='\ninput file with the count of cells used to make the measures',
                    default=None)
parser.add_argument('-m', '--marked_only', help='\n remove all non marked isotopologue', action='store_true')

args = parser.parse_args()


def table_results(df):
    # useful var

    global iso_mean, enr_mean, area_mean, cellnumber
    filenames = ['_CID.csv', '_enrichement.csv', '_area_means.csv']
    data_frame = f.data_reader(df)
    if args.marked_only:
        data_frame = data_frame[data_frame['isotopologue'] != 0]

    if args.output == './data_frame_results/':
        if not os.path.isdir(args.output):
            os.mkdir('./data_frame_results/')

    p = re.compile(r"(T\d*_\w*)_(.*)")  # sample recognition pattern
    samples = sorted(set(data_frame['sample']))  # samples list
    metabolites = sorted(set(data_frame['metabolite']))  # metabolites list
    if args.Cell_list:
        cellnumber = f.cell_list(args.Cell_list)
    else:
        cellnumber = None
    history = []
    df_create = False

    for i in samples:
        m = p.match(i)  # identify samples

        if m is not None and m.group(1) not in history:
            samp = m.group(1)
            history.append(samp)
            tmp = f.sample_extract(data_frame, samp)  # create subset from each sample
            if samp in ut.a_sample_list:  # avoiding AB samples with A samples
                tmp = tmp[tmp['sample'].str.match(ut.A_B_sample)]

            # creating empty dataframes
            if not df_create and tmp is not None:
                if args.Isotopologue is None and args.Areas is None and args.Enrichment is None:
                    return print("\n You need to choose an output, please type -h or --help for more information")

                    # transitional dataframe to create isotopologue  dataframe
                ind = tmp[tmp['sample'].str.match(ut.replicate_1)]
                iso_mean = pd.DataFrame(index=sorted(ind['metabolite']))  # isotopologue dataframe creation

                # create the two others dataframes (enrichment and mean Areas)
                enr_mean = pd.DataFrame(index=metabolites)
                area_mean = pd.DataFrame(index=metabolites)
                df_create = True

            if df_create:
                a_mean_list, e_mean_list, i_mean_list = [], [], []
                a_sd_list, e_sd_list, i_sd_list = [], [], []

                # creating the dataframe columns
                for u in metabolites:
                    m_t = f.metabolite_extract(tmp, u)
                    a_mean_list, a_sd_list = f.area_lists(m_t, a_mean_list, a_sd_list)
                    e_mean_list, e_sd_list = f.enrichment_lists(m_t, e_mean_list, e_sd_list)
                    i_mean_list, i_sd_list = f.iso_lists(m_t, i_mean_list, i_sd_list)

                # fulfill the dataframes
                iso_mean = f.iso_means(i_mean_list, i_sd_list, iso_mean, samp)
                enr_mean = f.enrichment_means(e_mean_list, e_sd_list, enr_mean, samp)
                area_mean = f.area_means(a_mean_list, a_sd_list, area_mean, samp, cellnumber)

    #  write the output files
    if args.marked_only:
        name = args.name + '_marked_only'
    else:
        name = args.name
    if args.Isotopologue:
        iso_mean.to_csv(args.output + name + filenames[0])
    if args.Enrichment:
        enr_mean.to_csv(args.output + name + filenames[1])
    if args.Areas:
        area_mean.to_csv(args.output + name + filenames[2])


if __name__ == "__main__":
    table_results(args.Data_frame)
