# -*- coding: utf-8 -*-

import argparse
import graph_functions as gr
import os

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--Areas', help='create the histograms of total area from each metabolite ')
parser.add_argument('-e', '--Enrichment', help='create histograms of isotopic marker enrichment')
parser.add_argument('-f', '--CID', help='create histograms of isotopologue fractions')
parser.add_argument('-c', '--leading_coefficient', help='create histograms with the leading coefficient from '
                                                        'the total area evolution in time '
                                                        'from each metabolite')
parser.add_argument('-o', '--output', help='output directory path', default="graphic_results/")
parser.add_argument('-n', '--name', help='file or experiment name', default='histogram')
args = parser.parse_args()


def graph_pipe():
    filename = args.name
    output = args.output

    if output == "graphic_results/":
        if not os.path.isdir(output):
            os.mkdir("./" + output)

    if args.Areas:
        gr.area_quanti_hist(args.Areas, output, filename)
    if args.Enrichment:
        gr.multi_samp_enr_hist(args.Enrichment, output, filename)
    if args.CID:
        gr.CID_hist(args.CID, output, filename)
    if args.leading_coefficient:
        gr.coeff_evolving(args.leading_coefficient, output, filename)


if __name__ == "__main__":
    graph_pipe()
