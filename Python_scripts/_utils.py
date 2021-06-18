# -*- coding: utf-8 -*-
import re

# the amino acid list
aa = ["Glutamate", "Alanine", "Aspartate", "Glutamine", "Glycine", "Serine", "Alanine", "Aspartic acid",
      "Arginine", "Cysteine", "Glutamic acid", "Histidine", "Isoleucine", "Lysine", "Leucine",
      "Phenylalanine", "Methionine", "Proline", "Tryptophan", "Threonine", "Tyrosine", "Valine"]

sample_graph_pattern = re.compile(r"(T\d*)_(\w*)_.*")

A_B_sample = re.compile(r"(T\d*_\w)_.*")  # identify A or B samples (not AB)

replicate_1 = re.compile(r"(T\d*_\w*)_1")  # identify one replicate

a_sample_list = ["T0_A", "T24_A", "T48_A"]  # A mutant sample list

graph_types = ['_enrichment.pdf', '_total_area.pdf', '_CID.pdf', '_leading_coefficient.pdf']