#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from iterativeClustering import clusterings_with_hors
from cluster import build_seqs_distance_matrix, distance_values
from Bio import SeqIO
from Bio import Phylo
from featureUtils import BED_file_to_features, feature_to_seq

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",type=str, help="centromere sequence fasta file",required=True)
parser.add_argument("--bed",type=str, help="bed file with the coordinates of all the monomers",required=True)
parser.add_argument("--out-tree",type=str, help="name and full path of the output phyloxml file",required=True)
parser.add_argument("--out-distHist",type=str, help="full path and name of the histogram of distances to save",required=True)
parser.add_argument("--out-distMatr",type=str, help="full path and name of the matrix of distances to save",required=True)
parser.add_argument("--log",type=str, help="name and full path of the log file",required=True)

args = parser.parse_args()

fasta=args.fasta
bed=args.bed
log=args.log
out_tree=args.out_tree
out_hist=args.out_distHist
out_matr=args.out_distMatr

old_stdout = sys.stdout
log_file = open(log,"w")
sys.stdout = log_file

references = {seq.id : seq for seq in SeqIO.parse(fasta, "fasta")}

monomers_as_features = BED_file_to_features(bed)

monomers_as_features[0]

[(feature.location.start, feature.__len__()) for feature in monomers_as_features if feature.__len__() < 70]

monomers_as_seqs = [feature_to_seq(feature, references) for feature in monomers_as_features]

monomer_dists = build_seqs_distance_matrix(monomers_as_seqs)

def max_len(strings):
    return max([len(s) for s in strings])

dist_values = distance_values(monomer_dists)
plt.hist(dist_values, bins=int(max(dist_values)))
plt.savefig(out_hist)

plt.matshow(monomer_dists)
plt.savefig(out_matr)

phyloXml, hor_tree_root, clusterings = clusterings_with_hors(monomers_as_seqs, seqs_as_features=monomers_as_features, distance_matrix=monomer_dists, min_len_loop=1, min_loop_reps=5)

Phylo.write(phyloXml, out_tree, format='phyloxml')

sys.stdout = old_stdout
log_file.close()
