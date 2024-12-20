#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import multiprocessing
import math
import numpy as np
import matplotlib.pyplot as plt
from parallel_distance import build_seqs_distance_matrix, build_seqs_distance_matrix_by_chunks, FileSystemChunkStore
from Bio import SeqIO
from Bio import Phylo
from featureUtils import BED_file_to_features, feature_to_seq
from matrix_utils import save_matrix_as_triu, load_matrix_from_triu

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",type=str, help="centromere sequence fasta file",required=True)
parser.add_argument("--bed",type=str, help="bed file with the coordinates of all the monomers",required=True)
parser.add_argument("--out-distMatr",type=str, help="full path and prefix of the distance matrix files",required=True)
#parser.add_argument("--out-tree",type=str, help="name and full path of the output phyloxml file",required=True)
#parser.add_argument("--out-histPlot",type=str, help="full path and name of the histogram of distances to save",required=True)
#parser.add_argument("--out-matrPlot",type=str, help="full path and name of the matrix of distances to save",required=True)
parser.add_argument("--log",type=str, help="name and full path of the log file",required=True)
parser.add_argument("--t",type=int, help="number of threads to use",required=True)

args = parser.parse_args()

fasta=args.fasta
bed=args.bed
log=args.log
out_distMatr=args.out_distMatr
#out_tree=args.out_tree
#out_hist=args.out_histPlot
#out_matr=args.out_matrPlot
threads=args.t

old_stdout = sys.stdout
log_file = open(log,"w")
sys.stdout = log_file

multiprocessing.cpu_count()

references = {seq.id : seq for seq in SeqIO.parse(fasta, "fasta")}

monomers_as_features = BED_file_to_features(bed)

monomers_as_seqs = [feature_to_seq(feature, references) for feature in monomers_as_features]

monomer_dists_par = build_seqs_distance_matrix_by_chunks(monomers_as_seqs, num_chunks=threads, chunk_store=FileSystemChunkStore(out_distMatr+"test_{row}_{col}"))

save_matrix_as_triu(monomer_dists_par, out_distMatr+'dist_matrix.npy')

#dist_values = distance_values(monomer_dists_par)
#plt.hist(dist_values, bins=int(max(dist_values)))
#plt.savefig(out_hist)

#plt.matshow(monomer_dists_par)
#plt.savefig(out_matr)

#phyloXml, hor_tree_root, clusterings = clusterings_with_hors(monomers_as_seqs, seqs_as_features=monomers_as_features, distance_matrix=monomer_dists_par, min_len_loop=1, min_loop_reps=5)

#Phylo.write(phyloXml, out_tree, format='phyloxml')

sys.stdout = old_stdout
log_file.close()

