#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from cluster import distance_values
from Bio import SeqIO
from Bio import Phylo
from showHOR import show_hor, show_hors, show_hor_tree
from featureUtils import BED_file_to_features, feature_to_seq, remove_overlapping_features
from parallel_distance import build_seqs_distance_matrix_by_chunks, FileSystemChunkStore
from hor_tree import phylogeny_to_hor_tree
from Bio.Phylo.PhyloXML import Phyloxml
from Bio.Phylo import PhyloXMLIO
from clustering_to_phylogeny import clustering_to_phylogeny
from mixed_direction_hors import find_inversion_loops
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",type=str, help="centromere sequence fasta file",required=True)
parser.add_argument("--bed",type=str, help="bed file with the coordinates of all the monomers",required=True)
parser.add_argument("--out-distMatr",type=str, help="full path and prefix of the distance matrix files",required=True)
parser.add_argument("--out-tree",type=str, help="name and full path of folder where the output phyloxml trees will be saved",required=True)
#parser.add_argument("--out-histPlot",type=str, help="full path and name of the histogram of distances to save",required=True)
#parser.add_argument("--out-matrPlot",type=str, help="full path and name of the matrix of distances to save",required=True)
parser.add_argument("--log",type=str, help="name and full path of the log file",required=True)
parser.add_argument("--t",type=int, help="number of threads to use",required=True)

args = parser.parse_args()

fasta=args.fasta
bed=args.bed
log=args.log
out_distMatr=args.out_distMatr
out_tree=args.out_tree
#out_hist=args.out_histPlot
#out_matr=args.out_matrPlot
threads=args.t

old_stdout = sys.stdout
log_file = open(log,"w")
sys.stdout = log_file

references = {seq.id : seq for seq in SeqIO.parse(fasta, "fasta")}

monomers_as_features = BED_file_to_features(bed)

monomers_as_features = remove_overlapping_features(
    features=monomers_as_features,
    expected_feature_size=171,
    max_allowed_overlap_fraction=0.25
)

monomers_as_seqs = [feature_to_seq(feature, references) for feature in monomers_as_features]

monomer_dists = build_seqs_distance_matrix_by_chunks(monomers_as_seqs, num_chunks=t, chunk_store=FileSystemChunkStore(out_distMatr+"matr_{row}_{col}"))

with open(out_distMatr+'dist_matrix.npy', 'wb') as f:
    np.save(f, monomer_dists)

curr_recursion_limit=sys.getrecursionlimit()
if (curr_recursion_limit<len(monomers_as_features)):
    sys.setrecursionlimit(len(monomers_as_features))

dist_values = distance_values(monomer_dists)
#plt.hist(dist_values, bins=int(max(dist_values)))

clustering_res = clustering_to_phylogeny(
    dist_matrix=monomer_dists,
    items_as_seq_features=monomers_as_features,
    seq_references=references
)
phylogeny = clustering_res.phylogeny

hor_tree = phylogeny_to_hor_tree(phylogeny, min_loops=4, allow_hor_overlap=False, discrete_sorted_levels=True)

inversion_loops = find_inversion_loops(seq_features=monomers_as_features, min_loops=4)
[str(loop_inSeq) for loop_inSeq in inversion_loops]

phyloXml = Phyloxml(phylogenies=[phylogeny, hor_tree.as_phyloxml], attributes=None)
Phylo.write(phyloXml, out_tree, format='phyloxml')

sys.stdout = old_stdout
log_file.close()
