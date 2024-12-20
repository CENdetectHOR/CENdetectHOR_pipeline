import argparse
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import Phylo
from featureUtils import BED_file_to_features, feature_to_seq
from Bio.Phylo.PhyloXML import Phyloxml
from hierarchical_clustering import hierarchical_clustering, build_phylogeny, get_matrices, save_phylogeny, save_matrices, load_phylogeny, load_matrices, get_clustering_matrices
from matrix_utils import load_matrix_from_triu

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",type=str, help="centromere sequence fasta file",required=True)
parser.add_argument("--bed",type=str, help="bed file with the coordinates of all the monomers",required=True)
parser.add_argument("--distMatr",type=str, help="full path of the distance matrix file",required=True)
parser.add_argument("--out-tree",type=str, help="name and full path of the output phyloxml file",required=True)

args = parser.parse_args()

fasta=args.fasta
bed=args.bed
distMatr=args.distMatr
out_tree=args.out_tree

references = {seq.id : seq for seq in SeqIO.parse(fasta, "fasta")}
monomers_as_features = BED_file_to_features(bed)

monomers_as_seqs = [feature_to_seq(feature, references) for feature in monomers_as_features]

dist_matrix = load_matrix_from_triu(distMatr)

phylogeny = hierarchical_clustering(dist_matrix=dist_matrix)
save_phylogeny(phylogeny=phylogeny, filename=out_tree)

phyloXml = Phyloxml(phylogenies=[build_phylogeny(simple_phylogeny=phylogeny, features=monomers_as_features)], attributes=None)
Phylo.write(phyloXml, out_tree+'.xml', format='phyloxml')