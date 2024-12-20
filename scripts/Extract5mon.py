#!/usr/bin/env python
# coding: utf-8

import numpy as np
import re
import random
import argparse
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",type=str, help="centromere sequence fasta file",required=True)
parser.add_argument("--bed",type=str, help="bed file with the absolute coordinates of the windows to analize",required=True)
parser.add_argument("--cons-file",type=str, help="fasta file with the consensus sequence of the satellite monomer",required=True)
parser.add_argument("--out-dir",type=str, help="full path of the directory where the output will be saved",required=True)

args = parser.parse_args()

file=args.fasta
wind_file=args.bed
consensus=args.cons_file
out=args.out_dir

def mon_extr(seq, mon_start, rev_mon_start):
    start_list=[]
    start_list_rev=[]
    mon_list=[]
    mon_list_r=[]
    for b,base in enumerate(seq):
        checkend=int(b)+6
        start2check=seq[int(b):int(checkend)]
        if start2check==mon_start:
            start_list.append(b)
        elif start2check==rev_mon_start:
            start_list_rev.append(b)
    for t,temp in enumerate(start_list):
        if t==0:
            prev=int(temp)
        else:
            mon2app=seq[int(prev):int(temp)]
            mon_list.append(mon2app)
            prev=temp
    fin_mon_list=[x for x in mon_list if len(x)==mon]
    for t_r,temp_r in enumerate(start_list_rev):
        if t_r==0:
            prev_r=int(temp_r)
        else:
            mon2app_r=seq[int(prev_r):int(temp_r)]
            mon_list_r.append(mon2app_r)
            prev_r=temp_r
    fin_mon_list_r=[x for x in mon_list_r if len(x)==mon]
    if len(fin_mon_list)<len(fin_mon_list_r):
        fin_mon_list=[y.reverse_complement() for y in fin_mon_list_r]
    if len(fin_mon_list)>5:
        random_seq=random.sample(fin_mon_list, 5)
    else:
        random_seq=[]
    return(random_seq)

with open(consensus, "r") as c:
    for b,line in enumerate(c):
        if b==0:
            continue
        elif b==1:
            cons=Seq(line)
print(cons)
mon=len(cons)
mon_start=cons[0:6]
print(mon_start)
rev_cons=cons.reverse_complement()
rev_mon_start=rev_cons[1:7]

seqs=[]
with open(wind_file, "r") as wind:
    for w in wind:
        w=w.strip().split("\t")
        #print(w)
        seq2check=[w[1],w[2]]
        seqs.append(seq2check)

mon2w=[]
full_seq=[]
with open(file, "r") as f:
    for l,line in enumerate(f):
        if l==0:
            line=line.strip()
            line=re.split(':|-|>',line)
            chrom=line[1]
            start=line[2]
            end=line[3]
        else:
            full_seq.append(line.strip())
full_seq="".join(full_seq)
#print(full_seq)
for s in seqs:
    rel_start=int(s[0])-int(start)
    rel_end=int(s[1])-int(start)
    seq=Seq(full_seq[rel_start:rel_end])
    ran_mon=mon_extr(seq, mon_start, rev_mon_start)
    print(ran_mon)
    mon2w.append(ran_mon)

with open(out, "w") as o:
    mon2w=[x for y in mon2w for x in y if len(y)>0]
    fin_ran_ls = random_seq = random.sample(mon2w, 5)
    for r, ran in enumerate(fin_ran_ls):
        o.write(">" + chrom + "_mon" + str(r + 1) + "\n")
        o.write(str(ran) + "\n")
