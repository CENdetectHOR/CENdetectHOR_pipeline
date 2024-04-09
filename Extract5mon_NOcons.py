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
parser.add_argument("--cons-file",type=str, help="full path of the directory where the output will be saved",required=True)

args = parser.parse_args()

file=args.fasta
wind_file=args.bed
cons=args.cons_file

def mon_extr_nocons(seq):
    s2c_list=[]
    for p,pos in enumerate(seq):
        endpos=int(p)+6
        s2c_2app=seq[p:endpos]
        s2c_list.append("".join(s2c_2app))
    for index in s2c_list:
        start_list=[]
        mon_list=[]
        for b,base in enumerate(seq):
            checkend=int(b)+6
            start2check=seq[int(b):int(checkend)]
            s2c="".join(start2check)
            if s2c==index:
                start_list.append(b)
        for t,temp in enumerate(start_list):
            if t==0:
                prev=int(temp)
            else:
                mon2app=seq[int(prev):int(temp)]
                mon_list.append(mon2app)
                prev=temp
        fin_mon_list=[x for x in mon_list if len(x)>min_mon and len(x)<max_mon]
        if len(fin_mon_list)> Nmin_mons:
            break
        else:
            continue
    fin_mon_list=[x for x in fin_mon_list if len(x)==mon]
    fin_mon_list=["".join(m) for m in fin_mon_list]
    random_seq=random.sample(fin_mon_list, 1)
    return(random_seq)

seqs=[]
with open(wind_file, "r") as wind:
    for w in wind:
        w=w.strip().split("\t")
        seq2check=[w[1],w[2]]
        seqs.append(seq2check)
        mon=int(w[3])

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

for n,s in enumerate(seqs):
    if n==0:
        len_seq=len(s)
        Nmax_mons=len_seq/mon
        Nmin_mons=Nmax_mons-(Nmax_mons/100*30)
        min_mon=mon-5
        max_mon=mon+5
        rel_start=int(s[0])-int(start)
        rel_end=int(s[1])-int(start)
        seq=Seq(full_seq[rel_start:rel_end])
        ran_mon=mon_extr_nocons(seq)
        mon2w.append(ran_mon)

with open(cons, "w") as o:
    mon2w=[x for y in mon2w for x in y if len(y)>0]
    for r,ran in enumerate(mon2w):
        o.write(">cons_"+str(r+1)+"\n")
        o.write(str(ran)+"\n")
