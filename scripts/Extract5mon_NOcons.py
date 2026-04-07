#!/usr/bin/env python
# coding: utf-8

import numpy as np
import re
import random
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--fasta",type=str, help="centromere sequence fasta file",required=True)
parser.add_argument("--bed",type=str, help="bed file with the absolute coordinates of the windows to analize",required=True)
parser.add_argument("--out-cons-file",type=str, help="full path of the directory where the output will be saved",required=True)
parser.add_argument("--mon-len",type=int, help="satellite monomer size, if known",required=False)

args = parser.parse_args()

file=args.fasta
wind_file=args.bed
cons=args.out_cons_file

def mon_extr_nocons(seq, min_mon, max_mon, Nmin_mons, monLen,k):
    seq = str(seq)
    kmer_pos = defaultdict(list)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_pos[kmer].append(i)
    for kmer, positions in kmer_pos.items():
        if len(positions) < 2:
            continue
        mon_list = []
        prev = positions[0]
        for pos in positions[1:]:
            mon_len = pos - prev
            if min_mon < mon_len < max_mon:
                mon_list.append(seq[prev:pos])
            prev = pos
        if len(mon_list) > int(Nmin_mons):
            fin_mon_list = [m for m in mon_list if len(m) == monLen]
            return fin_mon_list

seqs=[]
mon=[]
with open(wind_file, "r") as wind:
    for w in wind:
        w=w.strip().split("\t")
        seq2check=[w[1],w[2]]
        seqs.append(seq2check)
        mon.append(w[3])

if args.mon_len:
    monLen = args.mon_len
else:
    monLen=int(max(set(mon), key=mon.count))

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

for k in range(6,3,-1):
    mon2w = []
    for n,s in enumerate(seqs):
        len_seq=int(s[1])-int(s[0])
        Nmax_mons=len_seq/monLen
        Nmin_mons=Nmax_mons-(Nmax_mons/100*30)
        min_mon=monLen-5
        max_mon=monLen +5
        rel_start=int(s[0])-int(start)
        rel_end=int(s[1])-int(start)
        seq = full_seq[rel_start:rel_end]
        ran_mon=mon_extr_nocons(seq,min_mon,max_mon,Nmin_mons,monLen, k)
        #mon2w.append(ran_mon)
        if ran_mon:
            random_seq=random.sample(ran_mon, 1)
            mon2w.append(random_seq[0])
            break
    if len(mon2w) >= 1:
        break

with open(cons, "w") as o:
    #mon2w=[x for y in mon2w for x in y if len(y)>0]
    for r,ran in enumerate(mon2w):
        o.write(">cons_"+str(r+1)+"\n")
        o.write(str(ran))
