#!/usr/bin/env python
# coding: utf-8

from itertools import groupby
from operator import itemgetter
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import re
import argparse

def make_mers(seq,mer_len):
    d={}
    for w,win in enumerate(seq):
        for n,base in enumerate(win):
            kmer=win[n:(n+mer_len)]
            kmer="".join(kmer)
            if n==0:
                d[w]=[str(kmer.strip())]
            else:
                if str(kmer.strip()) not in d[w] and len(kmer)==mer_len:
                    d[w].append(kmer)
    return(d)

def dist_calc(d, seq):
    finlist=[]
    d2plot=[]
    diff=[]
    for k,key in enumerate(d):
        kmer_list=d[k]
        string=seq[k]
        for mer in kmer_list:
            ilist=[]
            len_site = len(mer)
            if mer in string:
                for i in range (0, len(string)):
                    if mer == string[i:i+len_site]:
                        ilist.append(i)
                    elif mer not in string:
                        print('No such a site in the string')
            finlist.append(ilist)
    for l in finlist:
        if len(l)>1:
            diff= [l[i+1] - l[i] for i in range(len(l)-1)]
            d2plot.append(diff)
    d2plot=[item for sublist in d2plot for item in sublist]
    return(d2plot)


parser = argparse.ArgumentParser()
parser.add_argument("--fasta",type=str, help="centromere sequence fasta file",required=True)
parser.add_argument("--window-size", type=int, help="splitting window size",required=False, default=1000)
parser.add_argument("--kmer-size", type=int, help="size of the kmer to search for periodicity",required=False, default=8)
parser.add_argument("--plot-fold",type=str, help="full path of the folder in which the output plots will be saved",required=True)
parser.add_argument("--wind-summ",type=str, help="output windows summary file",required=True)
parser.add_argument("--wind-bed",type=str, help="bed file with the windows absolute coordinates",required=True)

args = parser.parse_args()

file=args.fasta
windSize=args.window_size
kmer=args.kmer_size
outfold=args.plot_fold
windSumm=args.wind_summ
bed=args.wind_bed

seq=[]
with open(file, "r") as f:
    for l,line in enumerate(f):
        if l== 0:
            line = line.strip()
            line = re.split(':|-|>', line)
            chrom = line[1]
            start = int(line[2])
            print(start)
            end = line[3]
            print(line)
        else:
            seq.append(line.strip())
    seq=[x for y in seq for x in y]
    #dividing sequence in windows of 1kb 
    seq=[seq[i:i + windSize] for i in range(0, len(seq), windSize)]

d=make_mers(seq, kmer)

ind4merge=[1 if windSize-len(d[k])>100 else 0 for k in d]

seq=[''.join(x) for x in seq]

result=[]
for key, group in groupby(zip(ind4merge, seq), itemgetter(0)):
    if key ==1:
        l2app =[c for i, c in group]
        result.append(''.join(l2app))
    else:
        l2app=[c for i, c in group]
        result=result+l2app

with open(windSumm, "w") as o:
    lenlist=[]
    out2w=[]
    for r,res in enumerate(result):
        if r==0:
            rel_start=0
            l=len(res)
            rel_end=rel_start+l
            abs_start=int(start)+rel_start
            abs_end=int(start)+rel_end
            l2w=[str(abs_start),str(abs_end),str(l)]
            lenlist.append(l)
            out2w.append(l2w)
        else:
            prevline=r-1
            rel_start=int(out2w[prevline][1])
            l=len(res)
            rel_end=rel_start+l
            l2w = [str(rel_start), str(rel_end), str(l)]
            lenlist.append(l)
            out2w.append(l2w)
    [o.write('\t'.join(x)+'\n') for x in out2w]

#ind_list=[lenlist.index(x) for x in lenlist if x>5000]
ind_list=[]
thres=5*windSize
for v,val in enumerate(lenlist):
    if val>thres and v not in ind_list:
        ind_list.append(v)

l2bed_list=[]
for x in ind_list:
    new_seq=result[x]
    new_seq = [new_seq[i:i + windSize] for i in range(0, len(new_seq), windSize)]
    mer_dict=make_mers(new_seq, kmer)
    d2plot=dist_calc(mer_dict, new_seq)
    count_dist=Counter(d2plot)
    monlen=int(list(count_dist.most_common(1)[0])[0])
    plt.figure(figsize=(20,15))
    plt.hist(d2plot, bins = 963)
    plt.xticks(np.arange(0, 1000, 100))
    plt.savefig(outfold+"/"+chrom+"_"+str(out2w[x][0])+"-"+str(out2w[x][1])+".pdf")
    l2bed=[chrom, str(out2w[x][0]), str(out2w[x][1]), str(monlen)]
    l2bed_list.append(l2bed)
with open(bed, "w") as outbed:
    [outbed.write('\t'.join(x)+'\n') for x in l2bed_list]

