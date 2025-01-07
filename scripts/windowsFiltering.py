#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--bed",type=str, help="bed file with the whole genome set of windows to be filtered",required=True)
parser.add_argument("--cons",type=str, help="file with consensus sequence",required=False)
parser.add_argument("--out",type=str, help="filtered windows files",required=True)

args = parser.parse_args()

input_file=args.bed
outdir=args.out

def filter_1(data1):
    ls_cust = {}

    for line1 in data1:
        lenmon = int(line1[3])
        tname = f"mon_{lenmon}bps"

        if tname not in ls_cust:
            ls_cust[tname] = []
        ls_cust[tname].append(line1)

    return ls_cust

def comp_mon(unik_val):
    gr_multi = defaultdict(list)
    mon_comb = defaultdict(list)

    for i in unik_val:
        int_values = int(i)
        for j in unik_val:
            if j != i:
                int2comp = int(j)
                if (int2comp in range(int_values * 1 - 4, int_values * 1 + 4) or
                        int2comp in range(int_values * 2 - 4, int_values * 2 + 4) or
                        int2comp in range(int_values * 3 - 4, int_values * 3 + 4) or
                        int2comp in range(int_values * 4 - 4, int_values * 4 + 4) or
                        int2comp in range(int_values * 5 - 4, int_values * 5 + 4)):
                    gr_multi[int_values].append(int2comp)

    for key, multiples in gr_multi.items():
        #print(f"mon_{key}:", multiples)
        mon_comb[key].extend(multiples)
    return mon_comb

def main_mon(monomers):
    sums_mon = defaultdict(int)
    for k,v in monomers.items():
        sumlen=0
        vals=list(v)
        for record in vals:
            val2sum=record[4]
            sumlen=sumlen+val2sum
        sums_mon[k]=sumlen

    return sums_mon

def out_files(data1, mon_comb, output_dir):
    gr_data = defaultdict(list)
    for index, row in enumerate(data1):
        lenmon = int(row[3])
        for key, group in mon_comb.items():
            if lenmon == key or lenmon in group:
                group_key = f"mon_{key}bps"
                gr_data[group_key].append(row)

                break
    sums_mon = main_mon(gr_data)
    max_length = max(sums_mon.values())
    if args.cons is None:
        maxK=list(sums_mon.keys())[list(sums_mon.values()).index(max_length)]
    else:
        with open(args.cons, "r") as cons:
            for l,line in enumerate(cons):
                if l==1:
                    maxK=len(line.strip().split())
    for group_key, rows in gr_data.items():
        output_file = f"{output_dir}/{group_key}.bed"
        if group_key==maxK:
            output_file = f"{output_dir}/{group_key}_main.bed"
        with open(output_file, 'w') as o:
            for row in rows:
                row=row[0:4]
                o.write(("\t").join(row)+"\n")

data=[]
with open(input_file, "r") as f:
    for line in f:
        line=line.strip().split()
        start=line[1]
        end=line[2]
        width=int(end)-int(start)
        line.append(width)
        data.append(line)

len_mon_full=[int(x[3])for x in data]
unik = list(pd.Series(len_mon_full).unique())
filtered = filter_1(data)
combined = comp_mon(unik)
out_files(data, combined, outdir)
