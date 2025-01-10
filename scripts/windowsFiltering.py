#!/usr/bin/env python
# coding: utf-8

#import pandas as pd
#from collections import defaultdict
import pathlib
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--bed",type=str, help="bed file with the whole genome set of windows to be filtered",required=True)
parser.add_argument("--cons",type=str, help="file with consensus sequence",required=False)
parser.add_argument("--out",type=str, help="filtered windows files",required=True)

args = parser.parse_args()

input_file=args.bed
outdir=args.out
cons=Path(args.cons)

pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

def comp_mon(unik_val, mon_ls):
    gr_multi = {}
    mon_comb = {}

    for i in mon_ls:
        int_values = int(i)
        #print(i)
        for j in unik_val:
            #if j != i:
                #print(j)
             int2comp = int(j)
             if (int2comp in range(int_values * 1 - 4, int_values * 1 + 4) or
                 int2comp in range(int_values * 2 - 4, int_values * 2 + 4) or
                 int2comp in range(int_values * 3 - 4, int_values * 3 + 4) or
                 int2comp in range(int_values * 4 - 4, int_values * 4 + 4) or
                 int2comp in range(int_values * 5 - 4, int_values * 5 + 4)):
                 #print(int2comp)
                 [gr_multi[i].append(int2comp) if i in gr_multi else gr_multi.update({i:[int2comp]})]
     #   mon_comb[key].extend(multiples)
    return gr_multi

def main_mon(monomers):
    sums_mon = {}
    for k,v in monomers.items():
        sumlen=0
        vals=list(v)
        for record in vals:
            val2sum=int(record[4])
            sumlen=sumlen+val2sum
        sums_mon[k]=sumlen

    return sums_mon

def out_files(data1, mon_comb, output_dir,cons):
    gr_data = {}
    for index, row in enumerate(data1):
        lenmon = int(row[3])
        for key, group in mon_comb.items():
 #           if type(group) is list:
 #               if lenmon == key or lenmon in group:
 #                   group_key = key
 #                   [gr_data[group_key].append(row) if group_key in gr_data else gr_data.update({group_key:[row]})]
 #           else:
                if lenmon in group:
                    group_key = key
                    [gr_data[group_key].append(row) if group_key in gr_data else gr_data.update({group_key:[row]})]

  #  print(gr_data)
    sums_mon = main_mon(gr_data)
  # print(sums_mon)
    max_length = max(sums_mon.values())
    if not cons.exists():
        maxK=list(sums_mon.keys())[list(sums_mon.values()).index(max_length)]
        print(maxK)
    else:
        with open(cons, "r") as cons:
            for l,line in enumerate(cons):
                if l==1:
                    #line=line.strip().split()
                    line_split=[ch for ch in line]
                    maxK=len(line_split)
                    print(maxK)
    for group_key, rows in gr_data.items():
        output_file = f"{output_dir}/mon_{group_key}bps.bed"
        if group_key in range(maxK-4,maxK+4):
            print(group_key)
            output_file = f"{output_dir}/mon_{group_key}bps_main.bed"
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
unik = list(set(len_mon_full))
print(unik)
mon_l=[]
mon2c=[]
for m,mon in enumerate(unik):
    #if m ==0:
     #   mon2c.append(mon)
   # else:
    for t in unik:
        if mon in range(t-4,t+4):
            min_mon=min(mon,t)
            print(t)
            print(mon)
            print(min_mon)
            if min_mon not in mon_l:
                mon_l.append(min_mon)
                mon2c.append(min_mon)
        else:
            if mon not in mon_l:
                mon_l.append(mon)
                mon2c.append(mon)
print(mon2c)
print(mon_l)
combined = comp_mon(unik,mon_l)
print(combined)
out_files(data, combined, outdir,cons)
