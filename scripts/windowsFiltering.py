#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
from pathlib import Path


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument("--bed",type=str, help="bed file with the whole genome set of windows to be filtered",required=True)
parser.add_argument("--cons",type=str, help="file with consensus sequence",required=False)
parser.add_argument("--mon-len",type=int, help="satellite monomer size, if known",required=False)
parser.add_argument("--out",type=str, help="filtered windows files",required=True)

args = parser.parse_args()

input_file=args.bed
outdir=args.out
mon_len=args.mon_len
cons=Path(args.cons) if args.cons else None

Path(outdir).mkdir(parents=True, exist_ok=True)


# In[74]:


def comp_mon(unik_val, mon_ls):
    gr_multi = {}
    for base in mon_ls:
        base = int(base)
        gr_multi[base] = []
        for val in unik_val:
            val = int(val)
            for k in range(1, 6):  # 1x to 5x
                if abs(val - k * base) <= 4:
                    gr_multi[base].append(val)
                    break
    return gr_multi


# In[75]:


def main_mon(gr_data):
    sums = {}
    for k, rows in gr_data.items():
        total = sum(r[4] for r in rows)
        sums[k] = total
    return sums


# In[76]:


data = []
lengths = []
with open(input_file, "r") as f:
    for line in f:
        row = line.strip().split()
        start, end = int(row[1]), int(row[2])
        width = end - start
        row.append(width)
        data.append(row)
        lengths.append(int(row[3]))

unik = sorted(set(lengths))


# In[77]:


mon_l = []
current_cluster = [unik[0]]

for x in unik[1:]:
    if x - current_cluster[-1] <= 4:
        current_cluster.append(x)
    else:
        mon_l.append(round(sum(current_cluster) / len(current_cluster)))
        current_cluster = [x]

mon_l.append(round(sum(current_cluster) / len(current_cluster)))


# In[78]:


combined = comp_mon(unik, mon_l)
gr_data = {}
for row in data:
    length = int(row[3])
    for key, group in combined.items():
        group= [int(x) for x in group]
        if length ==int(key) or length in group:
            gr_data.setdefault(key, []).append(row)
        else:
            continue

sums_mon = main_mon(gr_data)


# In[81]:


if cons:
    with open(cons) as f:
        for i, line in enumerate(f):
            if i == 1:
                maxK = len(line.strip())
elif mon_len:
    maxK = mon_len
else:
    maxK = max(sums_mon, key=sums_mon.get)

print("Main monomer:", maxK)


# In[85]:


for key, rows in gr_data.items():
    if abs(key-maxK) <= 4 :
        outfile = f"{outdir}/mon_{key}bps_main.bed"
    else:
        outfile = f"{outdir}/mon_{key}bps.bed"

    with open(outfile, "w") as o:
        for r in rows:
            o.write("\t".join(map(str, r[:4])) + "\n")

