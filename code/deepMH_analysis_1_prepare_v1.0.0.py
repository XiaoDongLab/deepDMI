### deepMH - identify methylation heterogeneity from single-cell bisulfite sequencing using deep learning
# Copyright (C) 2022  Dong, et. al.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

### updates
# v1.0.0 - 2022.10.06 - initial version.

# in python3

import argparse
import sys
import string
import numpy as np
import datetime
import dill
import multiprocessing

version="v1.0.0"
date="2022.10.06"

###
parser=argparse.ArgumentParser(description="deepMH_analysis_1_prepare_" + version + ".py\nPrepare inputdata.\nDate: " + date + "; Version: " + version)
parser.add_argument("-i","--input",type=str, required=True, help="Input header of bed file")
parser.add_argument("-b","--num_bp",type=int, default=100, help="No. upstream and downstream basepairs to consider. Default=100")
parser.add_argument("-c","--num_cpg",type=int, default=10, help="No. flanking CpGs to consider. Default=10")
parser.add_argument("-r","--ref_genome",type=str, required=True, help="Path to refgenome.fa file")
parser.add_argument("-t","--num_cpu",type=int, default=24, help="No. CPUs to run in parallele. Default=24")
args=parser.parse_args()

### notes ###
input_refgenome=args.ref_genome
input_header=args.input
input_nbp=args.num_bp
input_ncpg=args.num_cpg
input_ncpu=args.num_cpu

### import data ###
genome_in = open(input_refgenome, 'r')
meth_in = open(input_header + '.bed', 'r')
dillfile = input_header + ".pkl"

### preprocessing the reference genome ###
genome = genome_in.read()
genome_in.close()
chromo = genome.split('>')

del genome

chromo_index = list()
i = 1
while i < len(chromo):
    tmp = chromo[i][0:200]
    chromo_index.append((tmp.split('\n')[0]).split(' ')[0])
    i = i + 1

basepair = list()
i = 1
while i < len(chromo):
    tmp = chromo[i]
    tmp = tmp.split('\n')
    tmp[0] = '0'
    tmp1 = str()
    for element in tmp:
        tmp1 = tmp1 + element
    basepair.append(tmp1)
    i = i + 1

del chromo

### function returns the neighboring nucleotide ###
# example of using it, in python
# triple('chrX',123)
# it would return a list, i.g. ['A','T','C'], corresponding to -1, 0, +1 position of your mutation
# note, 'chrX' should be provided in a proper way, you can check the variable 'chromo_index' of how to write it.
def triple(a_chr, b_pos, c_start=-input_nbp, d_end=input_nbp):
    if int(b_pos) + int(c_start) < 1:
        return "pos.exceedchrlength"
    a_chr = str(a_chr) # a_chr, string
    #if 'e+' in b_pos:
    #       b_pos = float(b_pos.split('e+')[0])*10**float(b_pos.split('e+')[1])
    b_pos = int(b_pos) # b_pos, int
    i=0
    pick = -1
    while i < len(chromo_index):
        if a_chr == chromo_index[i]:
            pick = i
            break
        i = i + 1
    if pick == -1:
        return "chr.notmatching"
    if b_pos + d_end + 1 >= len(basepair[pick]):
        return "pos.exceedchrlength"
    #result_triple = [basepair[i][b_pos-1], basepair[i][b_pos], basepair[i][b_pos+1]]
    result_triple = []
    j = c_start
    while j <= d_end:
        result_triple.append(basepair[pick][b_pos+j])
        j = j + 1
    return result_triple

### preprocessing methylation bed file ###
meth_data=[] 
while True:
    line=meth_in.readline()
    if len(line)==0:
        break
    line=line.split('\n')[0]
    if len(line)==0:
        break
    line=line.split('\t')
    meth_data.append(line)

meth_in.close()

meth_chr=[]
for e in meth_data:
    if e[0] in meth_chr:
        continue
    meth_chr.append(e[0])

meth_pos=[]
meth_val=[]
for e in meth_chr:
    meth_pos.append([])
    meth_val.append([])

for e in meth_data:
    meth_pos[meth_chr.index(e[0])].append(float(e[2]))
    meth_val[meth_chr.index(e[0])].append(e[3])

def triple_meth(a_chr, b_pos, c_start=-input_nbp, d_end=input_nbp, na_value="NA"):
    a_chr=str(a_chr)
    b_pos=int(b_pos)
    f1=[]
    f2=[]
    meth_datatmp=meth_pos[meth_chr.index(a_chr)]
    for e in range(0, len(meth_datatmp)):
        if int(meth_datatmp[e]) >= b_pos + c_start and int(meth_datatmp[e]) <= b_pos + d_end and meth_val[meth_chr.index(a_chr)][e] != "NA":
            f1.append(int(meth_datatmp[e]) - b_pos)
            f2.append(float(meth_val[meth_chr.index(a_chr)][e]))
    return([f1,f2])

##### test with two layers of data as input to model #####
# min CpG in bin: 5
def tf_format(e, tm, min_cpg=input_ncpg, binsize=input_nbp):
    tmp=tm
    if len(tmp[0]) <= min_cpg:
        return(np.array([0]))
    dd=np.array([0 for i in range(0,binsize*2+1)], dtype=np.float)
    dd.fill(np.nan)
    for f in range(0,len(tmp[0])):
        if tmp[0][f]!=0:
            dd[tmp[0][f] + binsize] = tmp[1][f]
    tmp2=triple(e[0], e[2])
    dd2=[]
    for f in tmp2:
        if f=="A":
            dd2.append(1)
        elif f=="C":
            dd2.append(2)
        elif f=="G":
            dd2.append(3)
        elif f=="T":
            dd2.append(4)
        else:
            dd2.append(0)
    dd2=np.array(dd2, dtype=np.float)
    return(np.array([dd, dd2], dtype="float32"))

def tf_meth(tm, min_cpg=input_ncpg):
    if len(tm[0]) <= min_cpg:
        return(np.array([0]))
    tm0=[abs(i) for i in tm[0]]
    tm0.sort()
    tm1=[]
    for i in range(1, min_cpg):
        for j in range(0, len(tm[0])):
            if tm0[i] == abs(tm[0][j]):
                tm1.append(tm[1][j]) # if -1, 0, 1 and min_cpg=1; I only selected -1 position
    return([tm0[1:min_cpg], tm1[0:(min_cpg-1)]])

def main(ele):
    if ele[3]=="NA":
        return(None, None, None, None, None)
    yy = triple_meth(ele[0],ele[2])
    if len(yy[0])<=min_cpg:
        return(None, None, None, None, None)
    xx = tf_format(ele, tm=yy)
    if len(xx) != 2:
        return(None, None, None, None, None)
    return(float(ele[3]), xx, np.array(tf_meth(yy), dtype="float32"), ele[0], ele[2])


min_cpg=input_ncpg
pool = multiprocessing.Pool(input_ncpu)
print(print(datetime.datetime.now()))
a1, a2, a3, a4, a5 = zip(*pool.map(main, meth_data))

print(print(datetime.datetime.now()))

del pool

a1=list(a1)
data_label = [i for i in a1 if type(i) != type(None)]

a2=list(a2)
data_feature = [i for i in a2 if type(i) != type(None)]

a3=list(a3)
data_meth = [i for i in a3 if type(i) != type(None)]

a4=list(a4)
data_label_chr = [i for i in a4 if type(i) != type(None)]

a5=list(a5)
data_label_pos = [i for i in a5 if type(i) != type(None)]

del a1, a2, a3, a4, a5

data_feature=np.array(data_feature)
data_label=np.array(data_label)
data_meth=np.array(data_meth)

del basepair

dill.dump_session(dillfile)
