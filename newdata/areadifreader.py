import sys
import os
import numpy as np
import math
import scipy as sp
import scipy.stats

gene2val = {}
vals = []
fpath = "areadif/select3.txt"
with open(fpath,"r") as infile:
    for line in infile:
        line = line.rstrip()
        tgene,tval = line.split("\t")
        tval = float(tval)
        gene2val[tgene] = tval
        vals.append(tval)
print np.mean(vals)

count = 0
for tval in vals:
    if tval >= 1.0:
       count +=1
print count        

#print scipy.stats.mstats.normaltest(np.array(vals))
print scipy.stats.shapiro(vals)     
        

