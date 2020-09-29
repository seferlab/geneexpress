import os
import sys
import math

fname1 = "127 LCM time course Data Not normalized.txt"
fname2 = "127 LCM time course Quantile Normalized logbased 2 transformed.txt"

with open(fname1,"r") as infile:  
    for line in infile:
        line = line.rstrip()
        vals = line.split("\r")
        splitted = vals[1].split("\t")
        items1 = [float(splitted[tind]) for tind in xrange(1,len(splitted))]

with open(fname2,"r") as infile:  
    for line in infile:
        line = line.rstrip()
        vals = line.split("\r")
        splitted = vals[1].split("\t")
        items2 = [float(splitted[tind]) for tind in xrange(1,len(splitted))]
              
print items1[0:20]
print [math.log(titem,2.0) for titem in items1[0:10]]
print [math.log(titem+1.0,2.0) for titem in items1[0:10]]
print items2[0:20]
print items1[8:20]

