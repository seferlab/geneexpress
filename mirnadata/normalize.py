import numpy as np
import scipy as sp
import sys
import os
import math
from copy import deepcopy

def normData(data,usetimes,normtype):
    """variance normalization
    Args:
       data:
       usetimes:
       normtype:
    Returns:
       moddata:   
    """
    assert normtype in ["all","pertime"]
    moddata = [{} for rind in xrange(len(data))]
    if normtype == "pertime":
       for ttime in usetimes:
           arrcount = len(data[0][ttime])
           Mdict,Mvaldict = {}, {}
           for arrind in xrange(arrcount):
               Mdict[arrind] = np.mean([rowlist[ttime][arrind] for rowlist in data])
               Mvaldict[arrind] = [rowlist[ttime][arrind] for rowlist in data]
           M = np.mean(Mdict.values())
           Vdict = {}
           for arrind in xrange(arrcount):
               Vdict[arrind] = np.var(Mvaldict[arrind])    
           V = np.mean(Vdict.values())        
           for rind,rowlist in enumerate(data):
               outs = []
               for arrind,titem in enumerate(rowlist[ttime]):
                   outs.append(float(titem-Mdict[arrind]+M)*math.sqrt(V)/math.sqrt(Vdict[arrind]))
               moddata[rind][ttime] = list(outs)
    elif normtype == "all":
       Mdict,Mvaldict = {}, {} 
       for ttime in usetimes:
           arrcount = len(data[0][ttime])
           Mdict[ttime], Mvaldict[ttime] = {},{}
           arrcount = len(data[0][ttime])
           for arrind in xrange(arrcount):
               Mdict[ttime][arrind] = np.mean([rowlist[ttime][arrind] for rowlist in data])
               Mvaldict[ttime][arrind] = [rowlist[ttime][arrind] for rowlist in data]
       M = np.mean([meanval for ttime in Mdict.keys() for meanval in Mdict[ttime].values()])
       Vdict = {}
       for ttime in usetimes:
           arrcount = len(data[0][ttime])
           Vdict[ttime] = {}
           for arrind in xrange(arrcount):
               Vdict[ttime][arrind] = np.var(Mvaldict[ttime][arrind])    
       V = np.mean([varval for ttime in Vdict.keys() for varval in Vdict[ttime].values()])       
       for rind,rowlist in enumerate(data):
           for ttime in usetimes:
               outs = []
               for arrind,titem in enumerate(rowlist[ttime]):
                   outs.append(float(titem-Mdict[ttime][arrind]+M)*math.sqrt(V)/math.sqrt(Vdict[ttime][arrind]))
               moddata[rind][ttime] = list(outs)         
    return moddata


def writeMiRNAData(outfname,gene2data,allgenes,timekeys):
    """writes mirna data
    Args:
       outfname,gene2data,allgenes:
       timekeys:
    Returns:
    """    
    with open(outfname,"w") as outfile:
       outfile.write("Gene\t{0}\n".format("\t".join(timekeys)))
       for gind in xrange(len(allgenes)):
           outfile.write("{0}\t{1}\n".format(allgenes[gind],"\t".join([str(titem) for titem in gene2data[gind]])))

           
def logTransform(mydata,gene2data,usetimes):
    """log transform data
    Args:
       mydata:
       gene2data:
       usetimes:
    Returns:
       loggene2data:
       logmydata:   
    """
    loggene2data,logmydata = [],[]
    for rowlist in mydata:
        newoutdict = {}
        baseval = np.median(rowlist[usetimes[0]])
        newoutdict[usetimes[0]] = [0.0 for tind in xrange(len(rowlist[usetimes[0]]))]
        for ttime in usetimes[1:]:
            newoutdict[ttime] = [math.log(titem/baseval,2.0) for titem in rowlist[ttime]]
        logmydata.append(deepcopy(newoutdict))
    for rowlist in logmydata:
        loggene2data.append([titem for ttime in usetimes for titem in rowlist[ttime]])
    return logmydata,loggene2data
           
           
def genData(fname,expfname,outfname):
    """generate data
    Args:
       fname:
       expfname:
       outfname:
    Returns:
       gene2data,timekeys:
    """
    gene2data = []
    timekeys = []
    allgenes = []
    count = 0
    emptyind = None
    with open(fname,"r") as infile:
       for line in infile:
           line = line.rstrip()
           splitted = line.split(",")
           if count == 0:
              timekeys = splitted[4:]
              emptyind = timekeys.index("")
              timekeys = timekeys[0:emptyind]
              count += 1
           else:
              assert splitted[0] not in allgenes
              genename = splitted[0]
              pref = genename.split("-")[0]
              if pref not in ["mghv","mcmv","mmu"]:
                 continue
              allgenes.append(genename)
              gene2data.append([float(titem) for titem in splitted[4:emptyind+4]])
    key2time = {keystr: float(keystr.replace("P","").replace("p","")) for keystr in timekeys}
    usetimes = sorted(set(key2time.values()))

    writeMiRNAData(outfname,gene2data,allgenes,timekeys)
               
    with open(expfname,"w") as outfile:
       outfile.write("Experiment Names\ttime inclusive\n")
       for keystr,tval in key2time.items():
           outfile.write("{0}\t{1}\n".format(keystr,tval))

    mydata = []
    for rowlist in gene2data:
        outdict = {}
        assert len(rowlist) == len(timekeys)
        for tind,titem in enumerate(rowlist):
            ttime = key2time[timekeys[tind]]
            outdict.setdefault(ttime,[])
            outdict[ttime].append(titem)
        mydata.append(deepcopy(outdict))  
    return gene2data,mydata,allgenes,timekeys,key2time,usetimes


fname = "LCM-mir.csv"
expfname = "expdes.txt"
outfname = "mirnadata.txt"
normoutfname = "norm_mirnadata.txt"
logoutfname = "logmirnadata.txt"
lognormoutfname = "lognorm_mirnadata.txt"
gene2data,data,allgenes,timekeys,key2time,usetimes = genData(fname,expfname,outfname)
#data time dict format, gene2data plain list
logdata,loggene2data = logTransform(data,gene2data,usetimes)
writeMiRNAData(logoutfname,loggene2data,allgenes,timekeys)
usetimes = sorted(data[0].keys())
for rowlist in data:
    assert len(set(rowlist.keys()) ^ set(usetimes)) == 0
normtype = "all" #"pertime"  
normdata = normData(data,usetimes,normtype)
normgene2data = [[] for tind in xrange(len(normdata))]
for rind,rowlist in enumerate(normdata):
    for ttime in usetimes:
        for titem in rowlist[ttime]:
            normgene2data[rind].append(titem)                      
writeMiRNAData(normoutfname,normgene2data,allgenes,timekeys)
exit(1)

normdata = normData(logdata,usetimes,normtype)
normgene2data = [[] for tind in xrange(len(normdata))]
for rind,rowlist in enumerate(normdata):
    for ttime in usetimes:
        for titem in rowlist[ttime]:
            normgene2data[rind].append(titem)                      
writeMiRNAData(lognormoutfname,normgene2data,allgenes,timekeys)


