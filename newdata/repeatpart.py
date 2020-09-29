#repeat analysis
import os
import sys
import math
import random
import operator
import scipy as sp
import numpy as np
sys.path.append("../lib")
import Utilities
from copy import deepcopy
import Tempselect


def assignTime(allpoints,usetimes,existset):
    """time assignemnt 
    """
    point2point = {}  
    for allpoint in allpoints:
        if allpoint in existset:
           continue 
        dists = {ttime:abs(ttime-allpoint) for ttime in usetimes}
        sorted_dists = sorted(dists.items(), key=operator.itemgetter(1))
        selectlist = [sorted_dists[0][0],sorted_dists[1][0]]
        closetime = random.choice(selectlist)
        point2point[allpoint] = closetime
    return point2point


def genExtendedData(allpoints,usetimes,existset,shiftdata,shiftmeddata):
    """
    """           
    point2point = assignTime(allpoints,usetimes,existset)
    addshiftdata = []
    addshiftmeddata = []
    for tind,tdata in enumerate(shiftdata):
        outdict = deepcopy(tdata)
        outdict2 = deepcopy(shiftmeddata[tind])
        for point in allpoints:
            if outdict.has_key(point):
               assert outdict2.has_key(point)
               continue  
            outdict[point] = list(outdict[point2point[point]])                 
            outdict2[point] = outdict2[point2point[point]]  
        addshiftdata.append(outdict)
        addshiftmeddata.append(outdict2)
    shiftdata = deepcopy(addshiftdata)
    shiftmeddata = deepcopy(addshiftmeddata)
    return shiftdata,shiftmeddata

def getSharedTimes(allpoints,usetimes):
    """
    """    
    existset = []
    for allpoint in allpoints:
        flag = False
        for tttime in usetimes:
            if abs(allpoint) <= 0.00001:
               flag = True
               break
        if flag:
           existset.append(allpoint)
    return existset


def runAnalysisExtraPoints(budget,shiftdata,shiftmeddata,usetimes,maptype,repcounts=[1,2]):
    """run analysis
    Args:
       budget:
       shiftdata,shiftmeddata:
       usetimes:
       maptype:
       repcounts: repeat counts
    Returns:
       errdict: 
    """
    assert maptype in ["det","stoc"]
    assert budget % repcounts[0] == 0 and budget % repcounts[1] == 0 and len(repcounts) >= 2
    start, end = min(usetimes),max(usetimes)
    samplepoints = {}
    for repcount in repcounts:
        usebudget = int(math.floor(float(budget)/repcount))
        div = float(end-start)/(usebudget+1)
        samplepoints[repcount] = [start+div*tind for tind in xrange(1,usebudget+1)]
        #div = float(end-start)/(usebudget-1)
        #samplepoints[repcount] = [start+div*tind for tind in xrange(0,usebudget)]
    weights = [1.0]*len(shiftdata)
    
    allpoints = set(ttime for repstr in samplepoints.keys() for ttime in samplepoints[repstr])
    existset = getSharedTimes(allpoints,usetimes)
    if maptype == "det":
       shiftdata,shiftmeddata = genExtendedData(allpoints,usetimes,existset,shiftdata,shiftmeddata)
       testdata = Utilities.getRemDataList(shiftmeddata,usetimes)   
    errdict = {repcount:[] for repcount in repcounts}
    samplecount = 1000
    for scount in xrange(samplecount):
        if maptype == "stoc":
           shiftdata,shiftmeddata = genExtendedData(allpoints,usetimes,existset,shiftdata,shiftmeddata)
           testdata = Utilities.getRemDataList(shiftmeddata,usetimes)
        for repcount in repcounts:
            singledata = Utilities.getSampleData(shiftdata,repcount)
            singlemeddata = Utilities.getMedian(singledata)
            curdata = Utilities.getCurData(singlemeddata,samplepoints[repcount])
            outsplines,y2knots,curavgsumval = Tempselect.splineFitInterpolateMulti(curdata,testdata,samplepoints[repcount],usetimes,weights)
            errdict[repcount].append(curavgsumval)
        if scount % 100 == 0:
           for keystr in errdict.keys():
               print keystr,np.mean(errdict[keystr])     
    return errdict        

def readData():
    """
    """
    expfname = "expdes.txt"     
    time2key,key2time = Utilities.readExpText(expfname)
    fname = "127 LCM time course Quantile Normalized logbased 2 transformed.txt"
    data,ind2time,gene2ind = Utilities.readMultiSampleData(fname,key2time)
    #modify data
    times = sorted(time2key.keys())
    usetimes = times[1:]
    usedata = deepcopy(data)
    for dpart in usedata:
        del dpart[times[0]]

    shiftdata = Utilities.normalizeShift(usedata,usetimes)
    shiftmeddata = Utilities.getMedian(shiftdata)
    return shiftdata, shiftmeddata,usetimes

def genSynData(shiftdata,shiftmeddata,noise):
    """
    """
    time2count = {ttime:len(shiftdata[0][ttime]) for ttime in shiftdata[0].keys()}
    newshiftdata = []
    for tdata in shiftmeddata:
        outdict = {}
        for ttime in tdata.keys():
            vals = [tdata[ttime]]
            for tind in xrange(time2count[ttime]-1):
                vals.append(tdata[ttime]+random.gauss(0,noise))
            outdict[ttime] = list(vals)
        newshiftdata.append(deepcopy(outdict))
    shiftdata = deepcopy(newshiftdata)
    shiftmeddata = Utilities.getMedian(shiftdata)        
    return shiftdata,shiftmeddata


budget = sys.argv[1]
runmode = sys.argv[2]
outfolder = sys.argv[3]
simulmode = sys.argv[4]
if len(sys.argv) == 6:
   noise = sys.argv[5]
assert runmode in ["real","syn"] and simulmode in ["det","stoc"]
budget = int(budget)
noise = float(noise)

if runmode == "real":
   fpath = "{0}/{1}.txt".format(outfolder,budget)
elif runmode == "syn":
   fpath = "{0}/{1}_{2}.txt".format(outfolder,budget,noise)
if os.path.exists(fpath):
   exit(1)
     
shiftdata, shiftmeddata,usetimes = readData()
#for tind,tdata in enumerate(shiftdata):
#    for ttime in shiftdata[tind].keys():
#        shiftdata[tind][ttime] = [titem*0.5 for titem in shiftdata[tind][ttime]]
#        shiftmeddata[tind][ttime] *= 0.5     
if runmode == "syn":
   shiftdata,shiftmeddata = genSynData(shiftdata,shiftmeddata,noise)
elif runmode == "real":
   for sind,shiftdataitem in enumerate(shiftdata):
       shiftdata[sind][7.5].append(shiftdataitem[7.5][1])     
repcounts = [1,2,3]

errdict = runAnalysisExtraPoints(budget,shiftdata,shiftmeddata,usetimes,"det",repcounts)
for keystr in errdict.keys():
    print keystr,np.mean(errdict[keystr])
with open(fpath,"w") as outfile:
    for keystr in errdict.keys():
        outfile.write("{0}\t{1}\n".format(keystr,np.mean(errdict[keystr])))    
exit(1)


if False:
 for sind,shiftdataitem in enumerate(shiftdata):
    shiftdataitem[7.5].append(shiftdataitem[7.5][1]) #??
    shiftmeddataitem = shiftmeddata[sind] 
    errdict = runAnalysisExtraPoints(12,[shiftdataitem],[shiftmeddataitem],usetimes,repcounts)
    denerr = np.mean(errdict[repcounts[0]])
    reperr = np.mean(errdict[repcounts[1]])
    #reperr2 = np.mean(errdict[repcounts[2]])
    if reperr >= denerr:
        print sind,denerr,reperr,gene2ind[sind]
 exit(1)
    
errdict = runAnalysisExtraPoints(12,shiftdata,shiftmeddata,usetimes,simulmode,repcounts)
#print errdict
for keystr in errdict.keys():
    print keystr,np.mean(errdict[keystr])
exit(1)

plotpath = "errorplot.png"
alldenerrs = []
allreperrs = []
for budget in [2*titem for titem in xrange(4,10)]:
    denerr,reperr = runAnalysis(budget,shiftdata,shiftmeddata,usetimes)
    print budget,denerr,reperr
    alldenerrs.append(denerr)
    allreperrs.append(reperr)
print alldenerrs
print allreperrs    




    
