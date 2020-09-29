#data utilities
import os
import sys
import math
import scipy as sp
import numpy as np
sys.path.append("../lib")
import Utilities
from copy import deepcopy
import random
import cPickle as cpickle


def readTraceFile(fpath):
    """read trace file
    Args:
       fpath:
    """
    with open(fpath,"rb") as infile:
        iterpointsall = cpickle.load(infile)
        y2knotslist = cpickle.load(infile)
        outsplineslist = cpickle.load(infile)
        avgerrorlist = cpickle.load(infile)
    return iterpointsall,y2knotslist,outsplineslist,avgerrorlist    


def saveExecutionTrace(savetracepath,iterpointsall,y2knotslist,outsplineslist,avgerrorlist):
    """saves execution trace
    """
    with open(savetracepath,"wb") as outfile:
        cpickle.dump(iterpointsall,outfile)
        cpickle.dump(y2knotslist,outfile)
        cpickle.dump(outsplineslist,outfile)
        cpickle.dump(avgerrorlist,outfile)

        
def getSubDataList(tdata,subtimes,mapval):
    """get data for subtimes in list
    """
    return [[tdata[yind][mapval[tpoint]] for tpoint in subtimes] for yind in xrange(len(tdata))]


def getSubDataDict(tdata,subtimes):
    """get data for subtimes in dict
    """ 
    return [{tpoint:tdata[yind][tpoint] for tpoint in subtimes} for yind in xrange(len(tdata))]


def readExpText(fname):
    """reads exp text
    """
    count = 0
    time2key,key2time = {}, {}
    with open(fname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if count != 0:
               keystr, time = line.split("\t")
               time = float(time)
               time2key.setdefault(time,set())
               time2key[time].add(keystr)
               key2time[keystr] = time
            else:
               count += 1
    return time2key,key2time           

def readMultiSampleData(fname,key2time):
    """reads multi sample data
    Args:
       fname:
       key2time:
    Returns:
       data:
    """
    gene2id, data = [], []
    ind2time = {}
    with open(fname,"r") as infile:
        for line in infile:
            uselines = line.split("\r")        
            for kind,keystr in enumerate(uselines[0].split("\t")[1:]):
                ind2time[kind+1] = key2time[keystr]

            for uind,useline in enumerate(uselines[1:]):
                splitted = uselines[uind+1].split("\t")
                gene = splitted[0]
                gene2id.append(gene)
                curvals = {}
                for ind,item in enumerate(splitted[1:]):
                    curtime = ind2time[ind+1]
                    curvals.setdefault(curtime,[])
                    curvals[curtime].append(float(item))
                data.append(deepcopy(curvals))
    return data,ind2time,gene2id


def readMultiSampleDataNewline(fname,key2time):
    """reads multi sample data
    Args:
       fname:
       key2time:
    Returns:
       data:
    """
    gene2id, data = [], []
    ind2time = {}
    count = 0
    with open(fname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            splitted = line.split("\t")
            if count == 0:
               for kind,keystr in enumerate(splitted[1:]):
                   ind2time[kind+1] = key2time[keystr]
               count +=1
               continue
            gene = splitted[0]
            gene2id.append(gene)
            curvals = {}
            for ind,item in enumerate(splitted[1:]):
                curtime = ind2time[ind+1]
                curvals.setdefault(curtime,[])
                curvals[curtime].append(float(item))
            data.append(deepcopy(curvals))
    return data,ind2time,gene2id


def readData(fname,expfname):
    """reads data
    Args:
       fname:
       expfname:
    Returns:
       shiftdata,shiftmeddata,usetimes,gene2ind,ind2time,time2key,key2time:   
    """
    time2key,key2time = readExpText(expfname)
    filestr = ""
    with open(fname,"r") as infile:
        for line in infile:
            filestr += line
    if filestr.find("\r") == -1:
       usedata,ind2time,gene2ind = readMultiSampleDataNewline(fname,key2time)
       usetimes = sorted(time2key.keys())
    else:
       data,ind2time,gene2ind = readMultiSampleData(fname,key2time)
       times = sorted(time2key.keys())
       usetimes = times[1:]
       usedata = deepcopy(data)
       for dpart in usedata:
           del dpart[times[0]]
                  
    shiftdata = normalizeShift(usedata,usetimes)
    shiftmeddata = getMedian(shiftdata)
    return shiftdata, shiftmeddata,usetimes,gene2ind,ind2time,time2key,key2time


def getMedian(shiftdata):
    """median dataset
    Args:
       shiftdata:
    Returns:
       shiftmeddata:   
    """
    shiftmeddata = []
    for cind,cdata in enumerate(shiftdata):
        shiftmeddata.append({ttime: np.median(cdata[ttime]) for ttime in cdata.keys()})
    return shiftmeddata

def getSampleData(usedata,count):
    """get a sample of data
    Args:
       usedata: 
       count: # of samples for each point
    Returns:
       sampledata:
    """
    flag = False
    for cind,cdata in enumerate(usedata):
        for ttime in cdata.keys():
            if len(cdata[ttime]) < count:
               flag = True
               break
    if flag:       
       print "not enough samples"
       return None
    
    sampledata = []
    for cind,cdata in enumerate(usedata):
        if type(cdata) == dict:
           sampledata.append({ttime:random.sample(cdata[ttime],count) for ttime in cdata.keys()})
        elif type(cdata) == list:
           sampledata.append([random.sample(cdata[ttime],count) for ttime in cdata.keys()])
    return sampledata


def getCurData(tdata,subtimes):
    """get data for subtimes
    Args:
       tdata:
       subtimes:
    """
    return [[tdata[yind][tpoint] for tpoint in subtimes] for yind in xrange(len(tdata))]


def getRemData(tdata,subtimes):
    """get data for rem times which is in dict form
    """ 
    return [{tpoint:tdata[yind][tpoint] for tpoint in subtimes} for yind in xrange(len(tdata))]


def getRemDataList(tdata,subtimes):
    """get data for rem times which is in dict form
    """ 
    return [{tpoint:[tdata[yind][tpoint]] for tpoint in subtimes} for yind in xrange(len(tdata))]


def normalizeShift(usedata,usetimes):
    """shifted normalization
    Args:
       usedata,usetimes:
    Returns:
       yvallist:   
    """
    ftime = sorted(usetimes)[0]
    yvallist = []
    for cind,cdata in enumerate(usedata):
        medval = np.median(cdata[ftime])
        cavgdata = {ttime:[titem-medval for titem in cdata[ttime]] for tind,ttime in enumerate(usetimes)}
        yvallist.append(deepcopy(cavgdata))
        assert abs(np.median(cavgdata[ftime])) <= 0.00000001
    return yvallist



