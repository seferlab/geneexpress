#Temporal selection algorithm
import sys
import math
import numpy as np
import random
import scipy as sp
import scipy.stats
import geneUtilities
import Utilities
import operator
import itertools
from copy import deepcopy
import scipy.interpolate
import cPickle as cpickle
import argparse
import os

def uniSelect(sorttimes,count):
    """init selection uniform
    Args:
       sorttimes,count:
    Returns:
       points: uniform points   
    """    
    #initialize
    curtimes = [sorttimes[0],sorttimes[-1]]
    time2pre = {time: tind+1 for tind,time in enumerate(sorttimes)}
    ratio = (len(sorttimes)-2)/(count-1)
    curind = 1
    unipoints = []
    for ttime in sorttimes[1:-1]:
        if len(unipoints) == count -2:
           break 
        if time2pre[ttime] >= curind*ratio+1:
           unipoints.append(ttime) 
           curind += 1
    return unipoints


def estQual(usespl,curremdata,rempoints,errtype="quad"):
    """estimate solution quality
    Args:
       usespl:
       curremdata:
       rempoints:
       errtype: error type
    Returns:
       errval: error value   
    """
    assert errtype in ["quad","abs"]
    errval = 0.0
    for rpoint in rempoints:
        estval = scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)
        if errtype == "quad":
           errval += sum([(estval-item)**2 for item in curremdata[rpoint]])
        elif errtype == "abs":
           errval += sum([abs(estval-item) for item in curremdata[rpoint]])
    return errval    


def initPoints(singledata,times,count,initmode):
    """initialize points
    Args:
       singledata: data with single value for each time point
       times: sorted times
       count:
       initmode:
    Returns:
       inittimes:   
    """
    #fixedpoints = [0.5, 2.5, 5.0, 10.0, 26.0]
    #fixedpoints = [0.5,7.0,28.0]
    fixedpoints = [times[0],times[-1]]
    assert initmode in ["change","equal"] and times == sorted(times)
    if initmode == "equal":
       points = uniSelect(times,count)
    elif initmode == "change":
       avgchange = {}
       for tind in xrange(1,len(times)-1):
           if times[tind] not in fixedpoints:
              ysumval = sum([abs(tyval[tind] - tyval[tind-1]) + abs(tyval[tind+1] - tyval[tind]) for tyval in singledata])
              avgchange[times[tind]] = ysumval
       points = [ttime[0] for ttime in sorted(avgchange.items(), key=lambda x: x[1], reverse=True)][0:count-2]
    inittimes = fixedpoints + points
    return sorted(inittimes)


def splineFitInterpolateMulti(curdata,remdata,curpoints,rempoints,weights):
    """spline fit for set of points among multiple reglambdas interpolating
    Args:
       curdata,remdata: 
       curpoints,rempoints:
       weights:
    Returns:
       outsplines,y2knots,sumval:
    """
    y2knots = []
    outsplines = []    
    sumval = 0.0
    for yind,curyvals in enumerate(curdata):
        minspl = scipy.interpolate.UnivariateSpline(curpoints, curyvals, s=0.0, k=3)
        for tpin,tpon in enumerate(curpoints):
            assert abs(scipy.interpolate.splev(tpon,minspl._eval_args, der=0, ext=0)-curyvals[tpin]) < 0.000001
        minval = weights[yind]*estQual(minspl,remdata[yind],rempoints,errtype="quad")
        sumval += minval
        outsplines.append(deepcopy(minspl))
        y2knots.append(list(minspl.get_knots()))
    rcount = sum([len(datam) for datam in remdata[0].values()])
    return outsplines,y2knots,sumval/(len(curdata)*rcount)


def runLOOCV(curyvals,curpoints,reglambdas,fixedpoints=[]):
    """leave one out cross validation
    Args:
       curyvals:
       curpoints:
       reglambdas:
       fixedpoints: do not remove them
    Returns:
       lambda2val:
    """
    degfit = min(3,len(curpoints)-2)
    lambda2val = {}
    for reglambda in reglambdas:
        tsumval = 0.0
        for lind,leavepoint in enumerate(curpoints):
            if leavepoint in fixedpoints:
               continue 
            curusepoints = list(curpoints)
            curusepoints.remove(leavepoint)
            useyvals = list(curyvals)
            del useyvals[lind]
            spl = scipy.interpolate.UnivariateSpline(curusepoints, useyvals, s=reglambda, k=degfit)
            if reglambda == 0.0:
               for tpin,tpon in enumerate(curusepoints):
                   assert abs(scipy.interpolate.splev(tpon,spl._eval_args, der=0, ext=0)-useyvals[tpin]) < 0.00001
            tsumval += (scipy.interpolate.splev(leavepoint,spl._eval_args, der=0, ext=0)-curyvals[lind])**2
        lambda2val[reglambda] = tsumval
    return lambda2val


def splineFitMulti(curdata,remdata,curdatamulti,curpoints,rempoints,weights,errormode,reglambdas=None):
    """spline fit for set of points among multiple reglambdas
    Args:
       curdata,remdata,curdatamulti:
       curpoints,rempoints:
       weights:
       errormode:
       reglambdas: lambdas to be used
    Returns:
       outsplines,y2knots,sumval:
    """
    y2knots = []
    outsplines = []
    bestlambdas = []    
    sumval = 0.0
    if reglambdas == None:
       #reglambdas = [0.001,0.01,0.1,0.25,0.5,1.0,2.0,5.0,10.0,15.0,20.0,25.0,50.0,75.0,100.0]
       reglambdas = [0.1,0.2,0.5,1.0,2.0,5.0,10.0]
    for yind,curyvals in enumerate(curdata):
        if errormode == "remain":
           fixedval = 0.0 
           for pind,ttime in enumerate(curpoints):
               trueval = curyvals[pind]
               fixedval += sum([(tval-trueval)**2 for tval in curdatamulti[yind][ttime]])
           fixedval *= weights[yind]
           
        if errormode == "remainloocv":
           sentfixedpoints = [min(curpoints), max(curpoints)]
           lambda2val = runLOOCV(curyvals,curpoints,reglambdas, sentfixedpoints)
           minlambda = sorted(lambda2val.items(), key=operator.itemgetter(1))[0][0]
           minspl = scipy.interpolate.UnivariateSpline(curpoints, curyvals, s=minlambda, k=3)
           minval = weights[yind]*estQual(minspl,remdata[yind],rempoints) 
        else:
           minval = 1000000000000000000000.0
           minspl,minlambda = None, None
           for reglambda in reglambdas:
               spl = scipy.interpolate.UnivariateSpline(curpoints, curyvals, s=reglambda, k=3)
               if reglambda == 0.0:
                  for tpin,tpon in enumerate(curpoints):
                      assert abs(scipy.interpolate.splev(tpon,spl._eval_args, der=0, ext=0)-curyvals[tpin]) < 0.00001
            
               tsumval = weights[yind]*estQual(spl,remdata[yind],rempoints)
               if errormode == "remain":
                  tsumval += fixedval
               elif errormode == "all":
                  tsumval += weights[yind]*estQual(spl,curdatamulti[yind],curpoints)
               if tsumval <= minval:
                  minval = tsumval
                  minspl = deepcopy(spl)
                  minlambda = reglambda
        sumval += minval
        outsplines.append(deepcopy(minspl))
        y2knots.append(list(minspl.get_knots()))
        bestlambdas.append(minlambda)
    rcount = sum([len(datam) for datam in remdata[0].values()])
    if errormode in ["all","remain"]:
       rcount += sum([len(datam) for datam in curdatamulti[0].values()])
    return outsplines,y2knots,sumval/(len(curdata)*rcount),bestlambdas
         

def splineFit(curdata,remdata,curpoints,rempoints,reglambda,weights):
    """spline fit for set of points
    Args:
       curdata,remdata: 
       curpoints,rempoints:
       reglambda:
       weights:
    Returns:
       outsplines,y2knots:
       sumval:
    """
    y2knots = []
    outsplines = []    
    sumval = 0.0
    for yind,curyvals in enumerate(curdata):
        spl = scipy.interpolate.UnivariateSpline(curpoints, curyvals, s=reglambda, k=3)
        if reglambda == 0.0:
           for tpin,tpon in enumerate(curpoints):
               assert abs(scipy.interpolate.splev(tpon,spl._eval_args, der=0, ext=0)-curyvals[tpin]) < 0.00001
        outsplines.append(deepcopy(spl))
        y2knots.append(list(spl.get_knots()))
        tsumval = weights[yind]*estQual(spl,remdata[yind],rempoints)
        sumval += tsumval
    rcount = sum([len(datam) for datam in remdata[0].values()])     
    return outsplines,y2knots,sumval/(len(curdata)*rcount) 


def runAnnealing():
    """runs simulated annealing
    Args:
    Returns:
    """
    return


def runExact(times,usedata,count,weights,initmode,multimode=1,errormode="remain",reglambdas=None,savetrace="greedytrace.pkl",fixedpoints=None):
    """runs exact method
    Args:
       times:
       usedata: which may not be shifted
       count:
       weights:
       initmode: initialization heuristic
       multimode: add/remove multiple points
       errormode:
       savetrace: save iters to file
       reglambdas:
       fixedpoints: 
    Returns:
       curtimes,y2knots:
       outsplines,avgerror:
    """
    assert errormode in ["remain","all","onlyremain","remainloocv"]
    if fixedpoints == None:
       fixedpoints = [times[0],times[-1]]
    mapval = {tval:tind for tind,tval in enumerate(times)}
    singledata,singledatadict = geneUtilities.makeSingleData(usedata,times)
    selecttimes = list(times)
    for tfixpoint in fixedpoints:
        selecttimes.remove(tfixpoint)
    if reglambdas == None:
       reglambdas = [0.1,0.2,0.5,1.0,2.0,2.5,5.0]
    assert count > len(fixedpoints) 
          
    outsplines, y2knots, avgerror, minsol = None,None,10000000000000000,None
    combs = list(itertools.combinations(selecttimes,count-len(fixedpoints)))
    random.shuffle(combs)
    allerrors = []
    for tempcurtimes in combs:
        print tempcurtimes
        curtimes = sorted(fixedpoints + list(tempcurtimes))
        rempoints = list(set(times) - set(curtimes))
        curdata = Utilities.getSubDataList(singledata,curtimes,mapval)
        remdata = Utilities.getSubDataDict(usedata,rempoints)
        curdatamulti = Utilities.getSubDataDict(usedata,curtimes)
        tempoutsplines,tempy2knots,tempavgerr,bestlambdas = splineFitMulti(curdata,remdata,curdatamulti,curtimes,rempoints,weights,errormode,reglambdas)
        allerrors.append(tempavgerr)
        if tempavgerr < avgerror:
           print curtimes,tempavgerr
           knotlens = {}
           for knots in tempy2knots:
               assert sorted(knots) == knots
               knotlens.setdefault(len(knots),0)
               knotlens[len(knots)] += 1
           print "knot distribution: ",knotlens
           
           lambdalens = {}
           for tlambda in bestlambdas:
               lambdalens.setdefault(tlambda,0)
               lambdalens[tlambda] += 1
           print "lambdadist ", lambdalens    
           
           avgerror = tempavgerr
           minsol = deepcopy(curtimes)
           outsplines = deepcopy(tempoutsplines)
           y2knots = deepcopy(tempy2knots)
           Utilities.saveExecutionTrace(savetrace.replace(".pkl","temp.pkl"),[minsol],[y2knots],[outsplines],allerrors)
    
    temptracepath = savetrace.replace(".pkl","temp.pkl")
    if os.path.exists(temptracepath):
       os.system("rm -rf {0}".format(temptracepath))        
    Utilities.saveExecutionTrace(savetrace,[minsol],[y2knots],[outsplines],allerrors)                      
    assert minsol == sorted(minsol)
    return minsol,y2knots,outsplines,avgerror


def runUniform(times,usedata,count,weights,multimode=1,errormode="remain",reglambdas=None,savetrace="greedytrace.pkl"):
    """selects uniform points without iteration
    Args:
       times: times to be used
       usedata: which may not be shifted
       count:
       weights:
       multimode: add/remove multiple points
       errormode:
       reglambdas:
       savetrace: save iters to file
    Returns:
       curtimes,y2knots:
       outsplines,avgerror:
    """
    assert errormode in ["remain","all"]
    fixedpoints = [times[0],times[-1]]
    mapval = {tval:tind for tind,tval in enumerate(times)} 
    singledata,singledatadict = geneUtilities.makeSingleData(usedata,times)
    curtimes = initPoints(singledata,times,count,"equal")
    rempoints = list(set(times) - set(curtimes))
    curdata = Utilities.getSubDataList(singledata,curtimes,mapval)
    remdata = Utilities.getSubDataDict(usedata,rempoints)
    curdatamulti = Utilities.getSubDataDict(usedata,curtimes)
    if reglambdas == None:
       reglambdas = [0.1,0.2,0.5,1.0,2.0,2.5,5.0]
    outsplines,y2knots,avgerror,bestlambdas = splineFitMulti(curdata,remdata,curdatamulti,curtimes,rempoints,weights,errormode,reglambdas)

    outsplineslist = [deepcopy(outsplines)]
    y2knotslist = [deepcopy(y2knots)]
    avgerrorlist = [avgerror]
    iterpointsall = [list(curtimes)]
    Utilities.saveExecutionTrace(savetrace,iterpointsall,y2knotslist,outsplineslist,avgerrorlist)
    assert curtimes == sorted(curtimes)     
    return curtimes,y2knots,outsplines,avgerror
           

def runGreedy(times,usedata,count,weights,initmode,multimode=1,errormode="remain",reglambdas=None,savetrace="greedytrace.pkl"):
    """runs greedy method
    Args:
       times:
       rawusedata: which may not be shifted
       count:
       weights:
       initmode: initialization heuristic
       multimode: add/remove multiple points
       errormode:
       reglambdas:
       savetrace: save iters to file
    Returns:
       curtimes,y2knots:
       outsplines,avgerror:
    """
    assert errormode in ["remain","all"]
    fixedpoints = [times[0],times[-1]]
    mapval = {tval:tind for tind,tval in enumerate(times)} 
    singledata,singledatadict = geneUtilities.makeSingleData(usedata,times)
    curtimes = initPoints(singledata,times,count,initmode)
    #if count == 5: #proteomics
    #   curtimes = [0.5, 4.0, 10.0, 23.0, 28.0]
    #elif count == 8: #proteomics
    #   curtimes = [0.5, 1.0, 4.0, 7.0, 10.0, 19.0, 23.0, 28.0]
       
    #   curtimes = [0.5,4.0,7.0,14.0,28.0]  
    #   curtimes = [0.5, 6.0, 9.5, 19.0, 28.0]
    #   curtimes = [0.5, 6.5, 9.0, 19.0, 28.0]
    #elif count == 8:
    #   curtimes = [0.5, 3.5, 7.0, 9.5, 13.5, 17.0, 23.0, 28.0]   
    
    rempoints = list(set(times) - set(curtimes))
    curdata = Utilities.getSubDataList(singledata,curtimes,mapval)
    remdata = Utilities.getSubDataDict(usedata,rempoints)
    curdatamulti = Utilities.getSubDataDict(usedata,curtimes)
    if reglambdas == None:
       reglambdas = [0.1,0.2,0.5,1.0,2.0,2.5,5.0]
    outsplines,y2knots,avgerror,bestlambdas = splineFitMulti(curdata,remdata,curdatamulti,curtimes,rempoints,weights,errormode,reglambdas) #initial fitting

    outsplineslist = [deepcopy(outsplines)]
    y2knotslist = [deepcopy(y2knots)]
    avgerrorlist = [avgerror]
    iterpointsall = [list(curtimes)]
    itcount = 0
    
    while True: #loop
        print "iter: {0}, error:{1} ".format(itcount,avgerror)
        print curtimes
        minsol = None
        for addpoints in itertools.combinations(rempoints,multimode):
            for delpoints in itertools.combinations(curtimes,multimode):
                if len(set(delpoints) & set(fixedpoints)) != 0:
                   continue 
                newrempoints, newpoints = list(rempoints), list(curtimes)
                [newrempoints.remove(addpoint) for addpoint in addpoints]
                newrempoints.extend(delpoints)
                [newpoints.remove(delpoint) for delpoint in delpoints]
                newpoints.extend(addpoints)
                newpoints = sorted(newpoints)

                newcurdata = Utilities.getSubDataList(singledata,newpoints,mapval)
                newremdata = Utilities.getSubDataDict(usedata,newrempoints)
                newcurdatamulti = Utilities.getSubDataDict(usedata,newpoints)
                tempoutsplines,tempy2knots,tempavgerr,bestlambdas = splineFitMulti(newcurdata,newremdata,newcurdatamulti,newpoints,newrempoints,weights,errormode,reglambdas)
                if tempavgerr < avgerror:
                   avgerror = tempavgerr
                   minsol = (addpoints,delpoints)
                   outsplines = list(tempoutsplines)
                   y2knots = list(tempy2knots)
            
        if minsol == None:
           break
       
        itcount += 1
        [rempoints.remove(tpoint) for tpoint in minsol[0]]
        rempoints.extend(minsol[1])
        [curtimes.remove(tpoint) for tpoint in minsol[1]]
        curtimes.extend(minsol[0])
        curtimes = sorted(curtimes)
        outsplineslist.append(deepcopy(outsplines))
        y2knotslist.append(deepcopy(y2knots))
        avgerrorlist.append(avgerror)
        iterpointsall.append(list(curtimes))
        Utilities.saveExecutionTrace(savetrace,iterpointsall,y2knotslist,outsplineslist,avgerrorlist)
        
    assert curtimes == sorted(curtimes)     
    return curtimes,y2knots,outsplines,avgerror


def getParser():
    """get command line parser
    """
    parser = argparse.ArgumentParser(description='Analyze clussters')
    parser.add_argument("-o", "--clustoutdir", type=str, default="clusters", help="output directory")
    parser.add_argument("-f", "--respath", type=str, default="output.txt", help="input file in xlsx format")
    #parser.add_argument("-f", "--respath", type=int, default=8, help="number of clusters. Default is run for 8 clusters")
    parser.add_argument("-m", "--multimode", type=int, default=1, help="number of points to change in each iterations. Default is 1.")
    return parser


def main(respath,multimode):
    """
    """
    curtimes,y2knots,outsplines,avgerror = runGreedy(times,rawusedata,count,weights,initmode,multimode)
    with open(respath,"w") as outfile:
        outfile.write("selected times: {0}\n".format("\t".join([str(ttime) for ttime in curtimes])))
        outfile.write("avg error: {0}\n".format(avgerror))
        for gind,gene in enumerate(allgenes):
            outfile.write("{0}\t{1}\n".format(gene,"\t".join([str(ttime) for ttime in curtimes])))
            

if False:            
 respath = "output.txt"
 multimode = 1
 parser = getParser()
 args = parser.parse_args()
 clustoutdir = args.clustoutdir
 if not os.path.exists(clustoutdir):
   os.makedirs(clustoutdir) 
 fname = args.fname
 clustnum = args.clustnum
 data,ind2userid,fields = readFile(fname)
 moddata,cat2assign = transData(data)

 main(respath,multimode)         
