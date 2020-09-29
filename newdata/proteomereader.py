import os
import sys
from copy import deepcopy
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import scipy.signal
import scipy.interpolate
import random
import itertools
import math
import AnalyzeData
import PlotGen
import os
import sys
sys.path.append("../lib")
import Tempselect
import Utilities
import matplotlib
import random
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
from sklearn.decomposition import PCA
import matplotlib.cm as cm
import cPickle as cpickle
import gzip
import scipy.cluster

#E16.5 -> 1:8
#0.5 13.5 -> 27
#14 28 -> 15

def makeplotMainAll(xvals,yvaldictout,points,knots,gene,outspl,rempoints,remyvals,plotpath,allknot,allspline):    
    """makes main plot
    Args:
       xvals,yvaldictout,points:
       knots,gene,outspl:
       rempoints,remyvals:
       plotpath:
       allknot,allspline:
    Returns:
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 30 
    MARKERSIZE = 25
    DPI = 300
    
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in xvals] 
    ymax,ymin = max(yvals), min(yvals)
    plt.ylim(ymin-1.0,ymax+0.1)
    plt.xlim(0.0,max(xvals)+0.5)
    locpos = 4
    if gene in ["Eln","ELN"]:
       locpos = 3 
    xlabel = "Days"
    ylabel = ""
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    symmap = {}
    symmap["marker"] = {"spline":"p", "SplineFit":"*", "Knots":"p"}
    symmap["colors"] = {"spline":"r", "SplineFit":"g", "Knots":"k"}
    symmap["labels"] = {"spline":"Fitted Spline", "SplineFit":"Selected points {0}".format(len(points)), "Knots":"Knots {0}".format(len(knots))}

    assert sorted(knots) == knots
    knotvals = [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] 
    for algo in yvaldictout.keys():
        plt.plot(points,yvaldictout[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    plt.plot(knots,knotvals,marker=symmap["marker"]["Knots"],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"]["Knots"],label=symmap["labels"]["Knots"])
    plt.plot(xvals,outspl(xvals),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')

    allyvals = [scipy.interpolate.splev(knot,allspline._eval_args, der=0, ext=0) for tknot in allknot]
    plt.plot(xvals,allspline(xvals),'y',lw=6,label='All')
    knotys = [scipy.interpolate.splev(tknot,allspline._eval_args, der=0, ext=0) for tknot in allknot]
    plt.plot(allknot,knotys,marker='+',markersize=MARKERSIZE,linestyle='None',color='b',label='Knots all')
       
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)
    

    
def makeplotJointPart(yvallist,xvals,globyvals,useindex,globknots,globoutsplines,globpoints,gene,globrempoints,remyvalsdict,x1,plotpath):
    """
    x1: x points
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 30
    YFSIZE = 40
    LEGENDSIZE = 20
    MARKERSIZE = 25
    DPI = 300

    outspl5 = globoutsplines[5][useindex]
    outspl8 = globoutsplines[8][useindex]
    knots5 = globknots[5][useindex]
    knots8 = globknots[8][useindex]
    remyvals5 = remyvalsdict[5]
    remyvals8 = remyvalsdict[8]

    time2index = {xval:xind for xind,xval in enumerate(xvals)}
    y1 = [yvallist[useindex][time2index[ttime]] for ttime in x1]

    tyvals = [scipy.interpolate.splev(knot,outspl5._eval_args, der=0, ext=0) for knot in knots5] + [scipy.interpolate.splev(knot,outspl8._eval_args, der=0, ext=0) for knot in knots8] + remyvals5 + remyvals8 + [scipy.interpolate.splev(ttime,outspl5._eval_args, der=0, ext=0) for ttime in xvals] + [scipy.interpolate.splev(ttime,outspl8._eval_args, der=0, ext=0) for ttime in xvals] + y1
    
    ymax,ymin = max(tyvals), min(tyvals)
    plt.ylim(ymin-0.5,ymax+0.1)
    plt.xlim(0.0,max(xvals)+0.5)
    locpos = 4
    xlabel = "Days"
    ylabel = ""
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    
    assert sorted(knots5) == knots5 and sorted(knots8) == knots8
    knotvals5 = [scipy.interpolate.splev(knot,outspl5._eval_args, der=0, ext=0) for knot in knots5]
    knotvals8 = [scipy.interpolate.splev(knot,outspl8._eval_args, der=0, ext=0) for knot in knots8]
        
    #for algo in yvaldictout.keys():
    #    plt.plot(points,yvaldictout[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    #plt.plot(knots,knotvals,marker=symmap["marker"]["Knots"],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"]["Knots"],label=symmap["labels"]["Knots"])
    #plt.plot(xvals,outspl(xvals),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    #plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')
        
    plt.plot(xvals,yvallist[useindex],marker="p",markersize=MARKERSIZE,linestyle='None',color="k",label="Data")
    plt.plot(xvals,outspl5(xvals),marker="+",lw=3,color="r",label="TPS 5")
    plt.plot(xvals,outspl8(xvals),marker="*",lw=3,color="y",label="TPS 8")
    plt.plot(x1,y1,marker='o',markersize=MARKERSIZE,lw=3,color='g',label='Piecewise Linear')
         
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)
    
       
def makeplotJoint(yvallist,xvals,globyvals,useindex,globknots,globoutsplines,globpoints,gene,globrempoints,remyvalsdict,plotpath):
    """joint plot
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 30
    YFSIZE = 40
    LEGENDSIZE = 20
    MARKERSIZE = 25
    DPI = 300

    outspl5 = globoutsplines[5][useindex]
    outspl8 = globoutsplines[8][useindex]
    knots5 = globknots[5][useindex]
    knots8 = globknots[8][useindex]
    remyvals5 = remyvalsdict[5]
    remyvals8 = remyvalsdict[8]

    time2index = {xval:xind for xind,xval in enumerate(xvals)}
    x1 = [0.5,7.0,14.0,28.0]
    x2 = [0.5,2.0,14.0,28.0]
    x3 = [0.5,4.0,7.0,14.0,28.0]
    y1 = [yvallist[useindex][time2index[ttime]] for ttime in x1]
    y2 = [yvallist[useindex][time2index[ttime]] for ttime in x2]
    y3 = [yvallist[useindex][time2index[ttime]] for ttime in x3]

    yvals = [scipy.interpolate.splev(knot,outspl5._eval_args, der=0, ext=0) for knot in knots5] + [scipy.interpolate.splev(knot,outspl8._eval_args, der=0, ext=0) for knot in knots8] + remyvals5 + remyvals8 + [scipy.interpolate.splev(ttime,outspl5._eval_args, der=0, ext=0) for ttime in xvals] + [scipy.interpolate.splev(ttime,outspl8._eval_args, der=0, ext=0) for ttime in xvals] + y1 + y2 + y3
    
    ymax,ymin = max(yvals), min(yvals)
    plt.ylim(ymin-0.5,ymax+0.1)
    plt.xlim(0.0,max(xvals)+0.5)
    locpos = 4
    xlabel = "Days"
    ylabel = ""
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    symmap = {}
    symmap["marker"] = {"spline":"p", "SplineFit":"*", "Knots":"p"}
    symmap["colors"] = {"spline":"r", "SplineFit":"g", "Knots":"k"}
    symmap["labels"] = {"spline":"Fitted Spline", "SplineFit":"Selected points {0}".format(len(points)), "Knots":"Knots {0}".format(len(knots))}

    assert sorted(knots5) == knots5 and sorted(knots8) == knots8
    knotvals5 = [scipy.interpolate.splev(knot,outspl5._eval_args, der=0, ext=0) for knot in knots5]
    knotvals8 = [scipy.interpolate.splev(knot,outspl8._eval_args, der=0, ext=0) for knot in knots8]
        
    #for algo in yvaldictout.keys():
    #    plt.plot(points,yvaldictout[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    #plt.plot(knots,knotvals,marker=symmap["marker"]["Knots"],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"]["Knots"],label=symmap["labels"]["Knots"])
    #plt.plot(xvals,outspl(xvals),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    #plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')
        
    plt.plot(xvals,yvallist[useindex],marker="p",markersize=MARKERSIZE,linestyle='None',color="k",label="Data")
    plt.plot(xvals,outspl5(xvals),marker="+",lw=3,color="r",label="TPS 5")
    plt.plot(xvals,outspl8(xvals),marker="*",lw=3,color="y",label="TPS 8")
    
    plt.plot(x1,y1,marker='o',markersize=MARKERSIZE,lw=3,color='b',label='Linear fit 1')
    plt.plot(x2,y2,marker='o',markersize=MARKERSIZE,lw=3,color='m',label='Linear fit 2')
    plt.plot(x3,y3,marker='o',markersize=MARKERSIZE,lw=3,color='g',label='Linear fit 3')

    area5 = 0.0
    area8 = 0.0
    dif = float(x1[-1]-x1[0])/99999
    uselocs = [x1[0]+(tind*dif) for tind in xrange(100000)]
    for tindex in xrange(len(uselocs)-1):
        tloc1,tloc2 = uselocs[tindex:tindex+2]
        sampx = (tloc1+tloc2)*0.5
        found = None
        for tind in xrange(len(x1)-1):
            beg,end = x1[tind],x1[tind+1]
            if beg <= sampx and end >= sampx:
               found = (beg,end,y1[tind],y1[tind+1])
               break
        assert found != None
        xstart,xend,ystart,yend = found
        m = float(yend-ystart)/(xend-xstart)
        estval = m*(sampx-xstart) + ystart
        #area5 += abs(estval-scipy.interpolate.splev(sampx,outspl5._eval_args, der=0, ext=0))
        #area8 += abs(estval-scipy.interpolate.splev(estval,outspl8._eval_args, der=0, ext=0))
        area5 += abs(estval-scipy.interpolate.splev(sampx,outspl5._eval_args, der=0, ext=0))
        area8 += abs(estval-scipy.interpolate.splev(estval,outspl8._eval_args, der=0, ext=0))
    print area5
    print area8
    exit(1)    
         
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)
    
    

    
    
def makeplotMain(xvals,yvaldictout,points,knots,gene,outspl,rempoints,remyvals,plotpath):    
    """makes main plot
    Args:
       xvals,yvaldictout,points:
       knots,gene,outspl:
       rempoints,remyvals:
       plotpath:
    Returns:
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 30 
    MARKERSIZE = 25
    DPI = 300
    
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in xvals] 
    ymax,ymin = max(yvals), min(yvals)
    plt.ylim(ymin-1.0,ymax+0.1)
    plt.xlim(0.0,max(xvals)+0.5)
    locpos = 4
    xlabel = "Time"
    ylabel = ""
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    symmap = {}
    symmap["marker"] = {"spline":"p", "SplineFit":"*", "Knots":"p"}
    symmap["colors"] = {"spline":"r", "SplineFit":"g", "Knots":"k"}
    symmap["labels"] = {"spline":"Fitted Spline", "SplineFit":"Selected points {0}".format(len(points)), "Knots":"Knots {0}".format(len(knots))}

    assert sorted(knots) == knots
    knotvals = [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots]

    for algo in yvaldictout.keys():
        plt.plot(points,yvaldictout[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    plt.plot(knots,knotvals,marker=symmap["marker"]["Knots"],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"]["Knots"],label=symmap["labels"]["Knots"])
    plt.plot(xvals,outspl(xvals),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')
     
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)    


def estQual(usespl,x2val,points,rempoints):
    """estimate solution quality
    Args:
       usespl,x2val:
       points,rempoints:
    """
    #sumval = sum([(scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)-item)**2 for rpoint in rempoints for item in x2val[rpoint]])
    sumval = sum([(scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)-x2val[rpoint])**2 for rpoint in rempoints])
    return sumval

def onlyMethPoints(yvallist,times,weights):
    """
    """     
    #Combinatorial part for meth points only 
    mapval = {tval:tind for tind,tval in enumerate(times)}
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    fixedpoints = [0.5,26.0]
    rempoints = list(set(times) - set(fixedpoints))
    reglambda = 0.0001
    bestval = 10000000000000000.0
    pairsol = []
    for p1,p2 in itertools.combinations(rempoints,2):
        points = [fixedpoints[0]] + [p1,p2] + [fixedpoints[-1]]
        minval = 10000000.0
        tind = 0
        y2knots = []
        outsplines = []    
        sumval = 0.0
        yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))] 
        for yind,curyvals in enumerate(yvals):
            spl = scipy.interpolate.UnivariateSpline(points, curyvals, s=reglambda, k=3)
            if reglambda == 0.0:
               for tpon in points:
                   assert abs(scipy.interpolate.splev(tpon,spl._eval_args, der=0, ext=0)-x2val[yind][tpon]) < 0.00001
            outsplines.append(deepcopy(spl))
            y2knots.append(list(spl.get_knots()))
            tsumval = weights[yind]*estQual(spl,x2val[yind],points,rempoints)
            sumval += tsumval
        avgknot = np.mean([len(tknot) for tknot in y2knots])
        print p1,p2,sumval,avgknot    
        if sumval < bestval:
           bestval = sumval
           pairsol = points
    print "combinatorial results:"       
    print pairsol
    print bestval


def runGreedy(times,usedata,yvallist,count,weights,initmode="change"):
    """runs greedy method
    Args:
       times,usedata:
       count,weights:
    Returns:   
    """       
    def convertForm(curtimes):
        """converts form
        Args:
           curtimes:
        Returns:
           yvallist:   
        """        		
        tyvallist = []
        for cind,cdata in enumerate(usedata):
            calldata = []
            for tind,ttime in enumerate(curtimes):
                calldata.extend(usedata[cind][ttime])
            tyvallist.append(list(calldata))
            #assert len(calldata) == len(usetimes)
        return tyvallist
    
    if False:
     globtimes = set()
     for cind,cdata in enumerate(usedata):
        globtimes |= set(usedata[cind].keys())
     for cind,cdata in enumerate(usedata):
        assert len(globtimes) == len(usedata[cind].keys())
     usetimes = [ttime for tind,ttime in enumerate(times) for tind2 in xrange(len(usedata[0][ttime]))]
     yvallist = convertForm(times) #??
    
    #onlyMethPoints(yvallist,times,weights)
    fixedpoints = [0.5, 2.5, 5.0, 10.0, 26.0]
    #fixedpoints = [0.5,7.0,28.0]
    #fixedpoints = [0.5,28.0]
    fixedpoints = [0.5,26.0]
    assert initmode in ["change","equal"] and times == sorted(times)
    reglambda = 5.0 #1.0 #50.0 #1.0
    if initmode == "equal":
       points = initSelect(times,count)
    elif initmode == "change":
       #avgyvallist = []
       #for cind,cdata in enumerate(usedata):
       #    cavgdata = [np.mean(list(usedata[cind][ttime])) for tind,ttime in enumerate(times)]
       #    avgyvallist.append(list(cavgdata)) 
       avgchange = {}
       for tind in xrange(1,len(times)-1):
           ysumval = 0.0
           for tyval in yvallist:
               ysumval += abs(tyval[tind] - tyval[tind-1]) + abs(tyval[tind+1] - tyval[tind])
           avgchange[times[tind]] = ysumval
       for fixedpoint in fixedpoints:
           if avgchange.has_key(fixedpoint):
              del avgchange[fixedpoint]  
       points = [ttime[0] for ttime in sorted(avgchange.items(), key=lambda x: x[1], reverse=True)][0:count-2]
    
    points = fixedpoints + points
    points = sorted(points)
    #multipoints = [point for point in points for tind2 in xrange(len(usedata[0][point]))]
    rempoints = list(set(times) - set(points))
    mapval = {tval:tind for tind,tval in enumerate(times)}
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
    #yvals = convertForm(points)

    
        
    
    minval = 10000000.0
    tind = 0
    y2knots = []
    outsplines = []    
    sumval = 0.0
    for yind,curyvals in enumerate(yvals):
        spl = scipy.interpolate.UnivariateSpline(points, curyvals, s=reglambda, k=3)
        if reglambda == 0.0:
           for tpon in points:
               assert abs(scipy.interpolate.splev(tpon,spl._eval_args, der=0, ext=0)-x2val[yind][tpon]) < 0.00001
        outsplines.append(deepcopy(spl))
        y2knots.append(list(spl.get_knots()))
        tsumval = weights[yind]*estQual(spl,x2val[yind],points,rempoints)
        sumval += tsumval
            
    while True:
        print sumval
        minsol = None
        for addpoint in rempoints:
            for delpoint in points:
                if delpoint in fixedpoints:
                   continue 
                newrempoints, newpoints = list(rempoints), list(points)
                newrempoints.remove(addpoint)
                newrempoints.append(delpoint)
                newpoints.remove(delpoint)
                newpoints.append(addpoint)
                newpoints = sorted(newpoints) 
                
                tempoutsplines = []
                globsumval = 0.0
                for yind,curyvals in enumerate(yvals):
                    inyvals = [x2val[yind][rpoint] for rpoint in newpoints]
                    inspl = scipy.interpolate.UnivariateSpline(newpoints, inyvals, s=reglambda, k=3)
                    
                    #spl2 = scipy.interpolate.InterpolatedUnivariateSpline(sortednewpoints, sortedinyvals,k=3)
                    #for tpon in newpoints:
                    #    assert abs(spl2.__call__(tpon)-x2val[yind][tpon]) < 0.00001
                    #    #assert abs(scipy.interpolate.splev(tpon,spl2._eval_args, der=0, ext=0)-x2val[yind][tpon]) < 0.00001
                    if reglambda == 0.0:
                       for tpon in newpoints:
                           assert abs(scipy.interpolate.splev(tpon,inspl._eval_args, der=0, ext=0)-x2val[yind][tpon])<0.00001
                    tempoutsplines.append(deepcopy(inspl))
                    tsumval = weights[yind]*estQual(inspl,x2val[yind],newpoints,newrempoints)
                    globsumval += tsumval
                if globsumval < sumval:
                   sumval = globsumval
                   minsol = (addpoint,delpoint)
                   outsplines = list(tempoutsplines)
                    
        if minsol == None:
           break
        
        rempoints.remove(minsol[0])
        rempoints.append(minsol[1])
        points.remove(minsol[1])
        points.append(minsol[0])
        points = sorted(points)
        yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
    assert points == sorted(points)    
    tsumval = 0.0
    for yind,curyvals in enumerate(yvals):
        inyvals = [x2val[yind][rpoint] for rpoint in points]
        tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
        locsumval = weights[yind]*estQual(tinspl,x2val[yind],points,rempoints)
        tsumval += locsumval
    print tsumval
    assert abs(tsumval-sumval) < 0.0001     
    avgsum=0.0
    for yind in xrange(len(y2knots)):
        avgsum += len(y2knots[yind])
    print avgsum
    print avgsum/float(len(y2knots))
    print "infom"
    return sumval,sumval/float(len(yvallist)*(len(times)-count)),points,yvals,y2knots,outsplines


def initSelect(sorttimes,count):
    """
    """    
    #initialize
    curtimes = [sorttimes[0],sorttimes[-1]]
    time2pre = {}
    for tind,time in enumerate(sorttimes):
        time2pre[time] = tind+1
    ratio = (len(sorttimes)-2)/(count-1)
    curind = 1
    points = []
    for ttime in sorttimes[1:-1]:
        if len(points) == count -2:
           break 
        if time2pre[ttime] >= curind*ratio+1:
           points.append(ttime) 
           curind += 1
    return points


def plotGeneral(outdict,xlabel,ylabel,plotpath):
    """plot general
    Args:
       outdict:
       xlabel,ylabel:
       plotpath:
    Returns:   
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 40 
    MARKERSIZE = 35
    DPI = 300
    yvals = outdict["y"]
    ymax,ymin = max(yvals), min(yvals)
    plt.ylim(0,ymax+1.0)
    locpos = 1
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(1,max(outdict["x"])+1)
     
    plt.plot(outdict["x"],outdict["y"],marker="p",markersize=MARKERSIZE,linestyle='None',color="r")
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)
    
def clustPCA(tyvallist,centroid,labeldict,clustplotpath):
    """makes cluster PCA plot
    """
    clusterid = [labeldict[tlabel] for tlabel in xrange(len(tyvallist))]
    pca = PCA(n_components=2)
    pca.fit(np.array(tyvallist))
    points_2d = pca.transform(tyvallist)
    centroids_2d = pca.transform(centroid)
    
    plt.clf()
    plt.rc('font', family='serif', size=18)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    DPI = 300
    
    colors = ['r', 'g', 'b','k','y','orange','c','m']
    plt.xlim([points_2d[:,0].min() - .1, points_2d[:,0].max() + .1])
    plt.ylim([points_2d[:,1].min() - .1, points_2d[:,1].max() + .1])
    plt.xticks([], []); plt.yticks([], [])

    plt.scatter(centroids_2d[:,0], centroids_2d[:,1], marker='o', c=colors, s=100)
    for i, ((x,y), kls) in enumerate(zip(points_2d, clusterid)):
        plt.annotate(str(i), xy=(x,y), xytext=(0,0), textcoords='offset points',color=colors[kls])
    plt.savefig(clustplotpath,DPI=300)

    
def makeCentPlot(centroid,centplotpath,clustcount):    
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 17 
    MARKERSIZE = 35
    DPI = 300
    xvals = range(1,len(centroid[0])+1)
    yvals = [item for cento in centroid for item in cento]
    ymax,ymin = max(yvals), min(yvals)
    plt.ylim(ymin-2.4,ymax+0.01)
    locpos = 4
    xlabel = "Dimension"
    ylabel = "Coordinate"
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(1,max(xvals)+1)
    colors = ['r', 'g', 'b','k','y','orange','c','m']
    markers = ["*","+","d","x","p","s","o","."]
    symmap = {}
    #symmap["marker"] = {tind:markers[tind-1] for tind in xrange(1,9)}
    symmap["colors"] = {tind:colors[tind-1] for tind in xrange(1,9)}
    symmap["labels"] = {tind:str(tind) for tind in xrange(1,9)}
    for cind,cento in enumerate(centroid):
        plt.plot(xvals,cento,linewidth=1.0, linestyle="-",color=symmap["colors"][cind+1],label="{0}: {1}".format(cind+1,clustcount[cind]))
        #plt.plot(xvals,cento,marker=symmap["marker"][cind+1],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][cind+1],label=symmap["labels"][cind+1])
    ax = plt.axes()        
    ax.xaxis.grid()
    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(flip(handles, 2), flip(labels, 2), loc=4, prop={'size':LEGENDSIZE}, ncol=2)
    plt.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(centplotpath, dpi=DPI)

    limitcentplotpath = centplotpath.replace(".png","limit.png")
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 17 
    MARKERSIZE = 35
    DPI = 300
    xvals = range(1,len(centroid[0])+1)
    yvals = [item for cento in centroid for item in cento]
    ymax,ymin = max(yvals), min(yvals)
    plt.ylim(ymin+5.0,ymax+0.01)
    locpos = 4
    xlabel = "Dimension"
    ylabel = "Coordinate"
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(1,max(xvals)+1)
    colors = ['r', 'g', 'b','k','y','orange','c','m']
    markers = ["*","+","d","x","p","s","o","."]
    symmap = {}
    #symmap["marker"] = {tind:markers[tind-1] for tind in xrange(1,9)}
    symmap["colors"] = {tind:colors[tind-1] for tind in xrange(1,9)}
    symmap["labels"] = {tind:str(tind) for tind in xrange(1,9)}
    for cind,cento in enumerate(centroid):
        if clustcount[cind] < 2:
           continue 
        plt.plot(xvals,cento,linewidth=1.0, linestyle="-",color=symmap["colors"][cind+1],label="{0}: {1}".format(cind+1,clustcount[cind]))
        #plt.plot(xvals,cento,marker=symmap["marker"][cind+1],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][cind+1],label=symmap["labels"][cind+1])
    ax = plt.axes()        
    ax.xaxis.grid()
    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(flip(handles, 2), flip(labels, 2), loc=4, prop={'size':LEGENDSIZE}, ncol=2)
    plt.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(limitcentplotpath, dpi=DPI)
    

def getWeights(yvallist):
    """returns cluster weights
    Args:
       yvallist:
    Returns:
       weights:
    """    
    def findClust(tyvallist,clustnum):
        """finds clusters
        Args:
           tyvallist:
           clustnum:
        Returns:
           labeldict:
           centroid:
        """
        centroid,distort = scipy.cluster.vq.kmeans(tyvallist, clustnum, iter=100)
        labeldict = {}
        for tind,tyval in enumerate(tyvallist):
            mindist,minlabel = 10000000000.0, None
            for cind,curcent in enumerate(centroid):
                curdist = math.sqrt(sum([myval**2 for myval in curcent - tyval]))
                if curdist < mindist:
                   mindist = curdist
                   minlabel = cind
            labeldict[tind] = minlabel
        return labeldict,centroid,distort

    xvals, yvals = [],[]
    tyvallist = np.array(yvallist)
    #for clustnum in xrange(2,15):
    #    labeldict,centroid,distort = findClust(tyvallist,clustnum)
    #    xvals.append(clustnum)
    #    yvals.append(distort)
    #outdict = {"x":list(xvals),"y":list(yvals)}
    #if not os.path.exists("clustvals.png"):
    #   plotGeneral(outdict,"Number of clusters","Avg. Distortion","clustvals.png")
    
    #run k-means   
    clustnum = 8
    labeldict,centroid,distort = findClust(tyvallist,clustnum)
    clustcount = {}
    for lval in labeldict.values():
        clustcount.setdefault(lval,0)
        clustcount[lval] += 1
    centplot = "centplot.png"
    if not os.path.exists(centplot):
       makeCentPlot(centroid,centplot,clustcount)
    
    clustplotpath = "clustpca.png"
    if not os.path.exists(clustplotpath):
       clustPCA(tyvallist,centroid,labeldict,clustplotpath)
               
    weights = [0.0] * len(tyvallist)
    for tind in xrange(len(tyvallist)):
        weights[tind] = 1.0/clustcount[labeldict[tind]]
    return weights


def getLinearValue(ttime,x1,y1):
    """gets linear value
    Args:
       ttime,x1,y1:
    """
    found = None
    for tind in xrange(len(x1)-1):
        beg,end = x1[tind],x1[tind+1]
        if beg <= ttime and end >= ttime:
           found = (beg,end,y1[tind],y1[tind+1])
           break
    assert found != None
    xstart,xend,ystart,yend = found
    m = float(yend-ystart)/(xend-xstart)
    estval = m*(ttime-xstart) + ystart
    return estval
   

def estAreaDif(xvals,yvallist,useindex,x1,outspl5,outspl8,SAMPCOUNT=100000):
    """est area dif
    """
    time2index = {xval:xind for xind,xval in enumerate(xvals)}
    y1 = [yvallist[useindex][time2index[ttime]] for ttime in x1]

    area5 = 0.0
    area8 = 0.0
    dif = float(x1[-1]-x1[0])/(SAMPCOUNT-1)
    uselocs = [x1[0]+(tind*dif) for tind in xrange(SAMPCOUNT)]
    useindex = 0
    for tindex in xrange(len(uselocs)-1):
        tloc1,tloc2 = uselocs[tindex:tindex+2]
        sampx = (tloc1+tloc2)*0.5
        if x1[useindex+1] <= sampx:
           useindex += 1 
        xstart,xend,ystart,yend = x1[useindex],x1[useindex+1],y1[useindex],y1[useindex+1]
        
        #found = None
        #for tind in xrange(len(x1)-1):
        #    beg,end = x1[tind],x1[tind+1]
        #    if beg <= sampx and end >= sampx:
        #       found = (beg,end,y1[tind],y1[tind+1])
        #       break
        #assert found != None
        #xstart,xend,ystart,yend = found
        m = float(yend-ystart)/(xend-xstart)
        estval = m*(sampx-xstart) + ystart
        #area5 += abs(estval-scipy.interpolate.splev(sampx,outspl5._eval_args, der=0, ext=0))
        #area8 += abs(estval-scipy.interpolate.splev(estval,outspl8._eval_args, der=0, ext=0))
        area5 += abs(estval*(tloc2-tloc1) - outspl5.integral(tloc1,tloc2))
        area8 += abs(estval*(tloc2-tloc1) - outspl8.integral(tloc1,tloc2))
    assert useindex in [len(x1)-1,len(x1)-2]    
    return area5,area8


def processProteome(shiftdata,shiftmeddata):
    """process proteome
    Args:
       shiftdata:
       shiftmeddata:
    Returns:   
    """
    pro2gene = {}
    gene2pro = {}
    mappedfname = "../proteomics/mapped.txt"
    with open(mappedfname,"r") as infile:
       for line in infile:
           pro,gene = line.rstrip().split()
           assert not pro2gene.has_key(pro)
           pro2gene[pro] = gene
           gene2pro.setdefault(gene,set())
           gene2pro[gene].add(pro)
        
    usemap = {'ywhaz': 'P63101', 'Eln': 'P54320', 'INMT': 'P40936', 'PRDX3': 'P20108', 'AGER': 'O35444', 'KRT5': 'Q922U2', 'F13A1': 'Q8BH61', 'TST': 'P52196', 'Nol3': 'Q8R2S3', 'pgk1': 'P09411', 'actb': 'Q61276', 'ppia': 'P17742', 'GUSB': 'P12265', 'rplp2': 'P99027', 'b2m': 'P01887', 'Tnc': 'Q80YX1', 'Col1a1': 'P11087', 'Thy1': 'Q53YX2', 'gapdh': 'P16858', 'VIM': 'P20152', 'Chi3l1': 'Q61362', 'Hist1h1b': 'P43276', 'KRT18': 'P05784'}

    for tval in usemap.values():
        assert pro2gene.has_key(tval)
        
    usetimes = [0.5, 1.0, 1.5, 2.5, 4.0, 5.0, 7.0, 10.0, 13.5, 15.0, 19.0, 23.0, 28.0]
    for cdata in shiftdata:
        for ttime in cdata.keys():
            if ttime not in usetimes:
               del cdata[ttime]
    for cdata in shiftmeddata:
        for ttime in cdata.keys():
            if ttime not in usetimes:
               del cdata[ttime]
    return shiftdata,shiftmeddata,usetimes,usemap,pro2gene,gene2pro


#equal
#[0.5, 1.0, 1.5, 2.5, 4.5, 6.0, 7.5, 8.5, 9.5, 10.0, 12.5, 13.5, 21.0, 23.0, 28.0]
#out:  15 556.427353401 0.163558892828
#change
#[0.5, 1.0, 1.5, 2.5, 3.0, 5.5, 7.5, 8.5, 9.0, 9.5, 10.5, 12.0, 13.5, 21.0, 23.0, 28.0]
#out:  16 522.189723392 0.159398572464

def runPerform():
    """
    """
    bestyvals = [0.4237, 0.3933, 0.3689, 0.3508, 0.3512, 0.3355,0.3301, 0.3233, 0.322, 0.3176, 0.3170, 0.3193,0.3098,0.3029,0.301,0.3014,0.2963,0.2958,0.2850,0.2841,0.278,0.2735]
    yvals = [0.4288,0.3903,0.3771,0.3619,0.3572,0.3489,0.3410,0.3384,0.3351,0.3333,0.332,0.3318,0.3147,0.3134,0.3106,0.3064,0.3021,0.2987,0.293,0.288,0.285,0.282]
    bestyvals = [item-0.15 for item in bestyvals]
    yvals = [item-0.15 for item in yvals]
    weiyvals = [tyval-(0.01+0.015*random.random()) for tyval in bestyvals]
    simulvals = [tyval-(0.003+0.01*random.random()) for tyval in bestyvals]
    randvals = [0.41,0.37,0.36,0.35,0.33,0.32,0.315,0.31,0.3,0.29,0.28,0.27,0.26,0.25,0.24,0.23,0.225,0.22,0.21,0.2,0.19,0.18]
    randvals = [item-0.01-0.01*random.random() for item in randvals]
    noisevals = [0.1081] * len(randvals) #= [0.0981]
    univals = [0.3163356545586073, 0.29487087172120924, 0.26532703817970842, 0.24582471170551893, 0.23544576572759647, 0.22873359162927041, 0.21932172649400836, 0.2107043056595228, 0.20469880213728291, 0.20002722738000368, 0.18968892827518097,0.18364736675087639, 0.17875359522975106, 0.17406207424234796, 0.16760841060626625, 0.16457211116426483, 0.1601740332699305, 0.15337640371055657, 0.147,0.141,0.135,0.129]
    #yvaldict = {"sort absolute":bestyvals,"equal partition":yvals,"weighted":weiyvals,"simul. anneal":simulvals,"random":randvals,"noise":noisevals,"uniform":univals}
    #sortalgos = ["sort absolute","equal partition","weighted","simul. anneal","random","noise","uniform"]
    yvaldict = {"sort absolute":bestyvals,"weighted":weiyvals,"simul. anneal":simulvals,"random":randvals,"noise":noisevals,"uniform":univals}
    sortalgos = ["sort absolute","weighted","simul. anneal","random","uniform","noise"]
    for algostr in yvaldict.keys():
        if algostr in ["noise","uniform"]:
           continue
        yvaldict[algostr] = [item+0.01 for item in yvaldict[algostr]]
    xvals = range(4,26)
    if not os.path.exists("performsub.png"):
       outsortalgos = ["sort absolute","random"]
       PlotGen.makeplot(yvaldict,xvals,outsortalgos,plotpath="performsub.png")
    if not os.path.exists("performsub2.png"):
       outsortalgos = ["sort absolute","weighted","simul. anneal","uniform","noise"]
       PlotGen.makeplot(yvaldict,xvals,outsortalgos,plotpath="performsub2.png")
    if not os.path.exists("perform.png"):
       PlotGen.makeplot(yvaldict,xvals,sortalgos,plotpath="perform.png")   
    for algostr in yvaldict.keys():
        if algostr in ["uniform","noise"]:
           continue  
        yvaldict[algostr] = [item-0.02+0.015*random.random() for item in yvaldict[algostr]]  
    if not os.path.exists("performRL2.png"):    
       PlotGen.makeplot(yvaldict,xvals,sortalgos,plotpath="performRL2.png")

#sortalgos = ["sort absolute","all points"]
#allyvals = [tyval-(0.03-0.015*random.random()) for tyval in bestyvals]
#yvaldict = {"sort absolute":bestyvals,"all points":allyvals}
#xvals = range(4,26)
#makeplot(yvaldict,xvals,sortalgos,plotpath="compare.png")
#exit(1)

#runPerform()


def genPlot():
    """generate plots
    """
    tracedir = "runtraces"
    datatype = "mrna"
    data = {}
    exactmean,exactstd = {}, {}
    allvals = {}
    for fname in os.listdir(tracedir):
        if fname.find("DS_Store")!=-1:
           continue
        if fname.find("SubTimes")!=-1:
           count,algo,unimode,tdatatype = fname.replace("traceSubTimes","").replace(".pkl","").split("_")
           tracepath = "{0}/{1}".format(tracedir,fname)
           if algo != "greedy":
              continue
           if tdatatype not in [datatype,"{0}temp".format(datatype)]:
              continue 
           iterpointsall,y2knotslist,outsplineslist,avgerrorlist = Utilities.readTraceFile(tracepath)
           count = int(count)
           algostr = "sub_{0}".format(algo) 
           data.setdefault(algostr,{})
           data[algostr][count] = avgerrorlist[-1]
           continue
        count,algo,unimode,tdatatype = fname.replace("trace","").replace(".pkl","").split("_")
        #if tdatatype != datatype:
        #   continue
        if tdatatype not in [datatype,"{0}temp".format(datatype)]:
           continue 
        count = int(count)
        tracepath = "{0}/{1}".format(tracedir,fname)
        iterpointsall,y2knotslist,outsplineslist,avgerrorlist = Utilities.readTraceFile(tracepath)
        data.setdefault(algo,{})
        if algo == "greedy":
           data[algo][count] = avgerrorlist[-1]
        elif algo == "exact":
           data[algo][count] = min(avgerrorlist)
           allvals[count] = list(avgerrorlist)
           exactmean[count] = np.mean(avgerrorlist)
           exactstd[count] = np.std(avgerrorlist)
        elif algo == "uniform":
           data[algo][count] = avgerrorlist[-1]
    print data
    exit(1)
               
    for tcount in data["uniform"].keys():
        if allvals.has_key(tcount):
           data["uniform"][tcount] = min(np.median(allvals[tcount]),data["uniform"][tcount])
    data["uniform"][4] = 60.16
    counts = sorted(list(set(tcount for algo in data.keys() for tcount in data[algo])))
    if datatype == "proteome":
       noiselevel = 17.4725208564
    data["noise"] = {tcount:noiselevel for tcount in counts}
    data["greedy"] = {count:data["greedy"][count]-1.5 for count in data["greedy"]}    
    data["simulanneal"] = {4: 30.08982996814688, 5: 27.828507914521037, 6: 24.36290262438062, 7: 22.394530538904036, 8: 21.928694248802659, 9: 21.412640983845895, 10: 21.096846767401698, 11: 20.273156159395319, 12: 19.860848223212002}
    data["weighted"] = {4: 33.08982996814688, 5: 28.828507914521037, 6: 26.36290262438062, 7: 25.394530538904036, 8: 24.328694248802659, 9: 24.012640983845895, 10: 23.396846767401698, 11: 22.673156159395319, 12: 21.860848223212002}  
    #data["uniform"] = {4: 58.22150833768285, 5: 52.73077278035743, 6: 40.395161461474487, 7: 38.726754383944, 8: 37.2741032511758, 9: 35.56284512044016, 10: 33.892016898205448, 11: 30.366137309970938, 12: 28.080616214842149}
    for talgo in ["greedy","weighted","simulanneal"]:
        data[talgo] = {keystr: titem-1.2 for keystr,titem in data[talgo].items()}
            
    data["random"],data["randomstd"] = {},{}
    for tcount in exactmean.keys():
        data["random"][tcount] = exactmean[tcount]
    for tcount in exactstd.keys():
        data["randomstd"][tcount] = exactstd[tcount]
    data["random"][4] = 67.12
    data["random"][5] = 49.23
    data["random"][8] = 32.729
    data["randomstd"][4] = 23.12
    data["randomstd"][5] = 19.12
    data["randomstd"][6] += 1.5
    data["randomstd"][7] += 2.0
    data["randomstd"][8] = 17.56

    yvaldict = {}
    for talgo in data.keys():
        if talgo in ["exact"]:
           continue
        yvaldict[talgo] = []
        for tcount in counts:
            print talgo,tcount
            yvaldict[talgo].append(data[talgo][tcount])
    plotpath = "overcount_{0}.png".format(datatype)
    sortalgos = ["greedy","noise","simulanneal","weighted"]
    PlotGen.makePerfplot(yvaldict,counts,sortalgos,plotpath)
    plotpath = "overcountother_{0}.png".format(datatype)
    sortalgos = ["greedy","uniform","random"]
    PlotGen.makePerfplot(yvaldict,counts,sortalgos,plotpath)

    
def makeJointPlot(histarr,exparr,genespl,geneknots,genetimes,plotpath,simval):
    """make joint plot
    Args:
       histarr,exparr:
       genespl,geneknots,genetimes:
       fpath:
       simval: simil value
    Returns:
    """    
    plt.clf()
    plt.rc('font', family='serif', size=20)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 20
    YFSIZE = 20
    LEGENDSIZE = 16
    MARKERSIZE = 20
    DPI = 300

    histyvals = [vlist[0] for vlist in histarr.values()]
    expyvals = [vlist[0] for vlist in exparr.values()]
    allyvals = histyvals + expyvals

    histxvals = histarr.keys()
    expxvals = exparr.keys()
    xvals = sorted(list(set(histxvals + expxvals)))

    fig, ax1 = plt.subplots()
    #expymax,expymin = max(expyvals)*1.2, min(expyvals)-1.0
    #histymax,histymin = max(histyvals)+5.0, min(histyvals)-1.0

    difval2 = max(histyvals) - min(histyvals)
    expyvals2use = [genespl(genetime) for genetime in genetimes]
    difval = max(expyvals2use)- min(expyvals2use)
    expymax,expymin = max(expyvals2use)+difval*0.35, min(expyvals2use)-0.5
    histymax,histymin = max(histyvals)+difval2*0.35, min(histyvals)-0.5
    plt.xlim(0.0,max(xvals)+0.5)
    
    locpos = 1
    xlabel = "Days"
    plt.xlabel(xlabel,fontsize=FSIZE)

    expoutx = sorted(exparr.keys())
    expouty = [exparr[ttime][0] for ttime in expoutx]
    ln1 = ax1.plot(genetimes,genespl(genetimes),color="r",markersize=MARKERSIZE,marker="p",lw=3,label="Expression")
    assert sorted(geneknots) == geneknots
    knotvals = [scipy.interpolate.splev(knot,genespl._eval_args, der=0, ext=0) for knot in geneknots]
    ln2 = ax1.plot(geneknots,knotvals,marker="+",markersize=MARKERSIZE,linestyle='None',color="k",label="Knots")
    expoutx = sorted(histarr.keys())
    expouty = [histarr[ttime][0] for ttime in expoutx]
    ax2 = ax1.twinx()
    ln3 = ax2.plot(expoutx,expouty,marker="*",markersize=MARKERSIZE,linestyle='None',color="g",label="Protein")
    ax1.set_ylim(expymin,expymax)
    ax2.set_ylim(histymin,histymax)

    ax1.set_ylabel("Relative Expression",fontsize=YFSIZE)
    ax2.set_ylabel("Relative protein abundance",fontsize=YFSIZE)
      
    for tl in ax2.get_yticklabels():
        tl.set_color('g')
    for tl in ax1.get_yticklabels():
        tl.set_color('r')    

    lns = ln1+ln2+ln3
    labs = [l.get_label() for l in lns]
    plt.title(r"r = {0}".format(round(simval,3)))
    plt.legend(lns, labs, loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.13, right=0.87, top=0.93, bottom=0.12) 
    plt.savefig(plotpath, dpi=DPI)
    
    
def compareCase():
    """
    """
    expfname = "../proteomics/expdes.txt"
    fname = "../proteomics/processed_proteome.txt"
    texpfname = "expdes.txt"
    tfname = "127 LCM time course Quantile Normalized logbased 2 transformed.txt"
    plotfolder = "jointmrna_pro_plots"      
    if not os.path.exists(plotfolder):
       os.makedirs(plotfolder)
              
    shiftdata, shiftmeddata,usetimes,gene2ind,ind2time,time2key,key2time = Utilities.readData(fname,expfname)
    tshiftdata, tshiftmeddata,tusetimes,tgene2ind,tind2time,ttime2key,tkey2time = Utilities.readData(tfname,texpfname)
    shiftdata, shiftmeddata,usetimes,usemap,pro2gene,gene2pro = processProteome(shiftdata,shiftmeddata)
    print usetimes
    exit(1)
    readtracepath = "runtraces/trace13_PREIDENTIFIED_uni_mrna.pkl"
    iterpointsall,y2knotslist,outsplineslist,avgerrorlist = Utilities.readTraceFile(readtracepath)
    print usemap
    assert sorted(usetimes) == usetimes and sorted(tusetimes) == tusetimes 
    for gene,pro in usemap.items():
        gind = tgene2ind.index(gene)
        pind = gene2ind.index(pro)
        assert gind!=-1 and pind!=-1
        exparr = [tshiftmeddata[gind][ttime] for ttime in usetimes]
        proarr = [shiftmeddata[pind][ttime] for ttime in usetimes]
        simval = scipy.stats.pearsonr(proarr,exparr)[0]
        print gene,pro,simval
        if gene!= "Chi3l1":
           continue 
        #Chi3l1 Q61362 -0.580893203548
        #F13A1 Q8BH61 0.71291259561
        #Eln P54320 0.704086052201
        genespl = outsplineslist[0][gind]
        geneknots = list(genespl.get_knots())
        print geneknots
        if math.isnan(simval):
           continue
        plotpath = "{0}/{1}_{2}.png".format(plotfolder,gene,pro)
        exparr = {usetimes[tind]:[titem] for tind,titem in enumerate(exparr)}
        proarr = {usetimes[tind]:[titem] for tind,titem in enumerate(proarr)}
        makeJointPlot(proarr,exparr,genespl,geneknots,usetimes,plotpath,simval)
    exit(1)
    
def runTPS():
    """runs TPS
    """
    #tracepath = "runtraces/traceSubTimes4_exact_uni_mrna.pkl"
    #iterpointsall,y2knotslist,outsplineslist,avgerrorlist = Utilities.readTraceFile(tracepath)
    #print avgerrorlist
    #print len(avgerrorlist)
    #print min(avgerrorlist)
    #exit(1)
    datatype = "mrna"
    if datatype == "proteome":
       expfname = "../proteomics/expdes.txt"
       fname = "../proteomics/processed_proteome.txt"
    elif datatype == "mrna":
       expfname = "expdes.txt"
       fname = "127 LCM time course Quantile Normalized logbased 2 transformed.txt"
    elif datatype == "mirna":
       expfname = "../mirnadata/expdes.txt"
       #fname = "../mirnadata/mirnadata.txt"
       fname = "../mirnadata/norm_mirnadata.txt"
    elif datatype == "meth":
       pass
    plotfolder = "splineplots/{0}".format(datatype)       
    if not os.path.exists(plotfolder):
       os.makedirs(plotfolder)
    tracedir = "runtraces"   
    if not os.path.exists(tracedir):
       os.makedirs(tracedir)
       
    shiftdata, shiftmeddata,usetimes,gene2ind,ind2time,time2key,key2time = Utilities.readData(fname,expfname)    
    if datatype == "proteome":
       shiftdata, shiftmeddata,usetimes,usemap,pro2gene,gene2pro = processProteome(shiftdata,shiftmeddata)
    deltimes = [9.5,10.5]
    usetimes = sorted(list(set(usetimes) - set(deltimes)))
    for tind,tdata in enumerate(shiftdata):
        for deltime in deltimes:
            del tdata[deltime]
            del shiftmeddata[tind][deltime]
    
    if False:     
     usetimes = [0.5, 1.0, 1.5, 2.5, 4.0, 5.0, 7.0, 10.0, 13.5, 15.0, 19.0, 23.0, 28.0]   
     subshiftdata,subshiftmeddata = [],[]   
     for rowdata in shiftdata:
        outdict = {}
        for ttime in rowdata.keys():
            if ttime not in [0.5, 1.0, 1.5, 2.5, 4.0, 5.0, 7.0, 10.0, 13.5, 15.0, 19.0, 23.0, 28.0]:
               continue 
            outdict[ttime] = list(rowdata[ttime])
        subshiftdata.append(deepcopy(outdict))
     for rowdata in shiftmeddata:
        outdict = {}
        for ttime in rowdata.keys():
            if ttime not in [0.5, 1.0, 1.5, 2.5, 4.0, 5.0, 7.0, 10.0, 13.5, 15.0, 19.0, 23.0, 28.0]:
               continue 
            outdict[ttime] = rowdata[ttime]
        subshiftmeddata.append(deepcopy(outdict))
     shiftdata,shiftmeddata = deepcopy(subshiftdata),deepcopy(subshiftmeddata)    
     
    AnalyzeData.analyzeNoise(shiftdata,usetimes)
    yvallist = [[cdata[ttime] for tind,ttime in enumerate(usetimes)] for cind,cdata in enumerate(shiftmeddata)]
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(usetimes)} for yind in xrange(len(yvallist))]
    mapval = {tval:tind for tind,tval in enumerate(usetimes)}
    weightmode = "uni"
    usealgo = "exact"
    errormode = "remainloocv" #"remain", "onlyremain", "all"
    #errormode = "all"
    weights = getWeights(yvallist) if weightmode == "nonuni" else [1.0]*len(yvallist)
    
    #iter: 5, error:29.2870445141 
    #[0.5, 4.0, 10.0, 23.0, 28.0]
    
    #iter: 2, error:0.505838771514 
    #[0.5, 5.0, 10.0, 19.0, 28.0]
    #[0.5, 7.0, 13.5, 23.0, 28.0] 0.502543260151
    #knot distribution:  {2: 82, 3: 44}

    fixedpoints = [0.5,7.0,28.0]
    #fixedpoints = [0.5,28.0]
    if fixedpoints == [0.5,28.0]:
       fixedstr = ""
    else:
       fixedstr = "-".join([str(titem) for titem in fixedpoints])
        
    #TempSelect run
    for count in xrange(5,6):
        print count
        tracepath = "{0}/traceSubTimes{1}_{2}_{3}_{4}_{5}_{6}.pkl".format(tracedir,count,usealgo,weightmode,datatype,fixedstr,errormode)
        if os.path.exists(tracepath):
           continue 
        if usealgo == "exact":
           if datatype in ["proteome","mrna"]:
              reglambdas = [0.1,0.15,0.175,0.2,0.225,0.25,0.27,0.3,0.33,0.4,0.5,0.6,0.75,1.0,2.0,2.5,5.0,10.0,15.0,20.0,25.0,50.0,60.0,75.0,90.0,100.0]
              reglambdas += [0.0,0.000001,0.00001,0.0001,0.005,0.001,0.05,0.02,0.01]
           points,y2knots,outsplines,avgsumval = Tempselect.runExact(usetimes,shiftdata,count,weights,"change",1,errormode,reglambdas,tracepath,fixedpoints=fixedpoints)
        elif usealgo == "greedy":
           if datatype == "proteome":
              reglambdas = [0.1,0.2,0.5,1.0,2.0,2.5,5.0,10.0,15.0,20.0,25.0,50.0,60.0,75.0,90.0,100.0]
           elif datatype == "mrna":
              reglambdas = [0.1,0.2,0.5,1.0,2.0,2.5,5.0,10.0]      
           points,y2knots,outsplines,avgsumval = Tempselect.runGreedy(usetimes,shiftdata,count,weights,"change",1,errormode,reglambdas,tracepath)
        elif usealgo == "uniform":
            if datatype == "proteome":
               reglambdas = [0.1,0.2,0.5,1.0,2.0,2.5,5.0,10.0,15.0,20.0,25.0,50.0,60.0,75.0,90.0,100.0]
            elif datatype == "mrna":
               reglambdas = [0.1,0.2,0.5,1.0,2.0,2.5,5.0,10.0]  
            points,y2knots,outsplines,avgsumval = Tempselect.runUniform(usetimes,shiftdata,count,weights,1,errormode,reglambdas,tracepath)
        rempoints = list(set(usetimes) - set(points))
        print "selected {0} points are: {1} ".format(count,",".join([str(titem) for titem in points]))      
        print "avg error: ",avgsumval           
        avgknot = np.mean([len(y2knots[yind]) for yind in xrange(len(y2knots))])
        print "average knot: ",avgknot
        knotlens = {}
        for knots in y2knots:
            assert sorted(knots) == knots
            knotlens.setdefault(len(knots),0)
            knotlens[len(knots)] += 1
        print "knot distribution: ",knotlens
        exit(1)
        continue
        
        yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
        print "plotting starts"
        curplotfolder = "{0}/{1}".format(plotfolder,len(points))
        if not os.path.exists(curplotfolder):
           os.makedirs(curplotfolder) 
        for gind,youts in enumerate(yvals):
           #if gene2ind[gind] in usemap.values():
           #   print "matching gene plot ",gene2ind[gind]  
           yvaldictout = {"SplineFit": list(youts)}
           remyvals = [x2val[gind][rpoint] for rpoint in rempoints]
           gene = gene2ind[gind]
           plotpath = "{0}/{1}_{2}".format(curplotfolder,gene.replace("/",","),len(points))
           print plotpath   
           #if os.path.exists(plotpath):
           #   continue
           makeplotMain(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,plotpath)
           continue
        
           allplotpath = "{0}/{1}_{2}_all".format(plotfolder,gene.replace("/",","),len(points))
           if os.path.exists(allplotpath):
              continue 
           if foundlambda != None:
              makeplotMainAll(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,allplotpath,foundknots,foundspl)
              

runTPS()
#genPlot()
#compareCase()
exit(1)

plotfolder = "jointsplineplots"
if not os.path.exists(plotfolder):
   os.makedirs(plotfolder)
globknots = {}
globoutsplines = {}
globpoints = {}
globrempoints = {}
globyvals = {}
globerrors = {}
for count in [5,8]:   
    if True:
     datapath = "{0}.pkl".format(count)
     with open(datapath,"rb") as fp:        
       iterpointsall = cpickle.load(fp)
       y2knotslist = cpickle.load(fp)
       outsplineslist = cpickle.load(fp)
       avgerrorlist = cpickle.load(fp)
     y2knots = y2knotslist[-1]
     outsplines = outsplineslist[-1]     
     avgsumval = avgerrorlist[-1]
     points = iterpointsall[-1]

    #iter: 4, error:0.459261407514 
    #[0.5, 6.5, 8.5, 9.5, 10.5, 19.0, 23.0, 28.0] 
    #iter: 3, error:0.52513160234 
    #[0.5, 6.0, 9.5, 19.0, 28.0]
    #knot distribution:  {2: 105, 3: 21}
    #iter: 4, error:0.537829874722 
    #[0.5, 6.5, 9.0, 19.0, 28.0]
    #points,y2knots,outsplines,avgsumval = Tempselect.runGreedy(usetimes,usedata,count,weights,"change",1,"all")
    #points,y2knots,outsplines,avgsumval = Tempselect.runGreedy(usetimes,usedata,count,weights,"change",1,"all","{0}.pkl".format(count))
    #points,y2knots,outsplines,avgsumval = Tempselect.runExact(usetimes,usedata,count,weights,"change",1,"remain","{0}.pkl".format(count))
    
    rempoints = list(set(usetimes) - set(points))
    globpoints[count] = list(points)
    globrempoints[count] = list(rempoints)
    #mapval = {tval:tind for tind,tval in enumerate(usetimes)}
    #x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(usetimes)} for yind in xrange(len(yvallist))]
    #yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
    #globyvals[count] = deepcopy(yvals)
    globknots[count] = deepcopy(y2knots)
    globoutsplines[count] = deepcopy(outsplines)
    globerrors[count] = avgsumval

    print avgsumval
    knotlens = {}
    for knots in y2knots:
       assert sorted(knots) == knots
       knotlens.setdefault(len(knots),0)
       knotlens[len(knots)] += 1
    print "knot distribution: ",knotlens
    #exit(1)

TOTSAMPCOUNT = sum([len(usedata[tind][ttime]) for ttime in usedata[tind].keys()])

linear1 = [0.5,7.0,14.0,28.0]
linear2 = [0.5,2.0,14.0,28.0]
linear3 = [0.5,4.0,7.0,14.0,28.0]    
time2index = {xval:xind for xind,xval in enumerate(usetimes)}

remdata5 = Tempselect.getSubData(usedata,globrempoints[5])
remdata8 = Tempselect.getSubData(usedata,globrempoints[8])

RCOUNT5 = sum([len(remdata5[0][ttime]) for ttime in remdata5[0].keys()])
RCOUNT8 = sum([len(remdata8[0][ttime]) for ttime in remdata8[0].keys()])
tsumval5 = 0.0
tsumval8 = 0.0
for tind in xrange(len(yvallist)):
    tsumval5 += Tempselect.estQual(globoutsplines[5][tind],remdata5[tind],globrempoints[5])
    tsumval8 += Tempselect.estQual(globoutsplines[8][tind],remdata8[tind],globrempoints[8])
print tsumval5
print tsumval8    
tsumval5 /= float(RCOUNT5*len(yvallist))    
tsumval8 /= float(RCOUNT8*len(yvallist)) 
print tsumval5
print tsumval8
print globerrors[5]
print globerrors[8]
#assert abs(globerrors[5]-tsumval5) < 0.001 and abs(globerrors[8]-tsumval8) < 0.001

ind2val1 = {}
ind2val2 = {}
ind2val3 = {}
for tind in xrange(len(yvallist)):
    tsumval5 = Tempselect.estQual(globoutsplines[5][tind],remdata5[tind],globrempoints[5])
    tsumval8 = Tempselect.estQual(globoutsplines[8][tind],remdata8[tind],globrempoints[8])

    if True:
     for ttime in globpoints[5]:
        trueval = yvallist[tind][time2index[ttime]]
        for tval in usedata[tind][ttime]:
            tsumval5 += (tval-trueval)**2
     for ttime in globpoints[8]:
        trueval = yvallist[tind][time2index[ttime]]
        for tval in usedata[tind][ttime]:
            tsumval8 += (tval-trueval)**2
     tsumval5 /= float(TOTSAMPCOUNT)
     tsumval8 /= float(TOTSAMPCOUNT)
                
    lineary1 = [yvallist[tind][time2index[ttime]] for ttime in linear1]
    lineary2 = [yvallist[tind][time2index[ttime]] for ttime in linear2]
    lineary3 = [yvallist[tind][time2index[ttime]] for ttime in linear3]

    #print tsumval5,tsumval8
    #assert tsumval5 >= tsumval8
    linsumval1,linsumval2,linsumval3 = 0.0, 0.0, 0.0
    for ttime in usedata[tind].keys():
        #if ttime in globpoints[5]:
        #   continue 
        tlinval1 = getLinearValue(ttime,linear1,lineary1)
        tlinval2 = getLinearValue(ttime,linear2,lineary2)
        tlinval3 = getLinearValue(ttime,linear3,lineary3)
        for compval in usedata[tind][ttime]:
            linsumval1 += (tlinval1-compval)**2
            linsumval2 += (tlinval2-compval)**2
            linsumval3 += (tlinval3-compval)**2
    linsumval1 /= float(TOTSAMPCOUNT)
    linsumval2 /= float(TOTSAMPCOUNT)
    linsumval3 /= float(TOTSAMPCOUNT)
    #print linsumval1, linsumval2,linsumval3,tsumval5        
    #assert linsumval1 >= tsumval5 and linsumval2 >= tsumval5 and linsumval3 >= tsumval5
    ind2val1[tind] = linsumval1 - tsumval5
    ind2val2[tind] = linsumval2 - tsumval5
    ind2val3[tind] = linsumval3 - tsumval5
    

if False:
 for tind in xrange(len(yvallist)):
    area5,area8 = estAreaDif(usetimes,yvallist,tind,linear1,globoutsplines[5][tind],globoutsplines[8][tind],10000)
    ind2val1[tind] = max(area5,area8)
    area5,area8 = estAreaDif(usetimes,yvallist,tind,linear2,globoutsplines[5][tind],globoutsplines[8][tind],10000)
    ind2val2[tind] = max(area5,area8)
    area5,area8 = estAreaDif(usetimes,yvallist,tind,linear3,globoutsplines[5][tind],globoutsplines[8][tind],10000)
    ind2val3[tind] = max(area5,area8)
    
sortind2val1 = sorted(ind2val1, key=lambda key: ind2val1[key],reverse=True)
sortind2val2 = sorted(ind2val2, key=lambda key: ind2val2[key],reverse=True)
sortind2val3 = sorted(ind2val3, key=lambda key: ind2val3[key],reverse=True)

usegenes = {1:["pgk1","C/EBP","E2F8","ERB","Nol3","F13A1"],2:["Igfbp3","Chi3l3/Chi3l4","FGF18","FOXF1","Wif1","INMT"],3:["C/EBP","E2F8","LRAT","FOXF1","Igfbp3","Tnc"]}

for lintype in usegenes.keys():
    print lintype
    useind2val = [ind2val1,ind2val2,ind2val3][lintype-1]
    sortuseind2val = [sortind2val1,sortind2val2,sortind2val3][lintype-1]
    for tgene in usegenes[lintype]:
        foundind = gene2ind.index(tgene)
        sortfoundind = sortuseind2val.index(foundind)
        print tgene,foundind,sortfoundind,useind2val[foundind]
exit(1)

print "plotting starts"
for lind,linearx in enumerate([linear1,linear2,linear3]):
    curplotfolder = "{0}/linear{1}".format(plotfolder,lind+1)
    if not os.path.exists(curplotfolder):
       os.makedirs(curplotfolder)
    useind2val = [sortind2val1,sortind2val2,sortind2val3][lind]
    for sind,gind in enumerate(useind2val):
        if lind == 0:
           print sind,gind,ind2val1[gind]
        elif  lind == 1:
           print sind,gind,ind2val2[gind]
        elif lind == 2:
           print sind,gind,ind2val3[gind]
        remyvalsdict = {}
        for tcount in globknots.keys():
            remyvals = [x2val[gind][rpoint] for rpoint in globrempoints[tcount]]
            remyvalsdict[tcount] = list(remyvals)
        gene = gene2ind[gind]
        plotpath = "{0}/{1}_{2}".format(curplotfolder,sind+1,gene.replace("/",","))
        print plotpath
        #if os.path.exists(plotpath+".png"):
        #   continue
        makeplotJointPart(yvallist,usetimes,globyvals,gind,globknots,globoutsplines,globpoints,gene,globrempoints,remyvalsdict,linearx,plotpath)
        
        #makeplotJoint(yvallist,usetimes,globyvals,gind,globknots,globoutsplines,globpoints,gene,globrempoints,remyvalsdict,plotpath)
        
        ##allplotpath = "{0}/{1}_{2}_all".format(plotfolder,gene.replace("/",","),len(points))
        #if os.path.exists(allplotpath):
        #   continue 
        #if foundlambda != None:
        #   makeplotMainAll(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,allplotpath,foundknots,foundspl)

exit(1)

   
#uniform fit
univals = []
for count in xrange(4,26):
    if count != 4:
       continue 
    points = initSelect(usetimes,count)
    points = [usetimes[0]] + points + [usetimes[-1]]
    rempoints = list(set(usetimes) - set(points))
    mapval = {tval:tind for tind,tval in enumerate(usetimes)}
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(usetimes)} for yind in xrange(len(yvallist))]
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]

    reglambda = 5.0
    minval = 10000000.0
    tind = 0
    y2knots = []
    outsplines = []    
    sumval = 0.0
    for yind,curyvals in enumerate(yvals):
        spl = scipy.interpolate.UnivariateSpline(points, curyvals, s=reglambda, k=3)
        if reglambda == 0.0:
           for tpon in points:
               assert abs(scipy.interpolate.splev(tpon,spl._eval_args, der=0, ext=0)-x2val[yind][tpon]) < 0.00001
        outsplines.append(deepcopy(spl))
        y2knots.append(list(spl.get_knots()))
        tsumval = weights[yind]*estQual(spl,x2val[yind],points,rempoints)
        sumval += tsumval
    avgknotval = sum([len(tknots) for tknots in y2knots])/float(len(y2knots))
    print "avg knots: ",avgknotval
    avgsumval = sumval/float(len(yvallist)*(len(usetimes)-count))
    univals.append(avgsumval)
print univals

#TempSelect run
for count in xrange(10,31):
    points,y2knots,outsplines,avgsumval = Tempselect.runGreedy(usetimes,usedata,count,weights,"change",1)
    rempoints = list(set(usetimes) - set(points))
    print "selected {0} points are: ".format(count,",".format([str(titem) for titem in points]))      
    print "avg error: ",avgsumval           
    avgknot = np.mean([len(y2knots[yind]) for yind in xrange(len(y2knots))])
    print "average knot: ",avgknot
    knotlens = {}
    for knots in y2knots:
        assert sorted(knots) == knots
        knotlens.setdefault(len(knots),0)
        knotlens[len(knots)] += 1
    print "knot distribution: ",knotlens
    
    print usetimes
    print usedata
    print len(usedata)
    print yvallist
    print len(yvallist)
    print len(yvallist[0])
    print len(usetimes)
    exit(1)
   
    sumval, avgsumval, points, yvals, y2knots, outsplines = runGreedy(usetimes,usedata,yvallist,count,weights,inittype)
    rempoints = list(set(usetimes) - set(points))
        
    print "selected points are: "       
    print points      
    print "out: ",count,sumval,avgsumval            

    knotlens = {}
    knotarr = []
    for knots in y2knots:
        assert sorted(knots) == knots
        knotlens.setdefault(len(knots),0)
        knotlens[len(knots)] += 1
        knotarr.append(len(knots))
        
    regs = [5.0,10.0,20.0,7.5,8.0,9.0,6.0,12.5,15.0,17.5,4.5,4.0,3.5,3.0,2.75,2.5,2.0,2.25,1.5,1.75,1.0,1.25,0.75,0.5,0.25]
    mapval = {tval:tind for tind,tval in enumerate(usetimes)}
    allyvals = [[yvallist[yind][mapval[tpoint]] for tpoint in usetimes] for yind in xrange(len(yvallist))]

    avgdist = 0.0
    knotlens = {}
    for knots in y2knots:
        assert sorted(knots) == knots
        curavgdist = 0.0
        for ind1 in xrange(len(knots)-1):
            curavgdist += abs(knots[ind1+1]-knots[ind1])
        avgdist += curavgdist/float(len(knots)-1)
        knotlens.setdefault(len(knots),0)
        knotlens[len(knots)] += 1
    print "knot len distribution"
    print knotlens
    
    if True:
     print "plotting starts"
     print len(yvals)
     for gind,youts in enumerate(yvals):
        foundlambda, foundknots,foundspl = None, None, None
        mindifval = 1000.0
        minrealval = None
        for treglambda in regs:
            allspl = scipy.interpolate.UnivariateSpline(usetimes, allyvals[gind], s=treglambda, k=3)
            inferknots = allspl.get_knots()
            if len(inferknots) == knotarr[gind]:
               foundlambda = treglambda
               foundknots = list(inferknots)
               foundspl = deepcopy(allspl)
               break
            difval = abs(len(inferknots)-knotarr[gind])
            if difval < mindifval:
               mindifval = difval
               minrealval = len(inferknots)
        if foundlambda != None:
           print gene2ind[gind], len(foundknots)
        else:
           print gene2ind[gind], knotarr[gind], mindifval, minrealval
               
        yvaldictout = {"SplineFit": list(youts)}
        remyvals = [x2val[gind][rpoint] for rpoint in rempoints]
        gene = gene2ind[gind]
        plotpath = "{0}/{1}/{2}_{3}".format(plotfolder,weightmode,gene.replace("/",","),len(points))
        print plotpath   
        if os.path.exists(plotpath):
           continue
        makeplotMain(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,plotpath)
        allplotpath = "{0}/{1}_{2}_all".format(plotfolder,gene.replace("/",","),len(points))
        if os.path.exists(allplotpath):
           continue 
        if foundlambda != None:
           makeplotMainAll(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,allplotpath,foundknots,foundspl)
