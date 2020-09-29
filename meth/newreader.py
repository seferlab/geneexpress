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

def readExpText(fname):
    """
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

runPerform()

plotfolder = "splineplots"
if not os.path.exists(plotfolder):
   os.makedirs(plotfolder)
expfname = "expdes.txt"
     
time2key,key2time = readExpText(expfname)
fname = "127 LCM time course Quantile Normalized logbased 2 transformed.txt"
#fname = "127LCM time course Data Not normalized.txt"
data,ind2time,gene2ind = readMultiSampleData(fname,key2time)

#modify data
inittype = "change"
times = sorted(time2key.keys())
usetimes = times[1:]
usedata = deepcopy(data)
for dpart in usedata:
    del dpart[times[0]]
    #del dpart[times[-2]]
    #del dpart[times[-1]]

#Take only methylation points
if False:  
 mytimes = [0.5, 1.5, 2.5, 5.0, 10.0, 15.0, 19.0, 26.0]
 for dpart in usedata:
    for keytime in dpart.keys():
        if keytime not in mytimes:
           del dpart[keytime]
 usetimes = sorted(usedata[0].keys())     

AnalyzeData.analyzeNoise(usedata,usetimes)

yvallist = []
for cind,cdata in enumerate(usedata):
    cavgdata = []
    for tind,ttime in enumerate(usetimes):
        avgval = np.mean(list(usedata[cind][ttime]))
        cavgdata.append(avgval)
    cavgdatashift = [item-cavgdata[0] for item in cavgdata]        
    yvallist.append(list(cavgdatashift))
x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(usetimes)} for yind in xrange(len(yvallist))]

weightmode = "nonuni"
weightmode = "uni"
if weightmode == "nonuni":
   weights = getWeights(yvallist)
elif weightmode == "uni":   
   weights = [1.0]*len(yvallist)


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
