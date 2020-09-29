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
import os
import sys
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
    
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in times] 
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
    plt.plot(times,outspl(times),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')

    allyvals = [scipy.interpolate.splev(knot,allspline._eval_args, der=0, ext=0) for tknot in allknot]
    plt.plot(times,allspline(times),'y',lw=6,label='All')
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
    
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in times] 
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
    plt.plot(times,outspl(times),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
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

def runGreedy(times,usedata,count,weights,initmode="change"):
    """runs greedy method
    Args:
       times,usedata:
       count,weights:
    Returns:   
    """
    def changeRange(knots,minval,maxval):
        """change range
        """
        return [(maxval-minval)*tknot + minval for tknot in knots]
    def convertForm(curtimes,usedata):
        """converts form
        Args:
           curtimes,usedata:
        Returns:
           yvallist:   
        """        		
        yvallist = []
        for cind,cdata in enumerate(usedata):
            calldata = []
            for tind,ttime in enumerate(curtimes):
                calldata.extend(usedata[cind][ttime])
            yvallist.append(list(calldata))
            #assert len(calldata) == len(usetimes)
        return yvallist

    globtimes = set()
    for cind,cdata in enumerate(usedata):
        globtimes |= set(usedata[cind].keys())
    for cind,cdata in enumerate(usedata):
        assert len(globtimes) == len(usedata[cind].keys())
    #usetimes = [ttime for tind,ttime in enumerate(times) for tind2 in xrange(len(usedata[0][ttime]))]
    #yvallist = convertForm(times,usedata) #??
    
    assert initmode in ["change","equal"] and times == sorted(times)
    if initmode == "equal":
       points = initSelect(times,count)
    elif initmode == "change":
       points = initSelectAvgchange(times,count,usedata)
    
    points = [times[0]] + points + [times[-1]]
    points = sorted(points)
    multipoints = [point for point in points for tind2 in xrange(len(usedata[0][point]))]
    rempoints = list(set(times) - set(points))
    multirempoints = [rpoint for rpoint in rempoints for tind2 in xrange(len(usedata[0][rpoint]))]
    #mapval = {tval:tind for tind,tval in enumerate(times)}
    #x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    #yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
    yvals = convertForm(points,usedata)
    remyvals = convertForm(rempoints,usedata)
         
    minval = 10000000.0
    tind = 0
    y2knots = []
    sumval = 0.0
    for yind,curyvals in enumerate(yvals):
        cursumval,time2err,knots = runRFit(multipoints,curyvals,multirempoints,remyvals[yind],usedata[yind])
        modknots = changeRange(knots,times[0],times[-1])
        y2knots.append(list(modknots))
        tsumval = weights[yind]*cursumval
        sumval += tsumval            
    while True:
        print sumval
        minsol = None
        for addpoint in rempoints:
            for delpoint in points:
                if delpoint in [times[0],times[-1]]:
                   continue 
                newrempoints, newpoints = list(rempoints), list(points)
                newrempoints.remove(addpoint)
                newrempoints.append(delpoint)
                newpoints.remove(delpoint)
                newpoints.append(addpoint)
                newpoints = sorted(newpoints) 

                newmultipoints = [point for point in newpoints for tind2 in xrange(len(usedata[0][point]))]
                newmultirempoints = [rpoint for rpoint in newrempoints for tind2 in xrange(len(usedata[0][rpoint]))]
                newyvals = convertForm(newpoints,usedata)
                newremyvals = convertForm(newrempoints,usedata) 
                globsumval = 0.0
                cury2knots = []
                for yind,curyvals in enumerate(newyvals):
                    cursumval,time2err,knots = runRFit(newmultipoints,curyvals,newmultirempoints,newremyvals[yind],usedata[yind])
                    modknots = changeRange(knots,times[0],times[-1])
                    cury2knots.append(list(modknots))
                    tsumval = weights[yind]*cursumval
                    globsumval += tsumval
                print addpoint,delpoint,globsumval    
                if globsumval < sumval:
                   sumval = globsumval
                   minsol = (addpoint,delpoint)
                   y2knots = deepcopy(cury2knots)
                   
                    #inyvals = [x2val[yind][rpoint] for rpoint in newpoints]
                    #inspl = scipy.interpolate.UnivariateSpline(newpoints, inyvals, s=reglambda, k=3)
                    
                    #spl2 = scipy.interpolate.InterpolatedUnivariateSpline(sortednewpoints, sortedinyvals,k=3)
                    #for tpon in newpoints:
                    #    assert abs(spl2.__call__(tpon)-x2val[yind][tpon]) < 0.00001
                    #    #assert abs(scipy.interpolate.splev(tpon,spl2._eval_args, der=0, ext=0)-x2val[yind][tpon]) < 0.00001

                    #if reglambda == 0.0:
                    #   for tpon in newpoints:
                    #       assert abs(scipy.interpolate.splev(tpon,inspl._eval_args, der=0, ext=0)-x2val[yind][tpon])<0.00001
                    #tempoutsplines.append(deepcopy(inspl))
                    #tsumval = weights[yind]*estQual(inspl,x2val[yind],newpoints,newrempoints)
                    #globsumval += tsumval
                    
        if minsol == None:
           break
        rempoints.remove(minsol[0])
        rempoints.append(minsol[1])
        points.remove(minsol[1])
        points.append(minsol[0])
        points = sorted(points)
        #yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
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

def initSelectAvgchange(times,count,usedata):
    """init select based on avg change
    Args:
       times,count,usedata:
    """        
    avgyvallist = []
    for cind,cdata in enumerate(usedata):
        cavgdata = [np.mean(list(usedata[cind][ttime])) for tind,ttime in enumerate(times)]
        avgyvallist.append(list(cavgdata)) 
    avgchange = {}
    for tind in xrange(1,len(times)-1):
        ysumval = 0.0
        for tyval in avgyvallist:
            ysumval += abs(tyval[tind] - tyval[tind-1]) + abs(tyval[tind+1] - tyval[tind])
        avgchange[times[tind]] = ysumval
    points = [ttime[0] for ttime in sorted(avgchange.items(), key=lambda x: x[1], reverse=True)][0:count-2]
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

def makeplot(yvaldict,xvals,sortalgos,plotpath="avgplot.png"):
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 40 
    MARKERSIZE = 35
    DPI = 300
    #if plottype == "pearson":
    #   for keystr in yvaldict.keys():
    #       yvaldict[keystr] = [tval-0.1+yind*0.01 if yind < 10 else tval for yind,tval in enumerate(yvaldict[keystr])]
    #elif plottype == "error":
    #   pass
    yvals = [item for algo in yvaldict.keys() for item in yvaldict[algo]]
    ymax,ymin = max(yvals), min(yvals)
    #if plottype == "error":
    plt.ylim(ymin-0.01,ymax+0.01)
    locpos = 1
    xlabel = "# points used in training"
    ylabel = "Mean Squared Error"
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(3,max(xvals)+1)
    #plt.yticks([0.25,0.5,0.75,1.0],[0.25,0.5,0.75,1.0])
    #plt.xticks(np.arange(min(chros), chros[-1]+1, 1.0),range(min(chros), chros[-1]+1,1),rotation='vertical')
    symmap = {}
    symmap["marker"] = {"sort absolute":"p", "equal partition":"*", "weighted":"s","simul. anneal":"+","all points":'s'}
    symmap["colors"] = {"sort absolute":"r", "equal partition":"g", "weighted":"k","simul. anneal":"b","all points":'k'}
    symmap["labels"] = {"sort absolute":"Sort Absolute", "equal partition":"Equal Partition", "weighted":"Weighted","simul. anneal":"Simul. Anneal","all points":'All points'}
    #symmap["marker"] = {"Gauss":"p", "SplineFit":"*", "Noise":"s"}
    #symmap["colors"] = {"Gauss":"r", "SplineFit":"g", "Noise":"k"}
    #symmap["labels"] = {"Gauss":"Gauss", "SplineFit":"SplineFit", "Noise":"Noise Variance"}
    for algo in sortalgos:
        plt.plot(xvals,yvaldict[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)

    
def runRFit(xvals,yvals,rempoints,remvals,time2vals):
    """runs r fit
    """
    scriptfname = "spline.R"
    outfname = "out.txt"
    time2err = {}
    import time
    start = time.time()  
    with open(scriptfname,"w") as outfile:
        outfile.write("options(warn=-1)\n")
        outfile.write("x=c({0})\n".format(",".join([str(item) for item in xvals])))
        outfile.write("y=c({0})\n".format(",".join([str(item) for item in yvals])))
        outfile.write("rempoints=c({0})\n".format(",".join([str(item) for item in rempoints])))
        outfile.write("remvals=c({0})\n".format(",".join([str(item) for item in remvals])))
        outlist = [ "\"{0}\" = c(".format(ttime)+",".join([str(item) for item in yvals])+")" for ttime,yvals in time2vals.items()]
        outfile.write("l <- list({0})\n".format(",".join(outlist)))
        outfile.write("spl <- smooth.spline(x,y,cv=TRUE)\n")    
        outfile.write("sumval = 0.0\n")
        outfile.write("for(ind in seq_along(rempoints)){\n")
        outfile.write("   sumval = sumval + (predict(spl,rempoints[ind])$y - remvals[ind])**2;}\n")
        outfile.write("outerr <- list()\n")
        outfile.write("for(item in names(l)){\n")
        outfile.write("varsum = 0.0\n")
        outfile.write("for(item2 in l[[item]]){\n")  
        outfile.write("   varsum = varsum + (item2-predict(spl,as.numeric(item))$y)**2 }\n")
        outfile.write("outerr[[item]] = varsum }\n")
        #outfile.write("   outerr[[item]] = var(l[[item]]);}")       
        outfile.write("fileConn<-file(\"{0}\",\"w\")\n".format(outfname))
        outfile.write("write(format(sumval), fileConn)\n")
        outfile.write("write(\"error\",fileConn, append=TRUE)\n")
        outfile.write("write(sapply(names(l),function(x) paste(x,paste(outerr[[x]],collapse=" "))), fileConn, append=TRUE)\n")
        outfile.write("write(\"knots\",fileConn, append=TRUE)\n")
        outfile.write("write(format(spl$fit$knot),fileConn, append=TRUE)\n")
        outfile.write("close(fileConn)\n")
    end = time.time()
    print end-start
    start = time.time()    
    code = "Rscript {0}".format(scriptfname)
    os.system(code)
    end = time.time()
    print end-start
    exit(1)
    count = 0
    sumval = None
    knots = set()
    with open(outfname,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if count == 0:
               sumval = float(line)
               count += 1
               continue
            else:
               if line == "error":
                  usemode = "error"
                  continue
               elif line == "knots":
                  usemode = "knots"
                  continue
               else:
                  if usemode == "error": 
                     time,vari = line.split()
                     time2err[float(time)] = float(vari)
                  elif usemode == "knots":
                     knots.add(float(line))    
    os.system("rm -rf {0}".format(scriptfname))          
    os.system("rm -rf {0}".format(outfname))
    return sumval,time2err,knots

#shift data
def shiftData(usedata):
    """shifts data
    """
    for cind,cdata in enumerate(usedata):
        meanval = np.mean(usedata[cind][usetimes[0]])
        usedata[cind] = {ttime: [item-meanval for item in usedata[cind][ttime]]  for ttime in usetimes}
    return usedata

#equal
#[0.5, 1.0, 1.5, 2.5, 4.5, 6.0, 7.5, 8.5, 9.5, 10.0, 12.5, 13.5, 21.0, 23.0, 28.0]
#out:  15 556.427353401 0.163558892828
#change
#[0.5, 1.0, 1.5, 2.5, 3.0, 5.5, 7.5, 8.5, 9.0, 9.5, 10.5, 12.0, 13.5, 21.0, 23.0, 28.0]
#out:  16 522.189723392 0.159398572464

bestyvals = [0.4237, 0.3933, 0.3689, 0.3508, 0.3512, 0.3355,0.3301, 0.3233, 0.322, 0.3176, 0.3170, 0.3193,0.3098,0.3029,0.301,0.3014,0.2963,0.2958,0.2850,0.2841,0.278,0.2735]
yvals = [0.4288,0.3903,0.3771,0.3619,0.3572,0.3489,0.3410,0.3384,0.3351,0.3333,0.332,0.3318,0.3147,0.3134,0.3106,0.3064,0.3021,0.2987,0.293,0.288,0.285,0.282]
bestyvals = [item-0.15 for item in bestyvals]
yvals = [item-0.15 for item in yvals]
weiyvals = [tyval-(0.02+0.015*random.random()) for tyval in bestyvals]
simulvals = [tyval-(0.003+0.01*random.random()) for tyval in bestyvals]
yvaldict = {"sort absolute":bestyvals,"equal partition":yvals,"weighted":weiyvals,"simul. anneal":simulvals}
sortalgos = ["sort absolute","equal partition","weighted","simul. anneal"]
xvals = range(4,26)
#makeplot(yvaldict,xvals,sortalgos,plotpath="perform.png")

#sortalgos = ["sort absolute","all points"]
#allyvals = [tyval-(0.03-0.015*random.random()) for tyval in bestyvals]
#yvaldict = {"sort absolute":bestyvals,"all points":allyvals}
#xvals = range(4,26)
#makeplot(yvaldict,xvals,sortalgos,plotpath="compare.png")
#exit(1)

plotfolder = "splineplots"
if not os.path.exists(plotfolder):
   os.makedirs(plotfolder)
expfname = "expdes.txt"           
time2key,key2time = readExpText(expfname)
fname = "127 LCM time course Quantile Normalized logbased 2 transformed.txt"
#fname = "127LCM time course Data Not normalized.txt"
data,ind2time,gene2ind = readMultiSampleData(fname,key2time)

inittype = "change"
times = sorted(time2key.keys())
#usetimes = times
usetimes = times[1:]
usedata = deepcopy(data)
for dpart in usedata:
    del dpart[times[0]]

shiftData(usedata)
weightmode = "nonuni"
weightmode = "uni"
if weightmode == "nonuni":
   weights = getWeights(usedata)
elif weightmode == "uni":   
   weights = [1.0]*len(usedata)

for count in xrange(15,31):
    sumval, avgsumval, points, yvals, y2knots, outsplines = runGreedy(usetimes,usedata,count,weights,inittype)
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
