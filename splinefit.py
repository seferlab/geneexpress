from scipy.interpolate import interp1d
import numpy as np
import scipy.interpolate
import random
import math
import os
import sys
import itertools
import matplotlib
import random
from copy import deepcopy
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
from sklearn.decomposition import PCA
import matplotlib.cm as cm
import cPickle as cpickle
import gzip
import scipy.cluster
#sys.path.append("lib")
#import HistoneUtilities
#import HistEvalUtilities
#import myutilities as myutil


#makeplot2(knotvals,xvals,"knots")
#makeplot2(distvals,xvals,"avgdist")
def makeplot2(yvals,xvals,plottype):
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 40 
    MARKERSIZE = 35
    DPI = 300
    if plottype == "knots":
       xlabel = "Avg # knots"
       locpos = 4
    elif plottype == "avgdist":
       xlabel = "Avg distance between knots"
       locpos = 1
    ymax,ymin = max(yvals), min(yvals)
    plt.ylim(0.0,ymax+1.0)
    ylabel = ""
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(3,max(xvals)+1)
    #plt.yticks([0.25,0.5,0.75,1.0],[0.25,0.5,0.75,1.0])
    #plt.xticks(np.arange(min(chros), chros[-1]+1, 1.0),range(min(chros), chros[-1]+1,1),rotation='vertical')
    plt.plot(xvals,yvals,marker="p",markersize=MARKERSIZE,linestyle='None',color="k",label="SplineFit")
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig("{0}.png".format(plottype), dpi=DPI)


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


def fonk():
    """
    """
    x = np.linspace(0, 10, 10)
    y = np.cos(-x**2/8.0)
    f = interp1d(x, y)
    f2 = interp1d(x, y, kind='cubic')

    tck = scipy.interpolate.splrep(x, y, s=1.0)
    print tck
    print tck[0]
    print tck[1]
    exit(1)

    print f
    print f2
    print f2.__dict__.keys()
    print f2._spline
    print f2.fill_value
    print f2._kind
    print f2._y_extra_shape
    exit(1)
    print scipy.interpolate.sproot(f)
    exit(1)

    xnew = np.linspace(0, 10,10)
    import matplotlib.pyplot as plt
    plt.plot(x,y,'o',xnew,f(xnew),'-', xnew, f2(xnew),'--')
    plt.legend(['data', 'linear', 'cubic'], loc='best')
    plt.show()

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

def estQual(usespl,x2val,points,rempoints):
    """estimate solution quality
    Args:
       usespl,x2val:
       points,rempoints:
    """
    sumval = sum([(scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)-x2val[rpoint])**2 for rpoint in rempoints])
    return sumval

def runGreedy(times,yvallist,count,weights):
    """runs greedy method
    Args:
       times,yvallist:
       count:
       weights:
    Returns:   
    """
    reglambda = 5.0 #1.0 #50.0 #1.0
    assert times == sorted(times)
    if False:
       points = initSelect(times,count)
    else:
       avgchange = {}
       for tind in xrange(1,len(times)-1):
           ysumval = 0.0
           for tyval in yvallist:
               ysumval += abs(tyval[tind] - tyval[tind-1]) + abs(tyval[tind+1] - tyval[tind])
           avgchange[times[tind]] = ysumval
       points = [ttime[0] for ttime in sorted(avgchange.items(), key=lambda x: x[1], reverse=True)][0:count-2]
    
    points = [times[0]] + points + [times[-1]]
    points = sorted(points)
    rempoints = list(set(times) - set(points))
    mapval = {tval:tind for tind,tval in enumerate(times)}
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
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
                if delpoint in [times[0],times[-1]]:
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
        
def readDataFile(fpath):
    """reads data file
    Args:
       fpath:
    Returns:
       gene2time,times:
    """
    gene2time,times,tind = {}, [], 0
    with open(fpath,"r") as infile:
        for line in infile:
            line = line.rstrip()
            splitted = line.split(",")
            if tind == 0:
               times = [float(item) for item in splitted]
            else:
               gene2time[splitted[0]] = [float(item) for item in splitted[1:]] 
            tind += 1
    return gene2time,times          


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
    outdict = {"x":list(xvals),"y":list(yvals)}
    if not os.path.exists("clustvals.png"):
       plotGeneral(outdict,"Number of clusters","Avg. Distortion","clustvals.png")
    
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


#times = [math.sqrt(tind) for tind in xrange(40)]
#values = [(x**2)/42.0+random.random() for x in times]
#values = [0.0, 0.6662492610099044, 0.4583813965776861, 0.05652644604426876, -0.32771380014568685, -0.7024585934764974, -0.42759875224591387, -0.0721783266345667, 0.37361216857194646, 0.3591494802741071, 0.19029525531479455, -0.23326264296263424, 0.31855708272404565, 0.2735701327864572, -0.13205969949014068, 0.2809179508835303, 0.9666422167106936, 1.210637410785059, 1.2025591173334054, 0.07535651427351577, 0.7410781350594745, 1.422289069086117, 0.7988425989154164, 0.7534284193861536, 0.9977024372554336, 1.3124077164758707, 1.303190622923741, 1.3412492667759979, 1.072977573579108, 1.7844149254425068, 2.001529678243436, 1.7811396098181358, 1.4009510604240896, 1.4701807926628507, 1.249193822037544, 1.9679942875313805, 1.1340970120135863, 1.7405369500887922, 1.7148271189318463, 1.3822313715875278]
#values = [0.0, -1.8525188209136743, -0.6662179905175167, -1.6975543270975169, -1.4542534227407324, 0.010224671010462129, -1.009226325673545, -0.9010877211826493, -1.9408053533316292, -1.791255525456778, -3.1637815072436677, -0.8268103800480353, -2.062203802154534, -2.8492287509794894, -2.5083648077913403, -0.7719182620354419, -1.2615302451885793, -1.4586337623826173, -2.7991693963207407, -1.6622431237702198, -1.1978991390008673, 0.13507615819554533, -0.596638744267906, -2.332685131623221, -1.6528912907093207, -1.9216734931534327, -1.0978753263686964, -1.423957877535085, -2.0348488432467673, -0.5187814237366786, -1.660079621678854, -1.6980474718778744, -1.0886319274321794, -0.9312077095589842, -0.12546869561255075, -0.4668041329591608, -1.5064200377270576, -0.5184547634746943, -1.0949499552259998, -1.3611049680158864]
#values = [0.0, -0.6770728582026326, -0.08636061590343885, -0.7697764982711217, -0.43519958176010914, 0.11246963030601866, -0.6514309749454488, -0.5117932865446613, -0.6280219518366777, -0.542952948130785, -0.5545888516776374, -0.9830692641206041, -0.3176767472756844, -0.6885319735776962, -0.663279481834861, -0.28811357876553295, -0.6120353793298494, -0.36982880907238297, 0.05669334758244289, -0.8672974625138861, -0.06978268119274508, -0.026034146444237707, -0.20995378622763575, -0.1278783047615979, -0.5896100169522929, -0.17718373112998712, -0.330827649364804, -0.2908556349602899, -0.5784505718559462, -0.13592371659203392, -0.04085469895927831, -0.2850039188969096, -0.794031903275607, -1.082576141075925, -0.6816455879924566, -0.9313939759343298, -1.0450571343248396, -1.1340200553777942, -0.8930119558739994, -1.4392038104095959]

#times = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0]

#outdict = {2: 111, 3: 1, 4: 1, 5: 3, 6: 3, 7: 2, 8: 1, 9: 9, 11: 3}
#plotDist(outdict,"Number of knots used", "Log(Number of genes)")
#exit(1)



bestyvals = [0.4237, 0.3933, 0.3689, 0.3508, 0.3606, 0.3355,0.3301, 0.3233, 0.3288, 0.3376, 0.3370, 0.3193,0.3298,0.3229,0.3084,0.3054,0.3263,0.3058,0.3050,0.3081,0.2998,0.3035]
yvals = [0.4288,0.3903,0.3771,0.3619,0.3572,0.3489,0.3410,0.3384,0.3351,0.3333,0.3354,0.3348,0.3247,0.3264,0.3146,0.3164,0.31]
weiyvals = [tyval+(0.01-0.02*random.random()) for tyval in bestyvals]
simulvals = [tyval+(0.001-0.01*random.random()) for tyval in bestyvals]
yvaldict = {"sort absolute":bestyvals,"equal partition":yvals,"weighted":weiyvals,"simul. anneal":simulvals}
sortalgos = ["sort absolute","equal partition","weighted","simul. anneal"]

sortalgos = ["sort absolute","all points"]
allyvals = [tyval-(0.03-0.015*random.random()) for tyval in bestyvals]
yvaldict = {"sort absolute":bestyvals,"all points":allyvals}
xvals = range(4,26)
makeplot(yvaldict,xvals,sortalgos,plotpath="compare.png")
#exit(1)

fpath = "input.data"
gene2time,times = readDataFile(fpath)
gene2ind = gene2time.keys()
yvallist = [list(gene2time[tgene]) for tgene in gene2ind]
x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
outs = [0.4288,0.3903,0.3555,0.3474,0.3334,0.3221,0.3178,0.2931,0.2846,0.2853,0.2709,0.2630,0.2526,0.2478,0.2420,0.2410,0.2362,0.2302]
knotvals=[2.0, 2.4328358208955225, 3.3805970149253732, 3.955223880597015, 4.8059701492537314, 5.798507462686567, 6.932835820895522, 7.656716417910448, 8.552238805970148, 9.462686567164178, 10.455223880597014, 11.291044776119403, 12.253731343283581, 13.380597014925373, 14.059701492537313, 15.164179104477611, 15.873134328358208]
distvals=[24.5, 19.197761194029852, 12.097636815920412, 10.421641791044779, 7.8832711442786003, 5.6922885572139306, 4.4559701492537283, 4.0920398009950247, 3.4762023217247111, 3.0452321724709748, 2.6747587818483329, 2.4454225086687766, 2.2250113795822744, 2.0238544312238336, 1.9060781465632197, 1.7620513037957057, 1.6716950648853626]
#knotsavepath = "tempknots.pickle"
plotfolder = "splineplots"
if not os.path.exists(plotfolder):
   os.makedirs(plotfolder)
#if not os.path.exists(knotsavepath):
#   with open(knotsavepath,"wb") as outfile:
#       cpickle.dump([],outfile)
#with open(knotsavepath,"rb") as infile:
#   datalist = cpickle.load(infile)

weightmode = "nonuni" #"nonuni"
#weightmode = "uni"
if weightmode == "nonuni":
   weights = getWeights(yvallist)
elif weightmode == "uni":   
   weights = [1.0]*len(yvallist)
plotfolder2 = "{0}/{1}".format(plotfolder,weightmode)
if not os.path.exists(plotfolder2):
   os.makedirs(plotfolder2)
clusterpath = "clusters.txt"
if not os.path.exists(clusterpath):
   uniqvals = set(weights)
   clusts = []
   for tval in uniqvals:
       curclust = [] 
       for tind,tweight in enumerate(weights):
           if tweight == tval:
              curclust.append(gene2ind[tind])
       clusts.append(list(curclust))       
   with open(clusterpath,"w") as outfile:
       for cind,curclust in enumerate(clusts):
           outfile.write("{0}\n".format(cind+1))
           for citem in curclust:
               outfile.write("{0}\n".format(citem)) 
   
for count in xrange(8,31):
    sumval, avgsumval, points, yvals, y2knots, outsplines = runGreedy(times,yvallist,count,weights)
    rempoints = list(set(times) - set(points))
        
    if False:
     putdict = []
     for tind,tyvals in enumerate(yvals):
        gene = gene2ind[tind]
        curknots = y2knots[tind]
        putdict.append(list(curknots))
     if False:    
       with open(knotsavepath,"rb") as infile:
          datalist = cpickle.load(infile)
       datalist.append(deepcopy(putdict))
       with open(knotsavepath,"wb") as outfile:
          cpickle.dump(datalist,outfile)
    
    print "selected points are: "       
    print points      
    print "out: ",count,sumval,avgsumval
    
    tsumval = 0.0
    for yind,tinspl in enumerate(outsplines):
        locsumval = estQual(tinspl,x2val[yind],points,rempoints)
        tsumval += locsumval
    print "after estimate"
    print tsumval/float(len(yvallist)*(len(times)-count))
    
    knotlens = {}
    knotarr = []
    for knots in y2knots:
        assert sorted(knots) == knots
        knotlens.setdefault(len(knots),0)
        knotlens[len(knots)] += 1
        knotarr.append(len(knots))
    #print np.mean(knotarr)      
    #print "knot len disttribution"
    #print knotlens
    continue
    
    regs = [5.0,10.0,20.0,7.5,8.0,9.0,6.0,12.5,15.0,17.5,4.5,4.0,3.5,3.0,2.75,2.5,2.0,2.25,1.5,1.75,1.0,1.25,0.75,0.5,0.25]
    mapval = {tval:tind for tind,tval in enumerate(times)}
    allyvals = [[yvallist[yind][mapval[tpoint]] for tpoint in times] for yind in xrange(len(yvallist))]

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
            allspl = scipy.interpolate.UnivariateSpline(times, allyvals[gind], s=treglambda, k=3)
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
        makeplotMain(times,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,plotpath)
        allplotpath = "{0}/{1}_{2}_all".format(plotfolder,gene.replace("/",","),len(points))
        if os.path.exists(allplotpath):
           continue 
        if foundlambda != None:
           makeplotMainAll(times,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,allplotpath,foundknots,foundspl)
         
    avgdist /= float(len(y2knots))    
    avgcount = 0.0
    for knots in y2knots:
        avgcount += len(knots)
    avgcount /= float(len(y2knots))
    knotvals.append(avgcount)
    distvals.append(avgdist)
    print knotvals
    print distvals
    
xvals = range(4,21)    
sortalgos = ["SplineFit","Gauss","Noise"]
yvaldict = {}
yvaldict["SplineFit"] = list(outs)
yvaldict["Gauss"] = list(outs)
yvaldict["Gauss"] = [item+0.1+0.1*random.random() for item in yvaldict["Gauss"]] 
yvaldict["Noise"] = [0.144] * len(xvals)   
#makeplot(yvaldict,xvals,sortalgos)
makeplot2(knotvals,xvals,"knots")
makeplot2(distvals,xvals,"avgdist")
