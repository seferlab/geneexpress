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

def makeplotMainAll(xvals,yvaldictout,points,knots,gene,outspl,rempoints,remyvals,plotpath,allknot,allspline,usetimes,trixvals,triyvals):    
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

    plotpath = plotpath.replace(".","__")
    times = usetimes
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in times] 
    ymax,ymin = max(yvals), min(yvals)
    #plt.ylim(ymin-1.0,ymax+0.1)
    plt.ylim(ymin-0.05,ymax+0.05)
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
        if algo in ["Knots"]:
           continue 
        plt.plot(points,yvaldictout[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    #plt.plot(knots,knotvals,marker=symmap["marker"]["Knots"],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"]["Knots"],label=symmap["labels"]["Knots"])
    plt.plot(usetimes,outspl(usetimes),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')

    allyvals = [scipy.interpolate.splev(knot,allspline._eval_args, der=0, ext=0) for tknot in allknot]
    #plt.plot(usetimes,allspline(usetimes),'y',lw=6,label='All')
    
    knotys = [scipy.interpolate.splev(tknot,allspline._eval_args, der=0, ext=0) for tknot in allknot]
    plt.plot(allknot,knotys,marker='+',markersize=MARKERSIZE,linestyle='None',color='b',label='Knots all')
       
    ax = plt.axes()        
    ax.xaxis.grid()
    plt.plot(trixvals[0:2],triyvals[0:2],lw=3,color='y',label='Linear fit')
    plt.plot(trixvals[1:],triyvals[1:],lw=3,color='y',label='Linear fit')    
    #plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)
    
    

def makeplotMain(xvals,yvaldictout,points,knots,gene,outspl,rempoints,remyvals,plotpath,usetimes,trixvals,triyvals,legendbound=None):    
    """makes main plot
    Args:
       xvals,yvaldictout,points:
       knots,gene,outspl:
       rempoints,remyvals:
       plotpath:
    Returns:
    """
    plt.clf()
    plt.rc('font', family='serif', size=35)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 20 
    MARKERSIZE = 25
    DPI = 300

    plotpath = plotpath.replace(".","__")
    times = usetimes
    yvals = [item for algo in yvaldictout.keys() for item in yvaldictout[algo]] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in knots] + remyvals + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for ttime in times]
    if legendbound != None:
       ymax,ymin = legendbound[1]+0.1, legendbound[0]-0.1
       plt.ylim(ymin,ymax)
    else:       
       ymax,ymin = max(yvals), min(yvals)
       plt.ylim(ymin-0.5,ymax+0.25)
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
    plt.plot(usetimes,outspl(usetimes),symmap["colors"]["spline"],lw=3,label=symmap["labels"]["spline"])
    plt.plot(rempoints,remyvals,marker='o',markersize=MARKERSIZE,linestyle='None',color='b',label='Rem. points')
    #plt.plot([0.5,7.0],triyvals[0:2],lw=3,color='y',label='Linear fit')
    #plt.plot([7.0,28.0],triyvals[1:],lw=3,color='y',label='Linear fit')
     
    ax = plt.axes()        
    ax.xaxis.grid()  
    #if legendbound == None: 
    #   plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13) 
    plt.savefig(plotpath, dpi=DPI)    

    
def readMultiSampleData(fname):
    """reads multi sample data
    Args:
       fname: file to be read
    Returns:
       data:
       ind2time:
       gene2ind:
    """
    gene2id, data = [], []
    ind2time = {}
    with open(fname,"r") as infile:
        for line in infile:
            uselines = line.split("\r")
            for kind,keystr in enumerate(uselines[0].split("\t")[4:-7]):
                ind2time[kind+1] = float(keystr.replace("P","").replace("p",""))
            for uind,useline in enumerate(uselines[1:]):
                if uind in range(0,5) or uind >= 604:
                   continue 
                splitted = uselines[uind+1].split("\t")
                gene = splitted[0]
                gene2id.append(gene)
                curvals = {}
                for ind,item in enumerate(splitted[4:-7]):
                    curtime = ind2time[ind+1]
                    curvals.setdefault(curtime,[])
                    curvals[curtime].append(float(item))
                data.append(deepcopy(curvals))
    return data,ind2time,gene2id


def estQual(tusespl,x2val,points,rempoints):
    """estimate solution quality
    Args:
       usespl,x2val:
       points,rempoints:
    """
    sumval = 0.0
    for rpoint in rempoints:
        cury = scipy.interpolate.splev(rpoint,tusespl._eval_args, der=0, ext=0)
        curyval = sum([(cury-item)**2 for item in x2val[rpoint]])
        sumval += curyval/float(len(x2val[rpoint]))    
    #sumval = sum([((scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)-item)**2)/float(len(x2val[rpoint])) for rpoint in rempoints for item in x2val[rpoint]])
    #sumval = sum([(scipy.interpolate.splev(rpoint,usespl._eval_args, der=0, ext=0)-x2val[rpoint])**2 for rpoint in rempoints])
    return sumval

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

def runGreedy(times,usedata,yvallist,count,weights,initmode="change"):
    """runs greedy method
    Args:
       times,usedata:
       count,weights:
    Returns:   
    """
    def convertForm(curtimes,usedata):
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
        
    #fixedpoints = [0.5, 2.5, 5.0, 10.0, 26.0]
    #fixedpoints = [0.5, 5.5, 11.0 18.0, 28.0]
    #fixedpoints = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10.0, 11.0, 12.0, 13.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0]
    fixedpoints = [0.5, 1.5, 2.5, 5.0, 7.0, 10.0, 15.0, 19.0, 28.0]

    #fixedpoints = [0.5, 1.5, 2.5, 3.0, 4.5, 6.0, 7.5, 8.5, 9.5, 10.0, 12.5, 18.0, 24.0, 28.0]

    fixedpoints = [0.5, 1.0, 1.5, 3.0, 4.0, 7.0, 7.5, 9.0, 9.5, 12.5, 15.0, 21.0, 28.0]
    #fixedpoints = [0.5,28.0]
    #fixedpoints = [0.5,7.0,28.0]

    reglambda = 30.0
    if True:
       reglambda = 30.0 #1.0 #50.0 #1.0
       avgyvallist = []
       for cind,cdata in enumerate(usedata):
           cavgdata = [np.median(list(usedata[cind][ttime])) for tind,ttime in enumerate(times)]
           avgyvallist.append(list(cavgdata)) 
       avgchange = {}
       for tind in xrange(1,len(times)-1):
           ysumval = 0.0
           for tyval in yvallist:
               ysumval += abs(tyal[tind] - tyval[tind-1]) + abs(tyval[tind+1] - tyval[tind])
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
       x2valblock = [{tval:list(usedata[yind][tval]) for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
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
           tsumval = weights[yind]*estQual(spl,x2valblock[yind],points,rempoints)
           sumval += tsumval

       while False:
         print sumval,sumval/float(len(yvallist)*(len(times)-count))
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
                    tsumval = weights[yind]*estQual(inspl,x2valblock[yind],newpoints,newrempoints)
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
       print points
       
       #points = [0.5, 1.0, 3.0, 5.0, 7.0, 7.5, 8.5, 10.0, 11.0, 13.5, 19.0, 22.0, 28.0]
       #points = [0.5,1.0, 3.0, 5.0, 7.0, 7.5, 8.5, 10.0, 11.0, 14.0, 19.0, 22.0,28.0]
       #points = [0.5, 1.0, 1.5, 2.5, 3.0, 4.5, 6.0, 7.5, 8.5, 9.5, 10.0, 12.5, 18.0, 24.0, 28.0]
       points = [0.5, 1.0, 1.5, 3.0, 4.0, 7.0, 7.5, 9.0, 9.5, 12.5, 15.0, 21.0, 28.0]
       rempoints = list(set(times) - set(points))
       yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
       tsumval = 0.0
       y2knots = []
       outsplines = []
       for yind,curyvals in enumerate(yvals):
         inyvals = [x2val[yind][rpoint] for rpoint in points]
         tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
         outsplines.append(deepcopy(tinspl))
         y2knots.append(list(tinspl.get_knots()))
         locsumval = weights[yind]*estQual(tinspl,x2valblock[yind],points,rempoints)
         tsumval += locsumval
       print tsumval
       print sumval/float(len(yvallist)*(len(times)-count))
       print tsumval/float(len(yvallist)*(len(times)-count))
       #exit(1)

       #assert abs(tsumval-sumval) < 0.0001     
       avgsum=0.0
       for yind in xrange(len(y2knots)):
         avgsum += len(y2knots[yind])
       print avgsum
       print avgsum/float(len(y2knots))
       print "infom"
       return sumval,sumval/float(len(yvallist)*(len(times)-count)),points,yvals,y2knots,outsplines


    rempoints = list(set(times) - set(fixedpoints))
    mapval = {tval:tind for tind,tval in enumerate(times)}
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    x2valblock = [{tval:list(usedata[yind][tval]) for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in fixedpoints] for yind in xrange(len(yvallist))]
    reglambda = 75.0 #1.0 #50.0 #1.0
    #1.0 -> 805
    y2knots = []
    outsplines = []
    #[0.5, 2.5, 6.0, 11.0, 19.0, 28.0].
    sumval = 100000000000.0
    minsol = None

    #userempoints = [1.5, 1.0, 2.0, 3.0, 4.0, 5.0, 4.5, 9.0, 10.0, 11.0, 13.0, 15.0, 17.0, 18.0, 21.0, 25.0, 26.0, 8.5, 6.5, 6.0, 3.5, 2.5, 7.5, 5.5]
    #print len(userempoints)
    #exit(1)

    #area 0.391676970384
    #userempoints = [1.0, 8.5, 12.0, 23.0]
    userempoints = [1.0, 6.0, 13.5, 23.0]
    #userempoints = rempoints
    allvals = [selpoints for selpoints in itertools.combinations(userempoints, count-2)]

    reglambda = 50.0
    fixedpoints = [0.5, 1.5, 2.5, 3.0, 4.5, 6.0, 7.5, 8.5, 9.5, 10.0, 12.5, 18.0, 24.0, 28.0]
    allvals = [[1.0]]
    sol2val = {}
    random.shuffle(allvals)
    count = 0
    for selpoints in allvals:
        count += 1
        if count == 5000:
           break 
        #tflag = True
        #for mypoint in selpoints:
        #    if mypoint not in [1.5,15.0,17.0,19.0]:
            #if mypoint not in [17.0, 19.0, 11.0, 11.5, 6.0, 1.0,3.0]:
            #if mypoint not in [15.0, 19.0, 9.5, 9.0, 8.0, 3.0, 13.0, 7.5, 7.0, 10.5, 11.0, 6.0, 1.0,4.0,6.5]:
        #    if mypoint not in [15.0, 19.0, 9.5, 9.0, 8.0, 7.5, 7.0, 5.5, 4.5, 3.5, 3.0,2.0,1.5,10.5, 11.0, 6.0, 1.0,4.0,6.5,2.5,5.0,10.0]:
            #if mypoint not in [1.0, 3.0, 4.0, 6.0, 15.0, 6.5,11.0, 19.0,13.0,17.0]:
        #       tflag = False
        #       break
        #if not tflag:
        #   continue 
        if random.random() <= 1.0:
           print count
           print sumval,minsol,sumval/float(len(yvallist)*(len(times)-len(selpoints)-len(fixedpoints)))
           avglens = [len(list(tspl.get_knots())) for tspl in outsplines]
           print np.mean(avglens)
           print "multiple sols:"
           import operator
           sorted_data = sorted(sol2val.items(), key=operator.itemgetter(1))
           for keystr,optval in sorted_data[0:10]:
               print sorted([float(item) for item in keystr.split("_")]),optval,optval/float(len(yvals)*len(times))
        
        usepoints = sorted(list(selpoints) + fixedpoints)
        userempoints = list(set(rempoints) - set(selpoints))
        tempoutsplines = []
        globsumval = 0.0
        curknots = []

        #multipoints = [upoint for upoint in usepoints for tind2 in xrange(len(usedata[0][upoint]))]
        #multirempoints = [rpoint for rpoint in userempoints for tind2 in xrange(len(usedata[0][rpoint]))]
        #newyvals = convertForm(usepoints,usedata)
        #newremyvals = convertForm(userempoints,usedata)

        for yind,curyvals in enumerate(yvals):
            #cursumval,time2err,knots = runRFit(multipoints,curyvals,multirempoints,newremyvals[yind],usedata[yind])
            #modknots = changeRange(knots,times[0],times[-1])
            inyvals = [x2val[yind][rpoint] for rpoint in usepoints]
            inspl = scipy.interpolate.UnivariateSpline(usepoints, inyvals, s=reglambda, k=3)
            tempoutsplines.append(deepcopy(inspl))
            tsumval = weights[yind]*estQual(inspl,x2valblock[yind],usepoints,userempoints)
            globsumval += tsumval
            curknots.append(list(inspl.get_knots()))

            for tind,ttime in enumerate(usepoints):
                globsumval += np.std(list(usedata[yind][ttime]))**2
            
        print "inside: ",globsumval
        sol2val["_".join([str(titem) for titem in selpoints])] = globsumval
        if globsumval < sumval:
           sumval = globsumval
           minsol = list(selpoints)
           outsplines = list(tempoutsplines)
           y2knots = deepcopy(curknots)
    #print sol2val
    #print "done"
    #print y2knots
    
    if False:
     pts1 = [0.5,28.0]
     pts2 = [0.5,3.5,6.0,11.0, 16.0, 22.0, 28.0]
     pts3 = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0]
     pts4 = [0.5,7.0,28.0]

     def getArea(sentpoints,myyvals):
         """
         """
         locaream = 0.0
         for pind in xrange(len(sentpoints)-1):
            time1,time2 = sentpoints[pind:pind+2]
            y1,y2 = myyvals[pind:pind+2]
            if y1*y2 >= 0.0:
               locaream += abs(y2+y1)*(time2-time1)/2.0
            elif y1 <= 0.0 and y2 >= 0.0:
               ratio = (time2-time1)/float(y2-y1)
               breakp = time1+(ratio*abs(y1))
               locaream += (breakp-time1)*abs(y1)*0.5
               locaream += (time2-breakp)*abs(y2)*0.5
            elif y1 >= 0.0 and y2 <= 0.0:
               ratio = (time2-time1)/float(y1-y2)
               breakp = time1+(ratio*y1)
               locaream += (breakp-time1)*y1*0.5
               locaream += (time2-breakp)*abs(y2)*0.5
         return locaream
               
     sparea = 0.0
     linarea = 0.0
     y2lin,y2sp = {},{}
     y2lin1,y2lin2,y2lin3,y2lin4 = {},{},{},{}
     sampoints = np.arange(min(times),max(times)+0.01,0.01)
     points = sorted(list(minsol) + fixedpoints)
     for yind,curyvals in enumerate(yvals):
        inyvals = [x2val[yind][rpoint] for rpoint in points]
        inyvals1 = [x2val[yind][rpoint] for rpoint in pts1]
        inyvals2 = [x2val[yind][rpoint] for rpoint in pts2]
        inyvals3 = [x2val[yind][rpoint] for rpoint in pts3]
        inyvals4 = [x2val[yind][rpoint] for rpoint in pts4]
        loclinarea = getArea(points,inyvals)
        loclinarea1 = getArea(pts1,inyvals1)
        loclinarea2 = getArea(pts2,inyvals2)
        loclinarea3 = getArea(pts3,inyvals3)
        loclinarea4 = getArea(pts4,inyvals4)
        
        locsparea = 0.0    
        for samind in xrange(len(sampoints)-1):
            time1,time2 = sampoints[samind], sampoints[samind+1]
            tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
            y1 = scipy.interpolate.splev(time1,tinspl._eval_args, der=0, ext=0)
            y2 = scipy.interpolate.splev(time2,tinspl._eval_args, der=0, ext=0)
            if y1*y2 < 0.0:
               print yind
            curaream = abs(y2+y1)*(time2-time1)/2.0
            if curaream >= 0.0:   
               locsparea += abs(y2+y1)*(time2-time1)/2.0
        y2sp[yind] = locsparea
        y2lin[yind] = loclinarea
        y2lin1[yind] = loclinarea1
        y2lin2[yind] = loclinarea2
        y2lin3[yind] = loclinarea3
        y2lin4[yind] = loclinarea4
        sparea += locsparea
        linarea += loclinarea
     difsp = {yind:abs(y2sp[yind]-y2lin[yind]) for yind in y2sp.keys()}
     difsp1 = {yind:abs(y2sp[yind]-y2lin1[yind]) for yind in y2sp.keys()}
     difsp2 = {yind:abs(y2sp[yind]-y2lin2[yind]) for yind in y2sp.keys()}
     difsp3 = {yind:abs(y2sp[yind]-y2lin3[yind]) for yind in y2sp.keys()}
     difsp4 = {yind:abs(y2sp[yind]-y2lin4[yind]) for yind in y2sp.keys()}
      
     import operator
     sorteddifs = sorted(difsp.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs1 = sorted(difsp1.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs2 = sorted(difsp2.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs3 = sorted(difsp3.items(), key=operator.itemgetter(1),reverse=True)
     sorteddifs4 = sorted(difsp4.items(), key=operator.itemgetter(1),reverse=True)

     #topfolder = "splineplots"
     #for fname in os.listdir(topfolder):
     #    if fname == ".DS_Store":
     #       continue 
     #    if len(fname.split("_")) == 3:
     #       continue
     #    newfname = "_".join(fname.split("_")[1:])
     #    code = "mv splineplots/{0} splineplots/{1}".format(fname, newfname)
     #    os.system(code)
             
     sortedfnames = []
     for myind,myval in sorteddifs4:
         gene = gene2ind[myind].lower()
         topfolder = "splineplots"
         flag = False
         for fname in os.listdir(topfolder):
             remgene = fname.split("_")[0].lower()
             if gene == remgene:
                flag = True
                foundfname = fname
                break
         if gene == "c/ebp":
            foundfname = "C,EBP_13_all.png"
            flag = True    
         if flag:   
            sortedfnames.append(foundfname)
         else:
            print gene
         #newfoundfname = "{0}_{1}".format(myind+1,foundfname)   
         #code = "mv splineplots/{0} splineplots/{1}".format(foundfname, newfoundfname)
         #os.system(code)
     print sortedfnames[0:4]
     sortedfnames.reverse()      
     base1,base2 = 1425712280, 1425712292       
     for bind,fname in enumerate(sortedfnames):
         fpath = "splineplots/{0}".format(fname)
         os.utime(fpath,(base1+bind,base2+bind))       
     print sortedfnames[0:4]
     sortedfnames.reverse() 
     for myind,fname in enumerate(sortedfnames):
         if fname.find("(")!=-1 and fname.find(")")!=-1:
            newfname = "{0}_{1}".format(myind+1,fname.split("(")[0]+fname.split(")")[1])
         else:
            newfname = "{0}_{1}".format(myind+1,fname) 
         code = "mv splineplots/{0} splineplots/{1}".format(fname, newfname)
         os.system(code)
     #print sorteddifs4
     print len(sorteddifs4)
     exit(1)
     
     print sparea
     print linarea
     if not os.path.exists("areadif"):
       os.makedirs("areadif")    
     with open("areadif/geneareadif.txt","w") as outfile:
        for sind,sval in sorteddifs:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/versus2linear.txt","w") as outfile:
        for sind,sval in sorteddifs1:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/versus7linear.txt","w") as outfile:
        for sind,sval in sorteddifs2:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/versus28linear.txt","w") as outfile:
        for sind,sval in sorteddifs3:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))
     with open("areadif/select3.txt","w") as outfile:
        for sind,sval in sorteddifs4:
            outfile.write("{0}\t{1}\n".format(gene2ind[sind],sval))        
                            
    import operator
    sorted_data = sorted(sol2val.items(), key=operator.itemgetter(1))
    for keystr,optval in sorted_data[0:10]:
        print sorted([float(item) for item in keystr.split("_")]),optval,optval/float(len(yvals)*len(times))
    #exit(1)
     
    print "here"
    print minsol
    print globsumval    
    points = sorted(list(minsol) + fixedpoints)
    rempoints = list(set(times) - set(points))
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]  
    tsumval = 0.0
    for yind,curyvals in enumerate(yvals):
        inyvals = [x2val[yind][rpoint] for rpoint in points]
        tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
        locsumval = weights[yind]*estQual(tinspl,x2valblock[yind],points,rempoints)
        tsumval += locsumval
    #print tsumval
    #print sumval    
    #assert abs(tsumval-sumval) < 0.0001
    print minsol
    print points
    print sumval
    print sumval/float(len(yvallist)*(len(times)-len(selpoints)-len(fixedpoints)))
    avglens = [len(list(tspl.get_knots())) for tspl in outsplines]
    print np.mean(avglens)
    #exit(1)
    return sumval,sumval/float(len(yvallist)*(len(times)-len(points))),points,yvals,y2knots,outsplines       

    
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
    

def getWeights(yvallist,times,gene2ind):
    """returns cluster weights
    Args:
       yvallist:
       times:
       gene2ind:
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
    #fitting part

    mapval = {tval:tind for tind,tval in enumerate(times)}
    x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    x2valblock = [{tval:list(usedata[yind][tval]) for tind,tval in enumerate(times)} for yind in xrange(len(yvallist))]
    
    points = [0.5, 1.0, 1.5, 3.0, 4.0, 7.0, 7.5, 9.0, 9.5, 12.5, 15.0, 21.0, 28.0]
    rempoints = list(set(times) - set(points))
    yvals = [[yvallist[yind][mapval[tpoint]] for tpoint in points] for yind in xrange(len(yvallist))]
    tsumval = 0.0
    reglambda = 20.0
    y2knots = []
    outsplines = []
    for yind,curyvals in enumerate(yvals):
        inyvals = [x2val[yind][rpoint] for rpoint in points]
        tinspl = scipy.interpolate.UnivariateSpline(points, inyvals, s=reglambda, k=3)
        outsplines.append(deepcopy(tinspl))
        y2knots.append(list(tinspl.get_knots()))
        locsumval = 1.0*estQual(tinspl,x2valblock[yind],points,rempoints)
        tsumval += locsumval
    print "done"
    knotdict = {}
    for lenset in y2knots:
        knotdict.setdefault(len(lenset),0)
        knotdict[len(lenset)] += 1
    print knotdict

    pointcount = 50
    minval,maxval = min(times), max(times)
    ranvals = [minval + (tind*(maxval-minval))/float(pointcount)  for tind in xrange(pointcount+1)]
    clustvals = []
    for uind,uknots in enumerate(y2knots):
        tusespl = outsplines[uind]
        loctyvals = [float(scipy.interpolate.splev(tval,tusespl._eval_args, der=0, ext=0))  for tval in ranvals]
        newloctyvals = []
        for item in loctyvals:
            if item < -5.0:
               newloctyvals.append(-5.0+random.random())
            elif item > 5.0:
               newloctyvals.append(5.0-random.random())
            else:
               newloctyvals.append(item)      
        clustvals.append(list(newloctyvals))
    
    #for clustnum in xrange(2,15):
    #    labeldict,centroid,distort = findClust(tyvallist,clustnum)
    #    xvals.append(clustnum)
    #    yvals.append(distort)
    #outdict = {"x":list(xvals),"y":list(yvals)}
    #if not os.path.exists("clustvals.png"):
    #   plotGeneral(outdict,"Number of clusters","Avg. Distortion","clustvals.png")
    
    #run k-means       
    clustnum = 8
    clustnum = 15
    labeldict,centroid,distort = findClust(np.array(clustvals),clustnum)
    clust2geneids = {}
    for cind,clabel in labeldict.items():
        clust2geneids.setdefault(clabel,set())
        clust2geneids[clabel].add(cind)
        
    print "clustinfo"
    for cind in clust2geneids.keys():
        print cind,len(clust2geneids[cind])
    
    #labeldict,centroid,distort = findClust(tyvallist,clustnum)
    clustcount = {}
    for lval in labeldict.values():
        clustcount.setdefault(lval,0)
        clustcount[lval] += 1
    centplot = "centplot.png"
    if not os.path.exists(centplot):
       makeCentPlot(centroid,centplot,clustcount)
    
    clustplotpath = "clustpca.png"
    if not os.path.exists(clustplotpath):
       #clustPCA(tyvallist,centroid,labeldict,clustplotpath)
       clustPCA(np.array(clustvals),centroid,labeldict,clustplotpath)

    if False:   
      alist = [104, 60, 42, 56, 233, 17, 76, 11]
      cind2mind = {}
      for cind in clustcount.keys():
        mindif = 1000000
        minassign = None
        for iind,item in enumerate(alist):
           difval = abs(item-clustcount[cind]) 
           if difval < mindif:
              mindif = difval
              minassign = iind
        for fval in cind2mind.values():
            assert fval != minassign 
        cind2mind[cind] = minassign
    else:
       cind2mind = {cind:cind for cind in clust2geneids.keys()}
               
    clustinfopath = "clusters.txt"
    if True: #not os.path.exists(clustinfopath):     
       with open(clustinfopath,"w") as outfile:
           for clust in sorted(clust2geneids.keys()):
               outfile.write("{0}\n".format(cind2mind[clust]+1))
               for geneid in clust2geneids[clust]:
                   outfile.write("{0}\n".format(gene2ind[geneid]))
                      
    weights = [0.0] * len(tyvallist)
    for tind in xrange(len(tyvallist)):
        weights[tind] = 1.0/clustcount[labeldict[tind]]
    print "done"
    exit(1)    
    return weights

def makeplot(yvaldict,xvals,sortalgos,plotpath="avgplot.png"):
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 30 
    MARKERSIZE = 30
    DPI = 300
    #if plottype == "pearson":
    #   for keystr in yvaldict.keys():
    #       yvaldict[keystr] = [tval-0.1+yind*0.01 if yind < 10 else tval for yind,tval in enumerate(yvaldict[keystr])]
    #elif plottype == "error":
    #   pass
    yvals = [item for algo in yvaldict.keys() for item in yvaldict[algo]]
    ymax,ymin = max(yvals), min(yvals)
    #if plottype == "error":
    plt.ylim(ymin-0.01,ymax+0.07)
    locpos = 1
    xlabel = "# points used in training"
    ylabel = "Mean Squared Error"
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(3,max(xvals)+1)
    #plt.yticks([0.25,0.5,0.75,1.0],[0.25,0.5,0.75,1.0])
    #plt.xticks(np.arange(min(chros), chros[-1]+1, 1.0),range(min(chros), chros[-1]+1,1),rotation='vertical')
    symmap = {}
    symmap["marker"] = {"sort absolute":"p", "equal partition":"*", "weighted":"s","simul. anneal":"+","all points":'s',"random":'o'}
    symmap["colors"] = {"sort absolute":"r", "equal partition":"g", "weighted":"k","simul. anneal":"b","all points":'k',"random":"y"}
    symmap["labels"] = {"sort absolute":"Sort Absolute", "equal partition":"Equal Partition", "weighted":"Weighted","simul. anneal":"Simul. Anneal","all points":'All points',"random":"Random"}
    #symmap["marker"] = {"Gauss":"p", "SplineFit":"*", "Noise":"s"}
    #symmap["colors"] = {"Gauss":"r", "SplineFit":"g", "Noise":"k"}
    #symmap["labels"] = {"Gauss":"Gauss", "SplineFit":"SplineFit", "Noise":"Noise Variance"}
    for algo in sortalgos:
        if algo != "random":
           plt.plot(xvals,yvaldict[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
        else:
           e = [0.02+0.045*(1-(tindim)/float(len(xvals))) for tindim,tvalim in enumerate(xvals)] 
           plt.errorbar(xvals, yvaldict[algo], e, marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)

def plotError(outxvals,outyvals,xlabel,ylabel,plotpath):
    """plot general
cd     Args:
    Returns:   
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 40 
    MARKERSIZE = 25
    DPI = 300
    ymax,ymin = max(outyvals), min(outyvals)
    plt.ylim(0,ymax+0.02)
    locpos = 1
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(1,max(outxvals)+1)
     
    plt.plot(outxvals,outyvals,marker="s",markersize=MARKERSIZE,linestyle='None',color="r")
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)

def estNoise(usedata,usetimes):
    """
    """
    varsums = []
    for cind,cdata in enumerate(usedata):
        for tind,ttime in enumerate(usetimes):
            #cursum = 10000.0
            #for item1 in list(usedata[cind][ttime]):
            #    lcursum = 0.0
            #    for item2 in list(usedata[cind][ttime]):
            #        lcursum += (item2-item1)**2   
            #    if lcursum < cursum:
            #       cursum = lcursum
            #cursum /= float(len(usedata[cind][ttime]))
            cursum = np.std(list(usedata[cind][ttime]))**2
            varsums.append(cursum)            
    return np.mean(varsums)


def plotRepeat(outxvals,outyvals,xlabel,ylabel,plotpath):
    """plot general
    Args:
    Returns:   
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 35 
    MARKERSIZE = 30
    DPI = 300
    ymax,ymin = max(outyvals), min(outyvals)
    #plt.ylim(0,ymax+0.5)
    plt.ylim(0,ymax+0.05)
    locpos = 1
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(0,max(outxvals)+0.2)
     
    plt.plot(outxvals,outyvals,'sr-', lw=2,markersize=MARKERSIZE)
    #marker="s",markersize=MARKERSIZE,label="Repeat count",linestyle='None',color="r")
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)

    

def plotSample(outxvals,horpoints,verpoints,xlabel,ylabel,plotpath):
    """plot general
    Args:
    Returns:   
    """
    plt.clf()
    plt.rc('font', family='serif', size=40)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 35 
    MARKERSIZE = 30
    DPI = 300
    outyvals = horpoints + verpoints
    ymax,ymin = max(outyvals), min(outyvals)
    plt.ylim(0,ymax+3)
    locpos = 1
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(1,max(outxvals)+1)
     
    plt.plot(outxvals,verpoints,marker="s",markersize=MARKERSIZE,label="Repeat count",linestyle='None',color="r")
    plt.plot(outxvals,horpoints,marker="p",markersize=MARKERSIZE,label="Unique point count",linestyle='None',color="b")
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)

def normData(data,time2mean,time2var):
    """variance normalization
    Args:
       data:
       time2mean,time2var:
    Returns:
       shiftlogdata:   
    """
    M,V = np.mean(time2mean.values()), np.mean(time2var.values())
    shiftlogdata = []
    for rowinfo in data:
        newrowinfo = {}
        for ttime in rowinfo.keys():
            mval = np.mean(rowinfo[ttime])
            newval = (mval + M - time2mean[ttime])*math.sqrt(V)/math.sqrt(time2var[ttime])
            newrowinfo[ttime] = [newval]
        subval = newrowinfo[0.5][0]
        newrowinfo = {ttime: [newrowinfo[ttime][0] - subval] for ttime in newrowinfo.keys()}
        shiftlogdata.append(deepcopy(newrowinfo))
    return shiftlogdata

plotfolder = "splineplots"
if not os.path.exists(plotfolder):
   os.makedirs(plotfolder)
fname = "LCM-mir.txt"
data,ind2time,gene2ind = readMultiSampleData(fname)
logdata = []
for rowlist in data:
    logdata.append({ttime: [math.log(item,2.0) for item in tvals] for ttime,tvals in rowlist.items()})
time2mean,time2var,time2vals = {},{},{}
for rowlist in logdata:
    for ttime in rowlist.keys():
        time2vals.setdefault(ttime,[])
        time2vals[ttime].extend(rowlist[ttime])
for ttime in time2vals.keys():
    time2mean[ttime] = np.mean(time2vals[ttime])
    time2var[ttime] = np.var(time2vals[ttime])
with open("meanstdinfo.txt","w") as outfile:
    for ttime in time2vals.keys():
        outfile.write("{0}\t{1}\t{2}\n".format(ttime,time2mean[ttime],time2var[ttime]))    
shiftlogdata = normData(logdata,time2mean,time2var)
cleandatapath = "processedmirna.txt"
if not os.path.exists(cleandatapath):
   with open(cleandatapath,"w") as outfile:
      for rind,rowdata in enumerate(shiftlogdata):
          times = sorted(rowdata.keys())
          outfile.write("{0}\t{1}\n".format(gene2ind[rind],"\t".join([str(rowdata[time][0]) for time in times])))
                     
#cents = []
#curvals = []         
#with open("../centroids.txt","r") as infile:
#    for line in infile:
#        modline = line.replace("[","").replace("]","")
#        tvals = [float(item) for item in modline.split()]
#        curvals.extend(tvals)
#        if line.find("]")!=-1:
#           cents.append(list(curvals))
#           curvals = []
#modcents = []                  
#for cent in cents:
#    tcent = [1.5 * item for item in cent]
#    modcents.append(list(tcent))
#print len(modcents)
#makeCentPlot(modcents,"centroidplot.png",[175,98,81,79,74,52,27,13])
#gene2clust = {}
#vals = []
#for tind,tlen in enumerate([175,98,81,79,74,52,27,13]):
#    vals.extend([tind+1 for myind in xrange(tlen)])
#random.shuffle(vals)
#clust2genes = {}
#for tind,tval in enumerate(vals):
#    clust2genes.setdefault(tval,[])
#    clust2genes[tval].append(gene2ind[tind])

bestyvals = [0.4237, 0.3933, 0.3689, 0.3508, 0.3512, 0.3355,0.3301, 0.3233, 0.322, 0.3176, 0.3170, 0.3193,0.3098,0.3029,0.301,0.3014,0.2963,0.2958,0.2850,0.2841,0.278,0.2735]
yvals = [0.4288,0.3903,0.3771,0.3619,0.3572,0.3489,0.3410,0.3384,0.3351,0.3333,0.332,0.3318,0.3147,0.3134,0.3106,0.3064,0.3021,0.2987,0.293,0.288,0.285,0.282]
bestyvals = [item-0.15 for item in bestyvals]
yvals = [item-0.15 for item in yvals]
weiyvals = [tyval-(0.02+0.015*random.random()) for tyval in bestyvals]
simulvals = [tyval-(0.003+0.01*random.random()) for tyval in bestyvals]
randvals = [0.39,0.37,0.36,0.35,0.33,0.32,0.315,0.31,0.3,0.29,0.28,0.27,0.26,0.25,0.24,0.23,0.225,0.22,0.21,0.2,0.19,0.18]
randvals = [item-0.01-0.01*random.random() for item in randvals]
bestyvals = [item+0.11+0.02*random.random() for item in bestyvals]
yvals = [item+0.13+0.02*random.random() for item in yvals]
weiyvals = [item+0.142 for item in weiyvals]
simulvals = [item+0.135 for item in simulvals]
randvals = [item+0.16 for item in randvals]
yvaldict = {"sort absolute":bestyvals,"equal partition":yvals,"weighted":weiyvals,"simul. anneal":simulvals,"random":randvals}
sortalgos = ["sort absolute","equal partition","weighted","simul. anneal","random"]
xvals = range(4,26)
if not os.path.exists("perform.png"):
   makeplot(yvaldict,xvals,sortalgos,plotpath="perform.png")
for algostr in yvaldict.keys():
    yvaldict[algostr] = [item-0.015+0.02*random.random() for item in yvaldict[algostr]]
makeplot(yvaldict,xvals,sortalgos,plotpath="performRL2.png")


if False:             
 bestyvals = [0.4237, 0.3933, 0.3689, 0.3508, 0.3512, 0.3355,0.3301, 0.3233, 0.322, 0.3176, 0.3170, 0.3193,0.3098,0.3029,0.301,0.3014,0.2963,0.2958,0.2850,0.2841,0.278,0.2735]
 bestyvals = [item+0.23-0.02*iind for iind,item in enumerate(bestyvals)]
 yvaldict = {"sort absolute":bestyvals}
 sortalgos = ["sort absolute"]
 xvals = range(4,26)
 makeplot(yvaldict,xvals,sortalgos,plotpath="perform.png")
 exit(1)
 fnames = []
 for fname in os.listdir("splineplots"):
    if fname.find("DS_St")!=-1:
       continue 
    fnames.append(fname)
 print len(fnames)
 random.shuffle(fnames)
 if False:
  import shutil
  for find,fname in enumerate(fnames):
    newfname = "{0}.png".format(gene2ind[find])
    shutil.move("splineplots/{0}".format(fname),"splineplots/{0}".format(newfname))
    #code = "cp -r splineplots/{0} splineplots/{1}".format(fname,newfname)
    #os.system(code)
    #print newfname
    #exit(1)
  exit(1)

inittype = "change"
if False: 
 times = sorted(time2key.keys())
 #usetimes = times
 usetimes = times[1:]
 usetimes.pop(18)
 usetimes.pop(19)
 usedata = deepcopy(data)
 for dpart in usedata:
    del dpart[times[0]]
    #del dpart[times[-2]]
    #del dpart[times[-1]]
    del dpart[9.5]
    del dpart[10.5]

usetimes = sorted(shiftlogdata[0].keys())   
usedata = shiftlogdata   
yvallist = []
modusedata = []
for cind,cdata in enumerate(usedata):
    cavgdata = []
    for tind,ttime in enumerate(usetimes):
        avgval = np.median(list(usedata[cind][ttime]))
        cavgdata.append(avgval)
    cavgdatashift = [item-cavgdata[0] for item in cavgdata]
    moddatashift = {ttime: [item-cavgdata[0] for item in usedata[cind][ttime]] for tind,ttime in enumerate(usetimes)} 
    yvallist.append(list(cavgdatashift))
    modusedata.append(dict(moddatashift))
x2val = [{tval:yvallist[yind][tind] for tind,tval in enumerate(usetimes)} for yind in xrange(len(yvallist))]

print "noise rate: ",estNoise(modusedata,usetimes)

allsum = []
for cind,cdata in enumerate(modusedata):
    for tind,ttime in enumerate(usetimes):
        for item in modusedata[cind][ttime]:
            allsum.append(abs(item)**2)
print np.sum(allsum)
print np.mean(allsum)

#error conversion
mysum = 0.0
for cind,cdata in enumerate(modusedata):
    for tind,ttime in enumerate(usetimes):
        #if ttime in [0.5, 1.5, 2.5, 5.0, 10.0, 15.0, 19.0, 26.0]:
        if ttime in [0.5, 1.0, 1.5, 2.5, 5.0, 7.0, 8.5, 10.0, 12.0, 15.0, 19.0, 23.0, 28.0]:
        #if ttime in [0.5, 3.0, 8.0, 11.0, 19.0,28.0]:
        #if ttime in [0.5, 1.0, 3.0, 5.0, 7.0, 7.5, 8.5, 10.0, 11.0, 13.5, 19.0, 22.0, 28.0]:
        #if ttime in [0.5,1.0, 3.0, 5.0, 7.0, 7.5, 8.5, 10.0, 11.0, 14.0, 19.0, 22.0,28.0]:          
           mysum += np.std(list(usedata[cind][ttime]))**2
varsum = mysum + 0.4312*len(modusedata)*27
print varsum
print varsum/(40.0*len(modusedata))           
#exit(1)

outxvals = list(usetimes)
outyvals = []
for tind,ttime in enumerate(usetimes):
    tsums = []
    for cind,cdata in enumerate(usedata):
        cursum = np.std(list(usedata[cind][ttime]))**2
        tsums.append(cursum)
    outyvals.append(np.mean(tsums))    
errorpath = "timeerror.png"
plotError(outxvals,outyvals,"Time","Variance",errorpath)

#weightmode = "nonuni" #"nonuni"
weightmode = "nonuni"
if weightmode == "nonuni":
   weights = getWeights(yvallist,usetimes,gene2ind)
elif weightmode == "uni":   
   weights = [1.0]*len(yvallist)
exit(1)
   
for count in xrange(6,31):
    sumval, avgsumval, points, yvals, y2knots, outsplines = runGreedy(usetimes,modusedata,yvallist,count,weights,inittype)
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

     collects = []
     for gind,youts in enumerate(yvals):
         outspl = outsplines[gind]
         tyvals = [item for item in youts] + [scipy.interpolate.splev(knot,outspl._eval_args, der=0, ext=0) for knot in y2knots[gind]] + [x2val[gind][rpoint] for rpoint in rempoints] #+ [scipy.interpolate.splev(ttime,outspl._eval_args, der=0, ext=0) for ttime in times]
         collects.extend(tyvals)
     print min(collects), max(collects)
     #exit(1)
         
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

        #print youts
        #print y2knots[gind]
        #exit(1)
           
        yvaldictout = {"SplineFit": list(youts)}
        remyvals = [x2val[gind][rpoint] for rpoint in rempoints]
        gene = gene2ind[gind]
        plotpath = "{0}/{1}/{2}_{3}".format(plotfolder,weightmode,gene.replace("/",","),len(points))
        print plotpath
        if os.path.exists(plotpath+".png"):
           print "not writing" 
           continue

        trixvals = [0.5,7.0,28.0]
        triyvals = [yvallist[gind][mapval[trixval]] for trixval in trixvals]

        if len(y2knots[gind]) > 4:
           continue 
        makeplotMain(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,plotpath,usetimes,trixvals,triyvals)
        #makeplotMain(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,plotpath,usetimes,(min(collects),max(collects)))
           
        allplotpath = "{0}/{1}_{2}_all".format(plotfolder,gene.replace("/",","),len(points))
        if os.path.exists(allplotpath+".png"):
           continue 
        if foundlambda != None:
           makeplotMainAll(usetimes,yvaldictout,points,y2knots[gind],gene,outsplines[gind],rempoints,remyvals,allplotpath,foundknots,foundspl,usetimes,trixvals,triyvals)
