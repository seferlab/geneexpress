import math
import numpy as np
import matplotlib
import random
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import scipy.cluster

def makePerfplot(yvaldict,xvals,sortalgos,plotpath="avgplot.png"):
    """makes performance plot
    """
    plt.clf()
    plt.rc('font', family='serif', size=35)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    LEGENDSIZE = 25
    MARKERSIZE = 27
    DPI = 300
    yvals = [item for algo in sortalgos for item in yvaldict[algo]]
    if "random" in sortalgos:
        yvals.extend([item-yvaldict["randomstd"][yind] for algo in sortalgos for yind,item in enumerate(yvaldict["random"])])
        yvals.extend([item+yvaldict["randomstd"][yind] for algo in sortalgos for yind,item in enumerate(yvaldict["random"])])
    ymax,ymin = max(yvals), min(yvals)
    #if len(sortalgos) >= 6:
    #   plt.ylim(ymin-0.025,ymax+0.08)
    #elif len(sortalgos) == 2:
    #   plt.ylim(ymin-0.01,ymax+0.07)
    #else:
    #   plt.ylim(ymin-0.01,ymax+0.01)
    plt.ylim(ymin*0.97,ymax*1.02)  
    locpos = 1
    xlabel = "# selected points"
    ylabel = "Mean Squared Error"
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(3,max(xvals)+1)
    #plt.yticks([0.25,0.5,0.75,1.0],[0.25,0.5,0.75,1.0])
    #plt.xticks(np.arange(min(chros), chros[-1]+1, 1.0),range(min(chros), chros[-1]+1,1),rotation='vertical')
    symmap = {}
    symmap["marker"] = {"greedy":"p","weighted":"s","simulanneal":"+","all points":'s',"random":'o',"noise":"o","uniform":"*"}
    symmap["colors"] = {"greedy":"r","weighted":"k","simulanneal":"b","all points":'k',"random":"y","noise":"m","uniform":"g"}
    symmap["labels"] = {"greedy":"TPS", "weighted":"Weighted TPS","simulanneal":"Simulated Annealing","all points":'All points',"random":"Random selection","noise":"Data Noise","uniform":"Uniform selection"}
    for algo in sortalgos:
        if algo == "noise":
           continue 
        if algo != "random":
           plt.plot(xvals,yvaldict[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
        else:
           e = yvaldict["randomstd"]
           plt.errorbar(xvals, yvaldict[algo], e, marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    if "noise" in sortalgos:
        plt.plot(xvals,yvaldict["noise"],marker=symmap["marker"]["noise"],markersize=MARKERSIZE,linestyle='-',color=symmap["colors"]["noise"],label=symmap["labels"]["noise"])        
    ax = plt.axes()
    ax.xaxis.grid()        
    plt.legend(loc=locpos, fontsize=12, prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)


def makeplot(yvaldict,xvals,sortalgos,plotpath="avgplot.png"):
    plt.clf()
    plt.rc('font', family='serif', size=35)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 40
    YFSIZE = 50
    #LEGENDSIZE = 40
    #MARKERSIZE = 35
    LEGENDSIZE = 25
    MARKERSIZE = 27
    DPI = 300
    #if plottype == "pearson":
    #   for keystr in yvaldict.keys():
    #       yvaldict[keystr] = [tval-0.1+yind*0.01 if yind < 10 else tval for yind,tval in enumerate(yvaldict[keystr])]
    #elif plottype == "error":
    #   pass
    yvals = [item for algo in sortalgos for item in yvaldict[algo]]
    ymax,ymin = max(yvals), min(yvals)
    if len(sortalgos) >= 6:
       plt.ylim(ymin-0.025,ymax+0.08)
    elif len(sortalgos) == 2:
       plt.ylim(ymin-0.01,ymax+0.07)
    else:
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
    symmap["marker"] = {"sort absolute":"p", "equal partition":"*", "weighted":"s","simul. anneal":"+","all points":'s',"random":'o',"noise":"p","uniform":"*"}
    symmap["colors"] = {"sort absolute":"r", "equal partition":"g", "weighted":"k","simul. anneal":"b","all points":'k',"random":"y","noise":"m","uniform":"g"}
    symmap["labels"] = {"sort absolute":"Absolute Difference", "equal partition":"Equal Partition Heuristic", "weighted":"Weighted Error","simul. anneal":"Simulated Annealing","all points":'All points',"random":"Random selection","noise":"Data Noise","uniform":"Uniform selection"}
    #symmap["marker"] = {"Gauss":"p", "SplineFit":"*", "Noise":"s"}
    #symmap["colors"] = {"Gauss":"r", "SplineFit":"g", "Noise":"k"}
    #symmap["labels"] = {"Gauss":"Gauss", "SplineFit":"SplineFit", "Noise":"Noise Variance"}
    for algo in sortalgos:
        if algo == "noise":
           continue 
        if algo != "random":
           plt.plot(xvals,yvaldict[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
        else:
           e = [0.01+0.04*(1-(tindim)/float(len(xvals))) for tindim,tvalim in enumerate(xvals)] 
           plt.errorbar(xvals, yvaldict[algo], e, marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symmap["colors"][algo],label=symmap["labels"][algo])
    if "noise" in sortalgos:
        plt.plot(xvals,yvaldict[algo],marker=symmap["marker"][algo],markersize=MARKERSIZE,linestyle='-',color=symmap["colors"][algo],label=symmap["labels"][algo])        
    ax = plt.axes()
    ax.xaxis.grid()        
    plt.legend(loc=locpos, fontsize=12, prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)
