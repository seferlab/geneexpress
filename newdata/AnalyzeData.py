import math
import numpy as np
import matplotlib
import random
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import scipy.cluster

def plotError(outxvals,outyvals,xlabel,ylabel,plotpath):
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
    LEGENDSIZE = 40 
    MARKERSIZE = 15
    DPI = 300
    ymax,ymin = max(outyvals), min(outyvals)
    plt.ylim(0,ymax+0.02)
    locpos = 1
    plt.xlabel(xlabel,fontsize=FSIZE)
    plt.ylabel(ylabel,fontsize=YFSIZE)
    plt.xlim(0,max(outxvals)+1)
     
    plt.plot(outxvals,outyvals,marker="s",markersize=MARKERSIZE,linestyle='None',color="r")
    ax = plt.axes()        
    ax.xaxis.grid()        
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.13)    
    plt.savefig(plotpath, dpi=DPI)

    
def analyzeNoise(usedata,usetimes):
    """noise analyze and plot gen
    Args:
       usedata:
       usetimes:
    Returns:
    """
    varsums = []
    for cind,cdata in enumerate(usedata):
        for tind,ttime in enumerate(usetimes):
            cursum = np.std(list(usedata[cind][ttime]))**2
            varsums.append(cursum)      
    print np.mean(varsums)
    
    outxvals = list(usetimes)
    outyvals = []
    outyvals2 = []
    for tind,ttime in enumerate(usetimes):
        tsums = []
        for cind,cdata in enumerate(usedata):
            cursum = 10000000000.0
            for item1 in list(usedata[cind][ttime]):
                lcursum = 0.0
                for item2 in list(usedata[cind][ttime]):
                    lcursum += (item2-item1)**2   
                if lcursum < cursum:
                   cursum = lcursum
            cursum /= float(len(usedata[cind][ttime]))**2
            tsums.append(cursum)
        outyvals.append(np.mean(tsums))
        outyvals2.append(np.mean([np.std(list(usedata[cind][ttime]))**2 for cind,cdata in enumerate(usedata)]))
    print outyvals2
    time2noise = {}
    for tind,ttime in enumerate(usetimes):
        time2noise[ttime] = outyvals[tind]
    sortedkeys = sorted(time2noise, key=lambda key: time2noise[key],reverse=True)
    for tkey in sortedkeys:
        print tkey,time2noise[tkey]   
    print np.mean(outyvals2)
    print "my noise is: ",np.mean(outyvals)
    errorpath = "timeerror.png"
    plotError(outxvals,outyvals,"Days","Noise (Variance)",errorpath)
