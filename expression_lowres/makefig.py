import os
import sys
import itertools
import numpy as np
import math
import networkx as nx
import random
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cPickle as cpickle
#sys.path.append("../../lib")
#sys.path.append("..")
#import EmbedUtilities
#import EvalUtilities
#import myutilities as myutil
#import PlotUtilities


def getSpatial():
    """
    """
    size = 10
    count = 5
    data = {ind:np.zeros((size,size),dtype=np.int) for ind in xrange(count)}
    sumdata = np.zeros((size,size),dtype=np.int)
    lens = {0:[(0,1),(4,8)],1:[(0,2),(4,9)],2:[(1,1),(3,8)],3:[(0,2),(4,5),(7,8)],4:[(0,1),(2,2),(4,7)]}
    for cind in xrange(count):
        for start,end in lens[cind]:
            data[cind][start:end+1,start:end+1] += 1
        p = 0.1
        for ind1,ind2 in itertools.combinations(range(size),2):
            if random.random() < p:
               data[cind][ind1,ind2] += 1
               data[cind][ind2,ind1] += 1
    for cind in data.keys():
        sumdata[0:size,0:size] += data[cind][0:size,0:size]
    return sumdata,data


def getTemporal():
    """
    """
    size = 10
    count = 4
    phases = ["G1","S","G2","M"]
    data = {ind:np.zeros((size,size),dtype=np.int) for ind in xrange(count)}
    tempdata = np.zeros((size,size),dtype=np.int)
    lens = {0:[(0,5),(7,8)],1:[(0,3),(6,8)],2:[(1,4),(6,9)],3:[(1,3),(4,6),(7,7)]}
    for cind in xrange(count):
        for start,end in lens[cind]:
            data[cind][start:end+1,start:end+1] += 1
        p = 0.1
        for ind1,ind2 in itertools.combinations(range(size),2):
            if random.random() < p:
               data[cind][ind1,ind2] += 1
               data[cind][ind2,ind1] += 1
    for cind in data.keys():
        tempdata[0:size,0:size] += data[cind][0:size,0:size]    
    p = 0.2
    for ind1,ind2 in itertools.combinations(range(size),2):
        if random.random() < p:
           tempdata[ind1,ind2] += 1
           tempdata[ind2,ind1] += 1
    return tempdata,data,phases


def makespaPlot(spadata,spasum):
    mycmap = plt.cm.Reds
    plt.clf()
    plt.rc('font', family='serif', size=25)
    #plt.subplots_adjust(left=0.01, right=0.99)

    maxval = np.amax(spasum) 
    fig, axes = plt.subplots(nrows=1,ncols=6)
    fig.set_size_inches(36,10)
    for cind in spadata.keys():
        ax = axes.flat[cind]
        im = ax.pcolor(spadata[cind],cmap=mycmap,alpha=0.8,edgecolors='k',linewidths=0.1,vmin=0, vmax=maxval)
        ax.set_xlabel("Class {0}".format(cind+1), fontsize = 40)    
    ax = axes.flat[5]
    im = ax.pcolor(spasum,cmap=mycmap,alpha=0.8,edgecolors='k',linewidths=0.1,vmin=0, vmax=maxval)
    ax.set_xlabel("Ensemble",fontsize = 40) 
    
    #plt.text(-10,-5,"+",fontsize = 60)c
    plt.subplots_adjust(left=0.02, right=0.98)
    plt.savefig("decompspa.png",DPI=800)

    
def makegenPlot(spadata,spasum):
    del spadata[4]
    mycmap = plt.cm.Reds
    plt.clf()
    plt.rc('font', family='serif', size=25)
    maxval = np.amax(spasum) 
    fig, axes = plt.subplots(nrows=1,ncols=5)
    fig.set_size_inches(54,10)
    for cind in spadata.keys():
        ax = axes.flat[cind]
        im = ax.pcolor(spadata[cind],cmap=mycmap,alpha=0.8,edgecolors='k',linewidths=0.1,vmin=0, vmax=maxval)
        ax.set_xlabel("Class {0}".format(cind+1), fontsize = 30)    
    ax = axes.flat[4]
    im = ax.pcolor(spasum,cmap=mycmap,alpha=0.8,edgecolors='k',linewidths=0.1,vmin=0, vmax=maxval)
    ax.set_xlabel("Ensemble",fontsize = 30)
    
    plt.text(-49.5,5,r"$\lambda_%s$" %1,fontsize = 60)
    plt.text(-38.2,5,r"$+\lambda_%s$" %2,fontsize = 60)
    plt.text(-26,5,r"$+\lambda_%s$" %3,fontsize = 60)
    plt.text(-13.9,5,r"$+\lambda_%s$" %4,fontsize = 60)
    plt.text(-1.6,5,"=",fontsize = 60)
    plt.subplots_adjust(left=0.05, right=0.95)
    plt.savefig("fdproblem.png",DPI=800)
    
    
def maketempPlot(tempdata,tempsum,phases):
    mycmap = plt.cm.Reds
    plt.clf()
    plt.rc('font', family='serif', size=25)
    #plt.subplots_adjust(left=0.01, right=0.99)

    maxval = np.amax(tempsum)
    fig, axes = plt.subplots(nrows=1,ncols=5)
    fig.set_size_inches(36,10)
    for cind in tempdata.keys():
        ax = axes.flat[cind]
        im = ax.pcolor(spadata[cind],cmap=mycmap,alpha=0.8,edgecolors='k',linewidths=0.1,vmin=0, vmax=maxval)
        ax.set_xlabel(phases[cind],fontsize = 60)  
    ax = axes.flat[4]
    im = ax.pcolor(spasum,cmap=mycmap,alpha=0.8,edgecolors='k',linewidths=0.1,vmin=0, vmax=maxval)
    ax.set_xlabel("Ensemble",fontsize = 40) 

    #plt.text(-10,-5,"+",fontsize = 40)
    #plt.arrow()
    plt.subplots_adjust(left=0.02, right=0.98)
    plt.savefig("decomptemp.png",DPI=800)


def explot(yvallist,usetimes,gene2ind):
    """
    """    
    time2ind,ind2time = {},{}
    for tind,ttime in enumerate(usetimes):
        time2ind[ttime] = tind
        ind2time[tind] = ttime
    plt.clf()
    plt.rc('font', family='serif', size=9)
    fig = plt.figure()
    plt.subplots_adjust(left=0.07, right=0.93,top=0.95,bottom=0.01)
    ax1 = plt.subplot2grid((10,11),(0, 0),colspan=3,rowspan=4)
    ax2 = plt.subplot2grid((10,11),(0, 4),colspan=3,rowspan=4)
    ax3 = plt.subplot2grid((10,11),(0, 8),colspan=3,rowspan=4)
    ax4 = plt.subplot2grid((10,11),(5, 0),colspan=3,rowspan=4)
    ax5 = plt.subplot2grid((10,11),(5, 4),colspan=3,rowspan=4)
    ax6 = plt.subplot2grid((10,11),(5, 8),colspan=3,rowspan=4)

    from matplotlib.patches import Rectangle,Circle
    #ax2.axes.get_xaxis().set_visible(False)
    #ax2.axes.get_yaxis().set_visible(False)
    #ax3.axes.get_xaxis().set_visible(False)
    #ax3.axes.get_yaxis().set_visible(False)
    
    for tax in [ax1,ax4]:
        tax.set_xlabel("Days")
        tax.set_ylabel("Expression")

    from matplotlib.patches import FancyArrow
    for xpos,ypos,xdif,ydif in [(0.3,0.75,0.03,0.0),(0.62,0.75,0.04,0.0),(0.87,0.535,0.0,-0.04),(0.35,0.3,-0.03,0.0),(0.67,0.3,-0.04,0.0)]:
        curarr = FancyArrow(xpos,ypos,xdif,ydif,head_width=0.01,head_length=0.01,fc='k',ec='k',transform=ax1.figure.transFigure,clip_on=False)
        ax1.add_patch(curarr)        
    #plt.text(-1.1,2.35,r"Initial points",fontsize = 10)
    #plt.text(0.3,2.35,r"Initial fit",fontsize = 10)
    #plt.text(-1.1,1.05,r"Final fit",fontsize = 10)
    #plt.text(0.3,1.05,r"Iterations 2-7",fontsize = 10)

    gind = 59
    if gind == 0:
       ax2.text(6.5,0.85,r"Initial points",fontsize = 10)
       ax3.text(8.5,1.55,r"Initial fit",fontsize = 10)
       ax6.text(3.5,1.05,r"Iterations 2-7",fontsize = 10)
       ax5.text(7.5,1.05,r"Final fit",fontsize = 10)

       plt.text(-90,1.35,"a)",fontsize = 13)
       plt.text(-46,1.35,"b)",fontsize = 13)
       plt.text(-2,1.35,"c)",fontsize = 13)
       plt.text(-90,-1.85,"f)",fontsize = 13)
       plt.text(-46,-1.85,"e)",fontsize = 13)
       plt.text(-2,-1.85,"d)",fontsize = 13)
    elif gind == 59:
       ax1.text(2.5,4.35,r"Expression profiles",fontsize = 10)
       ax4.text(2.5,4.35,r"Reconstructed profiles",fontsize = 10)
       ax2.text(6.5,3.25,r"Initial points",fontsize = 10)
       #ax3.text(8.5,72,r"Initial fit",fontsize = 10)
       ax3.text(8.5,3.3,r"Initial fit",fontsize = 10)
       ax6.text(3.5,3.25,r"Iterations 2-7",fontsize = 10)
       ax5.text(7.5,3.25,r"Final fit",fontsize = 10)

       plt.text(-90,4.55,"a)",fontsize = 13)
       plt.text(-46,4.55,"b)",fontsize = 13)
       plt.text(-2,4.55,"c)",fontsize = 13)
       plt.text(-90,-8.45,"f)",fontsize = 13)
       plt.text(-46,-8.45,"e)",fontsize = 13)
       plt.text(-2,-8.45,"d)",fontsize = 13)
        
    import sys
    sys.path.append("../lib")
    import Tempselect
    usedata = [{ttime:[rowdata[tind]] for tind,ttime in enumerate(usetimes)} for rowdata in yvallist]

    usecolors = ["r","b","g","m","y","c","k"]
    removeset = []
    for yind,rowdata in enumerate(yvallist):
        flag = False
        for item in rowdata:
            if item <= -4.0:
               flag = True
               break
        if flag:
           removeset.append(yind)    
        else:
           ax1.plot(usetimes,rowdata,color=usecolors[yind%len(usecolors)])
    count = 8
    weights = [1.0]*len(yvallist)
    
    #iterpointsall,y2knotslist,outsplineslist,avgerrorlist = Tempselect.runGreedyModified(usetimes,usedata,count,weights,"equal",1)
    savepath = "savetrace.pkl"
    with open(savepath,"rb") as infile:
        iterpointsall = cpickle.load(infile)
        y2knotslist = cpickle.load(infile)
        outsplineslist = cpickle.load(infile)
        avgerrorlist = cpickle.load(infile)

    import geneUtilities
    import Utilities
    from copy import deepcopy
    mapval = {tval:tind for tind,tval in enumerate(usetimes)}
    def getCurData(tdata,subtimes):
        """get data for subtimes
        """
        return [[tdata[yind][mapval[tpoint]] for tpoint in subtimes] for yind in xrange(len(tdata))]
    def getRemData(tdata,subtimes):
        """get data for rem times which is in dict form
        """ 
        return [{tpoint:tdata[yind][tpoint] for tpoint in subtimes} for yind in xrange(len(tdata))]
    singledata,singledatadict = geneUtilities.makeSingleData(usedata,usetimes)
    curtimes = Tempselect.initPoints(singledata,usetimes,count,"equal")
    rempoints = list(set(usetimes) - set(curtimes))
    curdata = getCurData(singledata,curtimes)
    remdata = getRemData(usedata,rempoints)   
    toutsplines,ty2knots,tavgerror = Tempselect.splineFitMulti(curdata,remdata,curtimes,rempoints,weights) #initial fitting
    print tavgerror
    avgknot = sum([len(ylen) for ylen in ty2knots])/float(len(ty2knots))
    print avgknot
    print ty2knots[59]
    iterpointsall[0] = list(curtimes)
    y2knotslist[0] = deepcopy(ty2knots)
    outsplineslist[0] = deepcopy(toutsplines) 
    avgerrorlist[0] = tavgerror

    #iterpointsall[0] = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    iterpointsall[0] = [0.5,11.0, 12.5, 14.0,15.0, 20.0, 24.0,28.0]
    #curtimes = iterpointsall[0]  
    #rempoints = list(set(usetimes) - set(curtimes))
    #curdata = getCurData(singledata,curtimes)
    #remdata = getRemData(usedata,rempoints)  
    #toutsplines,ty2knots,tavgerror = Tempselect.splineFitMulti(curdata,remdata,curtimes,rempoints,weights) #initial fitting

        
    inittimes = iterpointsall[0]
    #inittimes = Tempselect.initPoints(yvallist,usetimes,count,"change")
    #print inittimes
    vals1 = [yvallist[gind][time2ind[ttime]] for ttime in inittimes]
    ax2.plot(inittimes, vals1,marker='p',color="b",linestyle='None')
    ax2.set_xlim(0,30)
    
    import scipy.interpolate
    reglambda = 1.0
    curyvals = [yvallist[gind][time2ind[ttime]] for ttime in inittimes]
    outspl = scipy.interpolate.UnivariateSpline(inittimes, curyvals, s=2.0, k=3)
    ax3.plot(usetimes,outspl(usetimes),"r",lw=1)
    ax3.plot(inittimes,curyvals,marker='p',color="b",linestyle='None')
    ax3.plot(usetimes,yvallist[gind],marker='+',color="k",linestyle='None')
    ax2.set_ylim(ax3.get_ylim())
    
    iters = iterpointsall[0:7]
    itcolors = ["","r","b","g","m","y","c"]
    maxiter = len(iters)-1
    for itind in xrange(1,len(iters)):
        curyvals = [yvallist[gind][time2ind[ttime]] for ttime in iters[itind]]
        outspl = outsplineslist[itind][gind] 
        ax6.plot(iters[itind],outspl(iters[itind]),itcolors[itind],lw=1,label="{0}".format(itind+1))
        #ax6.plot(iters[itind],curyvals,marker='p',color=itcolors[itind],linestyle='None',label="{0}".format(itind+1))
    ax6.plot(usetimes,yvallist[gind],marker='+',color="k",linestyle='None',label="Original")    
    ax6.legend(loc=4,prop={'size':7})
    
    curyvals = [yvallist[gind][time2ind[ttime]] for ttime in iters[-1]]
    outspl = outspl = outsplineslist[maxiter][gind]
    ax5.plot(usetimes,outspl(usetimes),"r",lw=1,label="Reconstructed")
    ax5.plot(iters[-1],curyvals,marker='p',color="b",linestyle='None',label="Final points")
    ax5.plot(usetimes,yvallist[gind],marker='+',color="k",linestyle='None',label="Original")    
    ax5.legend(loc=4,prop={'size':7})
               
    for yind,rowdata in enumerate(yvallist):
        if yind in removeset:
           continue 
        curyvals = [rowdata[time2ind[ttime]] for ttime in iters[-1]]
        outspl = outsplineslist[maxiter][yind]
        ax4.plot(iters[-1],outspl(iters[-1]),color=usecolors[yind%len(usecolors)])
        
    #points,y2knots,outsplines,avgsumval = Tempselect.runGreedy(usetimes,usedata,count,weights,"change",1)

    #plt.text(-5.5,1.05,r"Final fit",fontsize = 10)
    #plt.text(0.3,1.05,r"Iterations 2-7",fontsize = 10)
    plt.savefig("algofig_{0}.png".format(gind),dpi=200,format='png')    
    exit(1)
               
    #rect1 = Rectangle((-200,-100), 400, 200, color='white')
    #rect2 = Rectangle((245,245), 150, 150, edgecolor='k',fill=False,facecolor=None)
    ax.add_patch(rect1)
    
    #rax = plt.subplot2grid((8,12),(2,8),colspan=4,rowspan=2)    
    #ax = plt.subplot2grid((8,8),(2, 0),colspan=5,rowspan=2)

    from matplotlib.patches import Rectangle,Circle
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    rax.axes.get_xaxis().set_visible(False)
    rax.axes.get_yaxis().set_visible(False)
    #rax.axis('off')
    rect1 = Rectangle((-200,-100), 400, 200, color='white')
    #rect2 = Rectangle((245,245), 150, 150, edgecolor='k',fill=False,facecolor=None)
    ax.add_patch(rect1)
    
    #rax.add_patch(rect2)
    plt.text(483,320,'Subpopulation 1',fontsize=16)
    plt.text(483,230,'Subpopulation 2',fontsize=16)
    plt.text(483,140,'Subpopulation 3',fontsize=16)
    plt.text(483,50,'Subpopulation 4',fontsize=16)
    plt.xlim([0, 400])
    plt.ylim([0, 400])
    allcolors = ['r','b','y','g']
    from matplotlib.patches import Circle 
    #cir1 = Circle((560,425),10,color='r',clip_on=False)
    #cir2 = Circle((560,393),10,color='g',clip_on=False)
    #cir3 = Circle((560,360),10,color='b',clip_on=False)
    #cir4 = Circle((560,330),10,color='y',clip_on=False)
    #cir1 = Circle((4500,3380),85,color='r',clip_on=False)
    #cir2 = Circle((4500,3130),85,color='g',clip_on=False)
    #cir3 = Circle((4500,2880),85,color='b',clip_on=False)
    #cir4 = Circle((4500,2630),85,color='y',clip_on=False)
    cir1 = Circle((2250,1700),40,color='r',clip_on=False)
    cir2 = Circle((2250,1575),40,color='g',clip_on=False)
    cir3 = Circle((2250,1450),40,color='b',clip_on=False)
    cir4 = Circle((2250,1325),40,color='y',clip_on=False) 
    ax.patches.append(cir1)
    ax.patches.append(cir2)
    ax.patches.append(cir3)
    ax.patches.append(cir4)
    plt.text(20,-85,'Ensemble Cell Population',fontsize=20)
    nodecount = 48
    rowcount,colcount = 12,4
    G = nx.Graph()
    colors = []
    nodes = []
    sizes = []
    node2pos = {}
    for node in xrange(nodecount):
        rowind = node/colcount
        colind = node-(rowind*colcount)
        G.add_node(node)
        nodes.append(node)
        sizes.append(250)
        colors.append(random.choice(allcolors))
        node2pos[node] = (32*rowind+15, 100*colind+40)
    nx.draw_networkx_nodes(G, node2pos, nodelist=nodes, node_size=sizes, node_color=colors, node_shape='o')
    
    mycmap = plt.cm.Reds
    mat1 = np.zeros((20,20),dtype=np.float)
    infos1 = [(0,3),(5,2),(10,3),(16,2)]
    for start,curlen in infos1:
        mat1[start:start+curlen,start:start+curlen] = 1.0
    for ind1,ind2 in itertools.combinations(range(20),2):
        if random.random() < 0.04:
           mat1[ind1,ind2] += 1.0
           mat1[ind2,ind1] += 1.0
    mat2 = np.zeros((20,20),dtype=np.float)
    infos2 = [(2,4),(7,6),(14,2),(17,2)]
    for start,curlen in infos2:
        mat2[start:start+curlen,start:start+curlen] = 1.0
    cnt = 0    
    for ind1,ind2 in itertools.combinations(range(20),2):
        if random.random() < 0.04:
           mat2[ind1,ind2] += 1.0
           mat2[ind2,ind1] += 1.0
           cnt += 1
    tax1 = plt.subplot2grid((30,30), (0, 0),colspan=4,rowspan=5)
    #tax1.imshow(mat1,cmap=mycmap)
    tax1.pcolor(mat1,cmap=mycmap,edgecolors='k',linewidths=0.02,vmin=0, vmax=2)
    tax1.get_xaxis().set_visible(False)
    tax1.get_yaxis().set_visible(False)
    tax2 = plt.subplot2grid((30,30), (0, 16),colspan=4,rowspan=5)
    #tax2.imshow(mat2,cmap=mycmap)
    tax2.pcolor(mat2,cmap=mycmap,edgecolors='k',linewidths=0.02,vmin=0, vmax=2)
    tax2.get_xaxis().set_visible(False)
    tax2.get_yaxis().set_visible(False)
    
    transFigure = fig.transFigure.inverted()
    c1 = transFigure.transform(ax.transData.transform([20,370]))
    #c2 = transFigure.transform(tax1.transData.transform([10,10]))
    from matplotlib.patches import FancyArrow
    curarr = FancyArrow(c1[0],c1[1],0.03,0.09,head_width=0.01,head_length=0.01, fc='k', ec='k',transform=ax.figure.transFigure,clip_on=False)
    ax.add_patch(curarr)
    
    c1 = transFigure.transform(ax.transData.transform([370,370]))
    from matplotlib.patches import FancyArrow 
    curarr = FancyArrow(c1[0],c1[1],0.03,0.09,head_width=0.01,head_length=0.01, fc='k', ec='k',transform=ax.figure.transFigure,clip_on=False)
    ax.add_patch(curarr)

    colsize = 25 
    hax1 = plt.subplot2grid((14,colsize),(8, 1),colspan=3,rowspan=3)
    hax2 = plt.subplot2grid((14,colsize),(8, 6),colspan=3,rowspan=3)
    hax3 = plt.subplot2grid((14,colsize),(8, 11),colspan=3,rowspan=3)
    hax4 = plt.subplot2grid((14,colsize),(8, 16),colspan=3,rowspan=3)
    haxsum = plt.subplot2grid((14,colsize),(8, 21),colspan=4,rowspan=3)
    curlabels = ["Subp. {0}".format(ind+1) for ind in xrange(4)] 
    haxlist = [hax1,hax2,hax3,hax4]
    for ind in tempdata.keys():
        hax = haxlist[ind]
        hax.set_xticklabels([])
        hax.set_yticklabels([])
        hax.set_xticks([])
        hax.set_yticks([])
        im = hax.pcolor(tempdata[ind],cmap=mycmap,edgecolors='k',alpha=0.3,linewidths=0.1) #,vmin=0, vmax=maxval)
        hax.set_xlabel(curlabels[ind],fontsize = 20) 
    haxsum.set_xticklabels([])
    haxsum.set_yticklabels([])
    haxsum.set_xticks([])
    haxsum.set_yticks([])
    haxsum.set_xlabel("Ensemble",fontsize = 20) 
    im = haxsum.pcolor(tempsum,cmap=mycmap,edgecolors='k',linewidths=0.1)
    for hax in haxlist:
        hax.text(3.0,3.0,'?',fontsize=50) #,color='grey')
    hax1.text(-4.0,5.0,r"$\lambda_%s$" %1,fontsize = 25)
    hax2.text(-8.0,5.0,r"$+\lambda_%s$" %2,fontsize = 25)
    hax3.text(-8.0,5.0,r"$+\lambda_%s$" %3,fontsize = 25)
    hax4.text(-8.0,5.0,r"$+\lambda_%s$" %4,fontsize = 25)
    haxsum.text(-5.0,5.0,r"$\cong$",fontsize = 35)

    pax1 = plt.subplot2grid((30,30), (0, 4),colspan=6,rowspan=5)
    pax2 = plt.subplot2grid((30,30), (0, 20),colspan=6,rowspan=5)
    import matplotlib.image as Image
    img = Image.imread('pic2.png') 
    pax1.imshow(img)
    pax1.axis('off')
    img2 = Image.imread('pic1.png') 
    #img2 = Image.imread('pymol2.png')    # Open image as PIL image object
    #rsize = img2.resize((img2.size[0]*1.5,img2.size[1]*1.5)) # Use PIL to resize
    rsizeArr = np.asarray(img2)  # Get array back
    pax2.imshow(rsizeArr)
    pax2.axis('off')

    eax1 = plt.subplot2grid((14,colsize),(12, 5),colspan=3,rowspan=2)
    eax2 = plt.subplot2grid((14,colsize),(12, 11),colspan=3,rowspan=2)
    eax3 = plt.subplot2grid((14,colsize),(12, 17),colspan=3,rowspan=2)
    for mtax in [eax1,eax2,eax3]: 
        mtax.set_xticklabels([])
        mtax.set_yticklabels([])
        mtax.set_xticks([])
        mtax.set_yticks([])
    mat3 = np.zeros((10,10),dtype=np.float)
    infos3 = [(0,3),(5,2),(8,2)]
    for start,curlen in infos3:
        mat3[start:start+curlen,start:start+curlen] = 1.0
    mat4 = np.zeros((10,10),dtype=np.float)
    infos4 = [(0,4),(2,3),(7,2)]
    for start,curlen in infos4:
        mat4[start:start+curlen,start:start+curlen] = 1.0
    infos5 = [(0,4),(4,3),(8,2)]    
    mat5 = np.zeros((10,10),dtype=np.float)
    for start,curlen in infos5:
        approx = int(math.ceil(curlen/2.0))
        for ind1 in xrange(curlen):
            for ind2 in xrange(curlen):
                if abs(ind1-ind2) < approx:
                   mat5[start+ind1,start+ind2] = 1.0
    eax1.pcolor(mat3,cmap=mycmap,edgecolors='k',linewidths=0.1)
    eax2.pcolor(mat4,cmap=mycmap,edgecolors='k',linewidths=0.1)
    eax3.pcolor(mat5,cmap=mycmap,edgecolors='k',linewidths=0.1)
    eax2.text(6.0,1.0,'X',fontsize=30,color='r')
    eax1.text(6.0,1.0,'$\checkmark$',fontsize=30,color='g')
    eax3.text(6.0,1.0,'$\checkmark$',fontsize=30,color='g')
    plt.savefig("decon_ex.png",dpi=400,format='png')

        
spasum,spadata = getSpatial()
tempsum,tempdata,phases = getTemporal()
for comp in tempdata.keys():
    for ind1 in xrange(np.shape(tempdata[comp])[0]):
        for ind2 in xrange(np.shape(tempdata[comp])[0]):
            assert tempdata[comp][ind1,ind2] == tempdata[comp][ind2,ind1]   

#makespaPlot(spadata,spasum)
#maketempPlot(tempdata,tempsum,phases)
#makegenPlot(spadata,spasum)

datapath = "../newdata/savedata.pkl"
fp = open(datapath,"rb")
yvallist = cpickle.load(fp)
usetimes = cpickle.load(fp)
gene2ind = cpickle.load(fp)

explot(yvallist,usetimes,gene2ind)
