import sys
import os
import itertools
import xlrd
import random
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model

def estCoefs(times,stdlist):
    """estimates coefs
    Args:
       times:
       stdlist:
    Returns:
       usecoefs:
    """
    usecoefs = []
    for tind1,time1 in enumerate(times):
        curcoef = []
        for tind2,time2 in enumerate(times):
            diftime = abs(time2-time1)
            curcoef.append(math.exp(-0.5*(diftime**2)/stdlist[tind2]))
        usecoefs.append(list(curcoef))
    return np.array(usecoefs)

def runGaussTwoStep(times,ydata,pointnum):
    """run gaussian mixture two step procedure
    Args:
       times:
       ydata:
       pointnum:
    Returns:
    """
    meant = np.mean(times)
    stdlist = [2.0*meant] * len(times)
    usecoefs = estCoefs(times,stdlist)
    optobj = 10000000000.0
    optweight,optstds = None,None
    while True:
        alphas, _, coefs = linear_model.lars_path(usecoefs, np.array(ydata), method='lasso', verbose=True)
        alpha2count = {}
        for tind in xrange(np.shape(coefs)[1]):
            zerocount = len([item for item in coefs[:,tind] if abs(item) > 0.0000000000000001])
            alpha2count[alphas[tind]] = zerocount
        print alpha2count
        for alpha1,alpha2 in itertools.combinations(alpha2count.keys(),2):
            if alpha1 < alpha2:
               if alpha2count[alpha1] < alpha2count[alpha2]:
                  print alpha1
                  print coefs[:,list(alphas).index(alpha1)]
                  print alpha2
                  print coefs[:,list(alphas).index(alpha2)]
                  exit(1)
        exit(1)
            
        print coefs
        print np.shape(coefs)
        print len(alphas)
        print len(times)
        
        sortindices =  sorted(range(len(alphas)), key=lambda k: alphas[k])
        sortalphas = [alphas[tind] for tind in sortindices]
        sentcoefs = np.transpose(np.array([coefs[:,tind] for tind in sortindices]))
        print np.shape(sentcoefs)
        for tind in xrange(np.shape(sentcoefs)[1]):
            print tind,sortalphas[tind],len([item for item in sentcoefs[:,tind] if abs(item) > 0.0000000000000001])
            if len([item for item in sentcoefs[:,tind] if abs(item) > 0.0000000000000000001]) == 0:
               print sentcoefs[:,tind]
               #exit(1)
        print len(times)    
        exit(1)
        
        print alphas
        print sortalphas
        print sortindices
        exit(1)
        sortcoefs = []
        #alphas,coefs = estPath(usecoefs,ydata)
        print alphas
        print coefs
        print alphas
        exit(1)
        print np.shape(coefs)
        for tind in xrange(np.shape(coefs)[0]):
            print tind,len([item for item in coefs[tind,:] if item!=0])
            print coefs[tind,:]
            #exit(1)
        exit(1)
        #estStd()
    return

def estPath(values):
    """estimates path
    Args:
       values: dict of x and y
    Returns:
       alphas: regularization lambdas
       coefs: coef matrix for features and alphas
    """
    X,y = values["x"], values["y"]
    alphas, _, coefs = linear_model.lars_path(X, y, method='lasso', verbose=True)
    return alphas,coefs

    print alphas
    print coefs
    print coefs[:,1]
    exit(1)
    xx = np.sum(np.abs(coefs.T), axis=1)
    xx /= xx[-1]

    plt.plot(xx, coefs.T)
    ymin, ymax = plt.ylim()
    plt.vlines(xx, ymin, ymax, linestyle='dashed')
    plt.xlabel('|coef| / max|coef|')
    plt.ylabel('Coefficients')
    plt.title('LASSO Path')
    plt.axis('tight')
    plt.show()
    plt.savefid("larspath.png")


#values = {}
#usevals = [1,3,5,6,8,13,4,3.2]
#values["x"] = np.array([usevals[tind:tind+3] for tind in xrange(len(usevals)-2)])
#values["y"] = np.array([-3,4,9,-2,13,-4])
#times = [0.5,1.0,1.5,2,3,4,5,6,6.5]
#ydata = [2.0*random.uniform(-1,1) for elem in times]
#pointnum = 3
#runGaussTwoStep(times,ydata,pointnum)
#exit(1)

#stdlist = [1.0] * len(times)
#runGauss(values)
#estCoefs(times,ydata,stdlist)
#cexit(1)

def writeDataFile(fpath,times,gene2time):
    """writes data file
    Args:
       fpath:
       times,gene2time:
    """
    with open(fpath,"w") as outfile:
        outfile.write(",".join([str(time) for time in times])+"\n")
        for tgene in gene2time.keys(): 
            outfile.write(tgene+","+",".join([str(tval) for tval in gene2time[tgene]])+"\n")

fname = "42 time points LCM nanostring log notmalized.xlsx"
workbook = xlrd.open_workbook(fname)
print workbook.sheet_names()
worksheet = workbook.sheet_by_name('RAW DATA')
num_rows = worksheet.nrows - 1
num_cells = worksheet.ncols - 1
labels = [worksheet.cell_value(0, tind) for tind in xrange(2,num_cells-1)]
print labels
print len(labels)
setlabels = set(labels)
print len(setlabels)
multis = []
for label in setlabels:
    curcount = 0
    for tlabel in labels:
        if label == tlabel:
           curcount += 1
    if curcount > 1:
       multis.append(label)
multi2ind = {}
for toplabel in multis:
    indset = set()
    for tind,tlabel in enumerate(labels):
        if toplabel == tlabel:
           indset.add(tind)
    assert len(indset) > 1       
    multi2ind[toplabel] = indset
curr_row = 0
gene2time = {}
varsum = 0.0
varcount = 0
while curr_row < num_rows:
    curr_row += 1
    genevals = [worksheet.cell_value(curr_row,tind) for tind in xrange(2,num_cells-1)]
    logvals = [math.log(tval,2.0) for tval in genevals]
    for mkey in multi2ind.keys():
        cvals = [logvals[sind] for sind in multi2ind[mkey]]
        varsum += np.std(cvals)**2
        varcount += 1
    #gene = str(worksheet.cell_value(curr_row,0))
    #gene2time[gene] = list(genevals)
print varsum
print varcount
print varsum/varcount
#exit(1)
                
fname = "42 time points LCM nanostring log notmalized.xlsx"
workbook = xlrd.open_workbook(fname)
print workbook.sheet_names()
worksheet = workbook.sheet_by_name('Sheet2')
num_rows = worksheet.nrows - 1
num_cells = worksheet.ncols
labels = [worksheet.cell_value(0, tind) for tind in xrange(2,num_cells)]
similars, seens = [], []
for tind,tlabel in enumerate(labels):
    curval = float(tlabel.replace("-a","").replace("P",""))
    if curval not in seens:
       seens.append(curval)
    else:
       similars.append(tind) 
times = sorted(list(set([float(label.replace("-a","").replace("P","")) for label in labels])))
curr_row = 0
gene2time = {}
while curr_row < num_rows:
    curr_row += 1
    genevals = [worksheet.cell_value(curr_row,tind) for tind in xrange(2,num_cells) if tind-2 not in similars]
    gene = str(worksheet.cell_value(curr_row,0))
    gene2time[gene] = list(genevals)
fpath = "input.data"    
writeDataFile(fpath,times,gene2time)
exit(1)
    
worksheet = workbook.sheet_by_name('Sheet1')
num_rows = worksheet.nrows - 1
num_cells = worksheet.ncols - 1
labels = [worksheet.cell_value(0, tind) for tind in xrange(1,num_cells+1)]
curr_row = 0
tgene2time = {}
while curr_row < num_rows:
    curr_row += 1
    genevals = [worksheet.cell_value(curr_row,tind) for tind in xrange(1,num_cells+1)]
    gene = str(worksheet.cell_value(curr_row,0))
    tgene2time[gene] = list(genevals)        
print len(tgene2time.values()[0])
exit(1)

print gene2time.values()[2]
print labels
labels = [float(tlabel.replace("P","")) if tlabel.find("-a")==-1 else float(tlabel.replace("P","").replace("-a",""))-0.25 for tlabel in labels]
times = [float(tlabel) for tlabel in labels]
for gene in gene2time.keys():
    print len(gene2time[gene])
fpath = "input.data"    
writeDataFile(fpath,times,gene2time)
exit(1)

minval,maxval = 1000000, -1000000
for gene in gene2time.keys():
    for tval in gene2time[gene]:
        if tval < minval:
           minval = tval
        if tval > maxval:
           maxval = tval
print minval
print maxval
print len(gene2time.keys())    

