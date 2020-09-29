#expression utilities
import numpy as np
import math
import scipy as sp
import random

def convertFormNOTUSED(curtimes):
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

def makeSingleData(usedata,usetimes):
    """makes single data
    Args:
       usedata:
       usetimes:
    Returns:
       singledata,singledatadict:
    """
    singledata, singledatadict = [],[]
    for cind,cdata in enumerate(usedata):
        cavgdata = [np.median(cdata[ttime]) for tind,ttime in enumerate(usetimes)]
        singledata.append(list(cavgdata))
        singledatadict.append({tval:cavgdata[tind] for tind,tval in enumerate(usetimes)})
    return singledata,singledatadict


def shiftData(rawusedata,mintime):
    """shift data if not done already
    Args:
       rawusedata:
       mintime:
    Returns:
       modusedata:
    """
    modusedata = []
    for cind,cdata in enumerate(rawusedata):
        begval = np.median(cdata[mintime])
        shiftdict = {ttime:[tval-begval for tval in cdata[ttime]] for ttime in cdata.keys()}
        modusedata.append(dict(shiftdict))
    return modusedata
