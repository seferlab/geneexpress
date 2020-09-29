if False:
     initoutpoints = list(set(times) - set([times[0],times[-1]]))
     selpoints = [sorttimes[0]] + random.sample(initoutpoints, count-2) + [sorttimes[-1]]
     tyvals = [x2val[0][spoint] for spoint in selpoints]
     trempoints = list(set(sorttimes) - set(selpoints))
     #tck,fp,ier = scipy.interpolate.splrep(points, yvals, s=reglambda)
     spl = scipy.interpolate.UnivariateSpline(selpoints, tyvals, s=reglambda, k=3)
     #print spl.__dict__.keys()
     #print spl._data
     while False:
        tind += 1
        if tind % 10000 == 0:
           print tind 
        initoutpoints = list(set(sorttimes) - set([sorttimes[0],sorttimes[-1]]))
        selpoints = [sorttimes[0]] + random.sample(initoutpoints, count-2) + [sorttimes[-1]]
        tyvals = [x2val[spoint] for spoint in selpoints]
        trempoints = list(set(sorttimes) - set(selpoints))
        #tck,fp,ier = scipy.interpolate.splrep(points, yvals, s=reglambda)
        spl = scipy.interpolate.UnivariateSpline(selpoints, tyvals, s=reglambda, k=3)
        sumval = estQual(spl,x2val,selpoints,trempoints)
        if sumval <= minval:
           minval = sumval
           print minval
        #print sumval
        continue
        
        exit(1)
        minsol = None
        for addpoint in rempoints:
            for delpoint in points:
                if delpoint in [sorttimes[0],sorttimes[-1]]:
                   continue 
                newrempoints, newpoints = list(rempoints), list(points)
                newrempoints.remove(addpoint)
                newrempoints.append(delpoint)
                newpoints.remove(delpoint)
                newpoints.append(addpoint)
                inyvals = [x2val[rpoint] for rpoint in newpoints]
                inspl = scipy.interpolate.UnivariateSpline(newpoints, inyvals, s=reglambda,k=3)
                cursumval = estQual(inspl,x2val,newpoints,newrempoints)
                if cursumval < minval:
                   minval = cursumval
                   minsol = (addpoint,delpoint)
        if minsol == None:
           break 
        rempoints.remove(minsol[0])
        rempoints.append(minsol[1])
        points.remove(minsol[1])
        points.append(minsol[0])
        yvals = [x2val[tpoint] for tpoint in points]
