
import sys

def transformData(fpath,outfpath,exppath):
    """transforms data
    Args:
       fpath,outfpath:
       exppath:
    """
    row2data = {}
    count = 0
    field2time = {}
    fields = None
    with open(fpath,"r") as infile:
        for line in infile:
            line = line.rstrip()
            parts = line.split(",")
            if count == 0:
               fields = list(parts[9:])
               for tpart in parts[9:]:
                   field2time[tpart] = float(tpart.split(" ")[0].split("e")[0].split("d")[0].split("c")[0].replace("P",""))
               count += 1
            elif count == 1:
               count += 1 
               continue    
            else:
               tpart = parts[1].split("(")[0].replace(" ","")
               assert len(tpart) == 6 
               assert not row2data.has_key(tpart)
               row2data[tpart] = [float(titem) for titem in parts[9:]]
    
    with open(exppath,"w") as outfile:
        outfile.write("Experiment Names\ttime inclusive\n")
        for tfield,ttime in field2time.items():
            outfile.write("{0}\t{1}\n".format(tfield,ttime))
    with open(outfpath,"w") as outfile:
        outparts = ["UniqueID"] + fields
        outfile.write("\t".join(outparts)+"\n")
        for tgene in row2data.keys():
            outparts = [tgene] + [str(titem) for titem in row2data[tgene]]
            outfile.write("\t".join(outparts)+"\n")
    exit(1)


datapath = "proteomedata.csv"
exppath = "expdes.txt"
outdatapath = "processed_proteome.txt"
transformData(datapath,outdatapath,exppath)
exit(1)

if False:
 parts = []
 fpath = "gennames.txt"
 with open(fpath,"r") as infile:
    for line in infile:
        line = line.rstrip()
        part = line.split("(")[0].replace(" ","")
        assert len(part) == 6
        parts.append(part)
 print len(parts)
 print len(set(parts))        

 outfpath = "genout.txt"
 with open(outfpath,"w") as outfile:
    for tpart in parts:
        outfile.write("{0}\n".format(tpart))

parts = []
outfpath = "genout.txt"
with open(outfpath,"r") as infile:
    for line in infile:
        line = line.rstrip()
        parts.append(line)
print len(parts)   

id2gene = {}        
mapfpath = "mapped.txt"
with open(mapfpath,"r") as infile:
    for line in infile:
        uid,gene = line.rstrip().split()    
        id2gene.setdefault(uid,set())
        id2gene[uid].add(gene)

print set(parts) - set(id2gene.keys())            
for tid in id2gene.keys():
    if len(id2gene[tid]) != 1:
       print tid,id2gene[tid] 
