"""""
Created to connect coords to the generated population files already.  This specific 
version was meant for the parviglumis populations and is ran twice - once for confidence
intervals and once for all individuals.
"""""
import sys

purePops = open(sys.argv[1]) # either allParvPops.out or pureParvPops.out
zeaAllInfo = open(sys.argv[2]) # ZeaAllInfo.14col.noMML
out = open(sys.argv[3], "w") # either allParvGeo.info or pureParvGeo.info

allInfoLST = {} #key:acc  val:coords        
for line in zeaAllInfo:
    if not line.startswith("#"):
        infoLine = line.split("\t")
        
        name = infoLine[2]
        coords = infoLine[5]+"\t"+infoLine[6]
        
        allInfoLST[name] = coords

mexLST = []
storage = []
count = 0
for pop in purePops:
    if not pop.startswith("#"):
        popLine = pop.split("\t")
        mexLST.append(popLine[0])
        count = count + 1

for  mexPop in mexLST:
    coords = allInfoLST[mexPop]
    outPops = mexPop + "\t" + coords + "\n"
    out.write(outPops)
    
out.close()
zeaAllInfo.close()
purePops.close()