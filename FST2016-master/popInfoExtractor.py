import sys

purePops = open(sys.argv[1])  #pureParvPops.list
zeaAllInfo = open(sys.argv[2]) # ZeaAllInfo.14col.noMML
out = open(sys.argv[3], "w") # parvCoords.out or mexCoords.out
NHCOut = open(sys.argv[4], "w") # NHCParvCoords.out or NHCMexCoords.out
hybrids = open(sys.argv[5], "w") # hybridsCoords.out

storage = []
allInfoLST = {} #key:acc  val:coords        
for line in zeaAllInfo:
    if not line.startswith("#"):
        infoLine = line.split("\t")
        
        name = infoLine[2]
        coords = infoLine[5]+"\t"+infoLine[6]
        type = infoLine[1]
        allInfoLST[name] = coords
        storage.append(infoLine)      
        
parvLST = []
count = 0
#For pure parv pops
for pop in purePops:
    if not pop.startswith("#"): #add "and count != 8" if Parv
        popLine = pop.split("\t")
        parvLST.append(popLine[0])
        count = count + 1
        
NCHLST = []
for parvPop in parvLST:
    coords = allInfoLST[parvPop]
    outPops = parvPop + "\t" + coords + "\n"
    out.write(outPops)
    
for info in storage:
    #For either mexicana or parviglumis pops depending on the == value
    if info[1] == "Zea_mays_mexicana": #change this to Zea_Parv..... 
        name = info[2]
        coords = info[5] + "\t"+info[6]
        #If not equal parvPop then it is an ambig indiv
        if not name == parvPop:
            if not name in NCHLST:
                NCHLST.append(name)
                dataLine = name +"\t"+coords+"\n"
                NHCOut.write(dataLine)
    #For getting hybrid populations
    if info[3] == "hybrid":
        name = info[2]
        coords = info[5] + "\t"+info[6]
        hybridData = name +"\t"+coords+"\n"
        hybrids.write(hybridData)
        
hybrids.close() 
NHCOut.close()
out.close()
zeaAllInfo.close()
purePops.close()