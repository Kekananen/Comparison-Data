#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 10:39:37 2019
This script finds the amount of homozygous and heterozygous biallelic counts 
in three populations within the hapmap file (PT, T43, and PTxT43). Then, it 
looks for the segragation within the population of PTxT43 with regards to the 
other two. This is achieved by looked though all variation in the PT and T43 
population as they consist of 5 individuals at a max if none are missing and
comparing each variation possibilty to the highest one in the PTxT43 until 
there is a match.  If there are no matches, it will be predicted based on the 
pair that is closest to the frequency of 12.5% in the population.

Flags will be added if the highest frequency in PTxT43 matches to certain 
cases outlined in comments to the rqtlOut file.

@author: Kathryn Kananen

Note: I have checked the SNPs_HAPMAP_DMz18-205.csv for duplicate IDs and 
there appears to be none - Kathryn 
"""
import sys, re

file = open(sys.argv[1], "r") #SNPs_HAPMAP_DMz18-205.csv    v3
BLASTFile = open(sys.argv[2], "r") #SNPs_BLASTv.4_TwoRows_DMz18-205.csv
rqtlOut = open(sys.argv[3], "w") #SNPs_HAPMAP_DMz18-205v2.csv
parentOut = open(sys.argv[4], "w") #parentsID.csv
PTorT43Out = open(sys.argv[5], "w") #SNPs_HAPMAP_DMz18-205v3.csv

"""""
Finds all the variation in a string of bialleles and returns a condensed list
of the various types. For example: ["AA, 4", "AC, 2"].

    Takes: a string of comma seperated bialleles and the locations
    Returns: a list of lists
"""""
def VariationGetter(alleleString, indexLST):
    alleleLST = alleleString.split(",")
    temp = []; outLST = []
    
    for index in indexLST:
        temp.append(alleleLST[index])
    outLST.append([[allele, temp.count(allele)] for allele in set(temp)])
    
    return outLST

"""""
Finds frequency of the types of bases and returns them as a list
    Takes: the output of the VariationGetter() function, string of seq
    Returns: a list frequencies in order of A T C G
"""""
def VariationFeqFinder(givenLST, string):
    #Get's rid of the Missing data from the calculations
    ACnt = 0; CCnt = 0; TCnt = 0; GCnt = 0;
    
    outLST = []; length = string
    for pair in givenLST[0]: 
        #Remove missing data
        if pair[0] == "NN":
            length = length - 2*pair[1]
        #Checks for homozygous pairs
        elif pair[0] == "CC":
            CCnt = CCnt + pair[1]*2 
        elif pair[0] == "GG":
            GCnt = GCnt + pair[1]*2 
        elif pair[0] == "AA":
            ACnt = ACnt + pair[1]*2 
        elif pair[0] == "TT":
            TCnt = TCnt + pair[1]*2 
        #Checks for heterozygous pairs
        elif pair[0][0] != pair[0][1]:
            #Section for A_
            if pair[0][0] == "A":
                ACnt = ACnt + pair[1]
            if pair[0][0] == "C":
                CCnt = CCnt + pair[1]
            if pair[0][0] == "G":
                GCnt = GCnt + pair[1]
            if pair[0][0] == "T":
                TCnt = TCnt + pair[1]
            #Section for _A
            if pair[0][1] == "A":
                ACnt = ACnt + pair[1]
            if pair[0][1] == "C":
                CCnt = CCnt + pair[1]
            if pair[0][1] == "G":
                GCnt = GCnt + pair[1]
            if pair[0][1] == "T":
                TCnt = TCnt + pair[1]
   
    if not length == 0:            
        outLST.append(ACnt/length); outLST.append(TCnt/length)  
        outLST.append(CCnt/length); outLST.append(GCnt/length)
    else:
        outLST.append(0); outLST.append(0)  
        outLST.append(0); outLST.append(0)

    return outLST

"""""
Makes a string of the all the variation and the frequencies of the variations
in the PTxT43 population.
    Takes: the output of the VariationGetter() and VariationFeqFinder()
    Returns: A string of combined information
"""""
def outLineBuilder(givenLST, givenFeqLST):
    outLine = ""
    homoCnt = 0; hetCnt = 0
    for pair in givenLST:
        for i in range(len(givenLST[0])):
            if not pair[i-1][0] == "NN":
                outLine = outLine +str(pair[i-1][1])+ pair[i-1][0]
                
                if pair[i-1][0][0] ==  pair[i-1][0][1]:
                    homoCnt = homoCnt + pair[i-1][1]           
                else:
                    hetCnt = hetCnt + pair[i-1][1]
                    
    outLine = outLine + "," + str(homoCnt) + "," + str(hetCnt) + ","
    
    for freq in givenFeqLST:
        outLine = outLine + str(freq) + ","
        
    return outLine[:len(outLine)-1]

"""""
Takes in a string of information from the outLineBuilder() and then finds the
most likely allele based on the other alleles found around it. 
    Takes: the output of outLineBuilder() and the snp
    Returns: probability of allele pair and the allele pair.
"""""
def MostLikelyAllele(givenString, mode):
    allInfoLST = givenString.strip().split(",")    
    #First find the pair is homozygous or heterozygous    
    posTypeLST = list(filter(None, re.split(r'(\d+)', allInfoLST[0]))) 
    popPercentDict = {}
    #get size of SNP
    sizeLST = re.findall(r'\d+', allInfoLST[0])
    length = 0
    for size in sizeLST:
        length = int(size) + length
    length = length*2    
    
    #Check for easy case first of all one type regardless of if its a het or
    #a homozygous pair.
    likelyAlle = []
    if len(posTypeLST) == 2:
        #This keeps missing data in the count.
        prob = (int(posTypeLST[0])*2)/length
        likelyAlle.append(posTypeLST[1] +","+ str(prob))
        popPercentDict[posTypeLST[1]] = prob
    elif len(posTypeLST) == 0:
        likelyAlle.append("NN,0")
        popPercentDict["NN"] = 0 
    #Check for multiple different types of alleles pairs in a snp and finds 
    #the probability that the allele is that allele pair in the snp
    else:
        cnt = 0; 
        for i in range(len(posTypeLST)):
            pair = "";          
            if not cnt%2==0:               
                pair = posTypeLST[i]                
                freq = (int(posTypeLST[i-1])*2)/length
                
                popPercentDict[pair] = freq 
        
            cnt+=1 
            
        #Finds the most likely allele pair
        maxfreq = max(popPercentDict.values())
        maxPair = max(popPercentDict, key=lambda k: popPercentDict[k])            
        likelyAlle.append(maxPair +","+ str(maxfreq))

    if mode == 0:
        return likelyAlle
    else:
        return popPercentDict
    
"""""
Decides if snp is homozygous or heterozygous depending on variation amount.
    Takes: The output from MostLikelyAllele(), mode 2.
    Returns: the following codes if the statement is true: 
        0 if all homozygous 
        1 if all heterozygous
        2 if combination of one homozygous and heterozygous type
        3 if combination of two or more homozygous types
        4 if combination of two or more heterozygous types
        5 if combination of case 3 and 4
        6 if none of the cases above
"""""
def HomoVHetFinder(variationDict):
    allelesLST = list(variationDict.keys())
    homoCnt = 0; hetCnt = 0
    if len(allelesLST) > 1:
        for i in range(len(allelesLST)):
            if allelesLST[i][0] == allelesLST[i][1]:
                homoCnt +=1
            elif allelesLST[i][0] != allelesLST[i][1]:
                hetCnt +=1
                      
        if homoCnt > 1 and hetCnt == 0:
            return 3
        elif homoCnt == 1 and hetCnt == 1:
            return 2
        elif (homoCnt > 1 and hetCnt >=1) or (homoCnt >= 1 and hetCnt >1):
            return 5
        elif homoCnt == 0 and hetCnt >= 1:
            return 4
        else:
            return 6        
    else:
        if allelesLST[0][0] == allelesLST[0][1]:
            return 0
        else:
            return 1
""""
Looks through the a list of the childs variation and sees if the parents match
with any of the pairs of the child. This continues until a match is found or 
the child is out of pairs in which the frequencies are looked through and the 
one that is closest to 12.5% is frequency is chosen to be the pair to be used 
for the function.
    Takes: A dictionary with pairs of alleles and their frequencies and two 
        dictionaries of just parent values.
    Returns: A string such as "parent1,donated allele,parent2,donated allele" 
"""""
def ParentFileBuilderHelper(childPairs,parent1Pairs,parent2Pairs):
    
    #The alleles and frequencies of the child 
    childAllelesLST = list(childPairs.keys())
    childFrequenceLST = list(childPairs.values())
    PTAllelesLST = list(parent1Pairs.keys())
    T43AllelesLST = list(parent2Pairs.keys())
    
    posPTPair = []; posT43Pair = []
    posPTFreq = []; posT43Freq = []
        
    outLine=""
    #look though for matches in PT population
    if len(posPTPair) == 0:
        #The PT is likely to be closest freq to 12.5%
        tempPTFreq = str(min(childFrequenceLST, key=lambda x: x-.125)) 
        cnt = 0 #While loop counter
        for i in range(len(childFrequenceLST)):
            base1 = childAllelesLST[i][0]
            base2 = childAllelesLST[i][1]
            #If frequencies match each others 
            if str(childFrequenceLST[i]) == str(tempPTFreq):
                while cnt <= len(childAllelesLST):                        
                    for PTPair in PTAllelesLST:
                        #Case 1: The child matches one of the PT
                        if((base1 == PTPair[0] or base1 == PTPair[1] or 
                           base2 == PTPair[0] or base2 == PTPair[1]) and 
                            len(posPTPair) == 0): 
                            posPTPair.append(childAllelesLST[i])
                            posPTFreq.append(childFrequenceLST[i])
                            break
                        break
                    cnt+=1    
                else:
                    if len(posPTPair) == 0:
                        posPTFreq.append(str(childFrequenceLST[i])+"*")
                        posPTPair.append(childAllelesLST[i])
                                                           
        #look though for matches in T43 population
        if len(posT43Pair) == 0:
            #The T43 is likely to be highest frequency 
            tempT43Freq = max(childFrequenceLST)
            cnt = 0 #While loop counter
            for i in range(len(childFrequenceLST)):
                base1 = childAllelesLST[i][0]
                base2 = childAllelesLST[i][1]
                #If frequencies match each others 
                if str(childFrequenceLST[i]) == str(tempT43Freq):
                    while cnt <= len(childAllelesLST): 
                        for T43Pair in T43AllelesLST:
                            #Case 1: The child matches one of the T43
                            if((base1 == T43Pair[0] or base1 == T43Pair[1] or 
                               base2 == T43Pair[0] or base2 == T43Pair[1]) and 
                                len(posT43Pair) == 0):                                 
                                posT43Pair.append(childAllelesLST[i])
                                posT43Freq.append(childFrequenceLST[i])
                                break
                            break
                        cnt+=1    
                    else:
                        if len(posT43Pair) == 0:
                            posT43Freq.append(str(childFrequenceLST[i])+"*")
                            posT43Pair.append(childAllelesLST[i])
    T43var = 0
    PTvar = 0                                            
    #Homozygous T43 pair
    if posT43Pair[0][0] == posT43Pair[0][1]:        
        for key in childPairs:
            if not str(posT43Freq[0]).endswith("*"):
                if key[0] == posT43Pair[0][0] and key[1] == posT43Pair[0][1]:
                    T43var = T43var + childPairs[key]
                elif((key[0] == posT43Pair[0][0] and 
                      not key[0] == posT43Pair[0][0] 
                      and not key[1] == posT43Pair[0][1]) or (
                              key[1] == posT43Pair[0][1] and 
                              not key[0] == posT43Pair[0][0] 
                              and not key[1] == posT43Pair[0][1])):
                    T43var = T43var + (childPairs[key]/2)
            else:
                T43var = posT43Freq[0]
                
            if not str(posPTFreq[0]).endswith("*"):
                if key[0] == posPTPair[0][0] and key[1] == posPTPair[0][1]:
                    PTvar = PTvar + childPairs[key]
                elif((key[0] == posPTPair[0][0] and 
                      not key[0] == posPTPair[0][0] 
                      and not key[1] == posPTPair[0][1]) or (
                              key[1] == posPTPair[0][1] and 
                              not key[0] == posPTPair[0][0] 
                              and not key[1] == posPTPair[0][1])): 
                    PTvar = PTvar + (childPairs[key]/2)
            else:
                PTvar = posPTFreq[0]
                
        # Heterozygous PT pair
        if not posPTPair[0][0] == posPTPair[0][1]:
            #T43 base1 matches PT base1
            if posT43Pair[0][0] == posPTPair[0][0]:
                outLine =(posPTPair[0][1] + "," + str(PTvar) 
                + "," + posT43Pair[0][0] + "," + str(T43var))
            #T43 base1 matches PT base2
            elif posT43Pair[0][0] == posPTPair[0][1]:
                outLine =(posPTPair[0][0] + "," + str(PTvar) 
                + "," + posT43Pair[0][0] + "," + str(T43var)) 
        else:
                outLine =(posPTPair[0][1] + "," + str(PTvar) 
                + "," + posT43Pair[0][0] + "," + str(T43var))                
    else:
        #Homozygous PT pair
        if posPTPair[0][0] == posPTPair[0][1]:
            #T43 base1 matches PT base1
            if posT43Pair[0][0] == posPTPair[0][0]:
                outLine =(posPTPair[0][1] + "," + str(PTvar) 
                + "," + posT43Pair[0][0] + "," + str(T43var))
            elif posT43Pair[0][1] == posPTPair[0][0]:  
                 outLine =(posPTPair[0][0] + "," + str(PTvar) 
                + "," + posT43Pair[0][0] + "," + str(T43var))
        else:
            print("Build this case into method") 
    
    return outLine
        

######## This section is the end of the functions ########
        
lineLST = [] #For output of the new HapMap file later
rowLST = [] #Names of the rows and populations
alleleLST = [] #The list of pairs of alleles

IDLST = []

IDDict = {} #The Allele ID's and the sequence information
IDDictInfo = {}

cnt = 0;
for line in file:    
    allInfoLST = line.strip().split(",")
    #This sections gathers up the row information
    if cnt == 4:
        for i in range(len(allInfoLST)):
            if not allInfoLST[i] == "*": 
                rowLST.append(allInfoLST[i].split(" ")[0])                
    #This section extracts the alleles and the SNP IDs
    if cnt > 5:
        alleleLine = ""
        
        IDDictInfo[str(allInfoLST[0])] = ""
        IDLST.append(allInfoLST[0]) 
        
        for i in range(len(allInfoLST)):
            if i > 20:
                alleleLine = alleleLine + allInfoLST[i] + ","
        #Grabbing the entire ID information instead of just the identifyer
        IDDict[allInfoLST[0]] = alleleLine[:len(alleleLine)-1] #Removes comma
        alleleLine = ""      
        
    cnt +=1
    lineLST.append(line)

for line in BLASTFile:
    if not line.startswith("*"):
        allInfo = line.strip().split(",")

        ID = allInfo[0]; chrom = allInfo[5]
        if(not chrom == "" and 
           not chrom[0] == "P" and 
           not chrom[0] == "B" and 
           not chrom[0] == "M"):
            
            chromNum = ""
            for letter in chrom:
                if letter.isdigit():
                    chromNum = chromNum + letter
                    
            chrom = chromNum            
        else:
            chrom = "0"
            
        #Since the two files are not the same line number wise, we need to 
        #match them to make sure the chromosomes are the correct ones.
        if ID in IDLST:
            IDDictInfo[ID] = allInfo[6] + "," + chrom

chromDict = {}
posDict = {}
for key in IDDictInfo:
    chromDict[key] = IDDictInfo[key].split(",")[1]
    posDict[key] = IDDictInfo[key].split(",")[0]
    
   
#Indeces of the different populations 
T43Index = [i for i,allele in enumerate(rowLST) if allele=="T43"]
PTIndex = [i for i,allele in enumerate(rowLST) if allele=="PT"]
PTxT43Index = [i for i,allele in enumerate(rowLST) if allele=="PTxT43"]

titleLine = ("AlleleID,Variation,Homo,Hetero," +
             "FreqA,FreqT,FreqC,FreqG,highFeqPair,freq,\n")
parentOut.write("ID,PTFreq,PTType,T43Freq,T43Type\n")

PTxT43outDict = {}; PTxT43VarDict = {}; PTxT43FlagDict = {}
T43outDict = {}; parentDict = {}
PToutDict = {} 
PTinPTxT43Dict = {}
#Starts to analyze the population homo vs. hetero.
for key in IDDict:   
    #PTPTxT43 Section and information processing
    PTxT43Variation = VariationGetter(IDDict[key], PTxT43Index) 
    PTxT43AlleFeq = VariationFeqFinder(PTxT43Variation, len(PTxT43Index)*2)
    PTxT43InfoLine = outLineBuilder(PTxT43Variation, PTxT43AlleFeq)
    PTxT43Alle = MostLikelyAllele(PTxT43InfoLine, 0)
    PTxT43VarDict[key] = PTxT43InfoLine.split(",")[0]
 
    #PT Section and information processing
    PTVariation = VariationGetter(IDDict[key], PTIndex)
    PTAlleFeq = VariationFeqFinder(PTVariation, len(PTIndex)*2)
    PTInfoLine = outLineBuilder(PTVariation, PTAlleFeq)
    PTAlle = MostLikelyAllele(outLineBuilder(PTVariation, PTAlleFeq), 0) 

    #T43 Section and information processing
    T43Variation = VariationGetter(IDDict[key], T43Index)
    T43AlleFeq = VariationFeqFinder(T43Variation, len(T43Index)*2)
    T43InfoLine = outLineBuilder(T43Variation, T43AlleFeq)
    T43Alle = MostLikelyAllele(outLineBuilder(T43Variation, T43AlleFeq), 0)
                
    #Information for out file
    PTxT43flag = HomoVHetFinder(MostLikelyAllele(PTxT43InfoLine, 1))
    PTxT43outDict[key] = PTxT43flag
    
    PTflag = HomoVHetFinder(MostLikelyAllele(PTInfoLine, 1))
    PToutDict[key] = PTflag
    
    T43flag = HomoVHetFinder(MostLikelyAllele(T43InfoLine, 1))
    T43outDict[key] = T43flag

    """"" --------------------------------------------------------------------
    This section is for a bunch of conditionals that creat the flags for the
    output file.  Each has been tested using a counter that no longer exists
    in this file and has the number of individuals in comments next to the 
    case for each conditional statement. This section relates directly with 
    the SNPs_HapMap_DMz18-205.csv only.  Disregard numbers if this file is not 
    being used. - Kathryn
    
    The indivs are no longer up to date that are next to them in the comments
    """""
    PTxT43 = PTxT43Alle[0].split(",")[0]
    T43 = T43Alle[0].split(",")[0]
    PT = PTAlle[0].split(",")[0]
    
    #Holds all the individuals of the T43 and PT for a single SNP.
    T43Dict = MostLikelyAllele(T43InfoLine, 1)
    posT43 = T43Dict.keys()
    PTDict = MostLikelyAllele(PTInfoLine, 1)
    posPT = PTDict.keys()
    
    PTxT43Dict = MostLikelyAllele(PTxT43InfoLine, 1)
    
    PTinPTxT43Dict[key]=ParentFileBuilderHelper(PTxT43Dict,PTDict,T43Dict).split(",")[1]
    parentDict[key]=ParentFileBuilderHelper(PTxT43Dict,PTDict,T43Dict)
    
    #These are the cases where the predicted parent alleles match or the PT 
    #wasn't found (there were no samples in the population so it was predicted)

    #See if the PT allele is in the PTxT43 population SNP        
    PTPredicted =(parentDict[key].split(",")[0] + parentDict[key].split(",")[0])
    T43Predicted =(parentDict[key].split(",")[2] + parentDict[key].split(",")[2])
    HetPredicted1 = (parentDict[key].split(",")[2] + parentDict[key].split(",")[0])
    HetPredicted2 = (parentDict[key].split(",")[0] + parentDict[key].split(",")[2])

    #This is just the variation in the PTxT43 population
    NNCnt = 0;
    hetCnt = 0
    homoCnt = 0
    for firstLevel in PTxT43Variation:
        for variation in firstLevel:
            #remove the NN data from the CNT
            if variation[0] == "NN":
                NNCnt = int(variation[1])
            elif variation[0][0] == variation[0][1]:
                homoCnt = homoCnt + int(variation[1])                
            elif variation[0][1] != variation[0][0]:
                hetCnt = hetCnt + int(variation[1])
            else:
                print("ERROR 5")

#       #removes odd case without PT in PTxT43 population
#       #Only look at those individuals in the PTxT43 population
        noneParentLine = (IDDict[key].split(",")[:len(IDDict[key].split(","))-11]
                [:len(IDDict[key].split(","))-10])
        del noneParentLine[len(noneParentLine)-2]
        noneParentLine = ",".join(noneParentLine)

        #the allele from the T43 population
        majorCnt = noneParentLine.split(",").count(T43Predicted)
        #the allele from the PT population        
        minorCnt = noneParentLine.split(",").count(PTPredicted)
        missingCnt = noneParentLine.split(",").count("NN")
        newHetCnt =(noneParentLine.split(",").count(HetPredicted1) + 
                 noneParentLine.split(",").count(HetPredicted2))
        
        #this ensures that with the missing data subtracted from only the 
        #PTxT43 (282) population that lines that are less than or equal to 
        #only the T43 count are removed so any lines retaining variation 
        # are kept.  In this case 16085 individuals are kept with the 
        # filtered_HapMap file (April 14 2019)
        if not (len(noneParentLine.split(","))-missingCnt <= majorCnt): 
            #This is where adjustments to how many possible errors are 
            #acceptable, this is also where the chance of cutting out rare 
            #alleles is greatest so the balance between type II and I error
            #is to be found here. Everything after this statement is fine 
            #tuning and will allow for smaller alterations.
            if (newHetCnt<=4 or minorCnt<=1):
                #Case 1: PTxT43 biallele is homozygous                  
                if PTxT43[0] == PTxT43[1]:                      #(10467 indivs)
                    #Case 1A: match to PT parent and not to T43                          
                    if PTxT43 in posPT and not PTxT43 in posT43:    #(0 indivs)
                        while True:
                            for i in range(len(posT43)):            #(0 indivs)
                                pair = list(posPT)[i]
                                #example AA[PT] x TT[T43] = AA [PTxT43]
                                if pair[0] == pair[1]:              #(1 indivs)
                                    parentDict[key]=ParentFileBuilderHelper(PTxT43Dict,PTDict,T43Dict)
                                    PTxT43FlagDict[key] = 1; break
                                elif pair == "NN":                  #(0 indivs)
                                    PTxT43FlagDict[key] = 0; break
                                elif pair[0] != pair[1]:            #(0 indivs)
                                    PTxT43FlagDict[key] = 0; break
                                else:
                                    print("ERROR1"); break
                            break
                    #Case 1B: match to PT parent and T43
                    elif PTxT43 in posPT and PTxT43 in posT43:   #(7123 indivs)
                        PTxT43FlagDict[key] = 0
                        #Since this case is impossible to tell the donor allele from the
                        #parent as it is a cross such as AA x AA = AA
                    #Case 1C: match to T43 parent and not to PT
                    elif PTxT43 in posT43 and not PTxT43 in posPT:#(3352 indivs)
                        while True:
                            for i in range(len(posPT)):                        
                                PTPair = list(posPT)[i]
                                #Case: PT is homozygous varient only
                                if PTPair[1] == PTPair[0] and len(posPT) == 1:
                                    parentDict[key]=ParentFileBuilderHelper(PTxT43Dict,PTDict,T43Dict)                        
                                    PTxT43FlagDict[key] = 1; break
                                #Case: PT is heterozygous varient only
                                elif PTPair[1] != PTPair[0] and len(posPT) == 1:
                                    parentDict[key]=ParentFileBuilderHelper(PTxT43Dict,PTDict,T43Dict)
                                    PTxT43FlagDict[key] = 1; break
                                #Case: PT is missing
                                elif PTPair == "NN":
                                    PTxT43FlagDict[key] = 0; break                    
                                else:
                                    PTxT43FlagDict[key] = 0; break 
                                    #This resolves all cases since if the cross is 
                                    #AA x __ = AA, the only options are AA x AA = AA
                                    #AA x AT = AA, or AA x TT = AA in which case it is 
                                    #viewed as a mismatch type error with a 0                        
                            break                             
                    #Case 1D: match to neither parent 
                    elif PTxT43 not in posPT and not PTxT43 in posT43:#(1 indivs)
                        while True:
                            #Check if one of the values is missing
                            for i in range(len(posPT)):
                                if(list(posPT)[i]) == "NN":
                                    PTxT43FlagDict[key] = 0; break
                                else: 
                                    PTxT43FlagDict[key] = 0
                                    break
                            break
                            for i in range(len(posT43)):
                                if(list(posT43)[i]) == "NN":
                                    PTxT43FlagDict[key] = 0; break
                                else:
                                    PTxT43FlagDict[key] = 0
                                    break
                            break
                        else:
                            PTxT43FlagDict[key] = 0               
                    else:
                        PTxT43FlagDict[key] = 0;
                #Case 2: PTxT43 biallele is heterozygous                       
                elif PTxT43[0] != PTxT43[1]:#(2 indivs)
                    while True:
                        #Check if one of the values is missing
                        for i in range(len(posPT)):
                            if(list(posPT)[i]) == "NN":
                                PTxT43FlagDict[key] = 0; break
                            else: 
                                PTxT43FlagDict[key] = 1; 
                                parentDict[key]=ParentFileBuilderHelper(PTxT43Dict,PTDict,T43Dict)
                                break
                        break
                        for i in range(len(posT43)):
                            if(list(posT43)[i]) == "NN":
                                PTxT43FlagDict[key] = 0; break
                            else: 
                                PTxT43FlagDict[key] = 1;
                                parentDict[key]=ParentFileBuilderHelper(PTxT43Dict,PTDict,T43Dict)
                                break
                        break
                    else:
                        PTxT43FlagDict[key] = 0
                else:
                    PTxT43FlagDict[key] = 0
            else:
                PTxT43FlagDict[key] = 0
        else:
            PTxT43FlagDict[key] = 0
        
#Make new hapmap file and output
cnt = 0
converCnt = 0
finalLineLST = []
for line in lineLST: 
    allInfoLST = line.strip().split(",")
    IDsLST = list(PTxT43outDict.keys())
    
    if(cnt <= 5):
        finalLineLST.append(line)
    else:   
        for ID in IDsLST:
            if allInfoLST[0] == ID and not PTxT43FlagDict[ID] == 0:
                finalLineLST.append(line)               
    cnt+=1
 
finCnt = 0
for line in finalLineLST:
    allInfoLST = line.strip().split(",")

    outLine = "" 
    if finCnt <= 4:
        allInfoLST.insert(5,"*")
        allInfoLST.insert(6,"*")
        allInfoLST.insert(7,"*")
        allInfoLST.insert(8,"*")
        allInfoLST.insert(9,"*")
        allInfoLST.insert(10,"*")
        allInfoLST.insert(11,"*")
        for item in allInfoLST:
            outLine = outLine + str(item) + ","
       
        rqtlOut.write(outLine[:len(outLine)-1]+"\n")
        PTorT43Out.write(outLine[:len(outLine)-1] + "\n")
            
    titleLine = ""
    if finCnt == 5:
        allInfoLST.insert(5,"HomoVHetVariation")
        allInfoLST.insert(6,"HomoVHetFlag(PTxT43)")
        allInfoLST.insert(7,"HomoVHetFlag(PT)")
        allInfoLST.insert(8,"HomoVHetFlag(T43)")
        allInfoLST.insert(9,"PTParentFreq")
        allInfoLST.insert(10,"ChromPosition")
        allInfoLST.insert(11,"ChromNum")
        for item in allInfoLST:
            titleLine = titleLine + str(item) + ","
            
        rqtlOut.write(titleLine+"\n")
        PTorT43Out.write(titleLine  + "\n")
     
    #This section extracts the alleles and the SNP IDs
    alleleInfoLine = ""
    if finCnt > 5:
        IDsLST = list(PTxT43outDict.keys())
        for ID in IDsLST:
           if ID == allInfoLST[0]:
               allInfoLST.insert(5,PTxT43VarDict[ID])
               allInfoLST.insert(6,PTxT43outDict[ID])
               allInfoLST.insert(7,PToutDict[ID])
               allInfoLST.insert(8,T43outDict[ID])               
               allInfoLST.insert(9,PTinPTxT43Dict[ID])
               allInfoLST.insert(10,posDict[ID])                
               allInfoLST.insert(11,chromDict[ID])
                     
        for item in allInfoLST:
            alleleInfoLine = alleleInfoLine + str(item) + ","
            
        if not allInfoLST[12] == "0": 
            alleleLine = ""
            if converCnt > 5:
                #These are duplicated so a single donor allele given, say A, turns 
                #into be an AA in this case for each of the parents
                PT =(parentDict[allInfoLST[0]].split(",")[0] + 
                     parentDict[allInfoLST[0]].split(",")[0])
                T43 =(parentDict[allInfoLST[0]].split(",")[2] + 
                      parentDict[allInfoLST[0]].split(",")[2])
                for i in range(len(allInfoLST)):
                    if i < 28:
                        alleleLine = alleleLine + str(allInfoLST[i]) + ","
                    else:
                        if allInfoLST[i] == PT and not allInfoLST[i] == T43:
                            alleleLine = alleleLine + "BB" + ","
                        elif allInfoLST[i] == T43 and not allInfoLST[i] == PT:
                            alleleLine = alleleLine + "AA" + ","                        
                        elif allInfoLST[i] == "NN" or allInfoLST[i] == "NN":
                            alleleLine = alleleLine + "NA" + ","   
                        elif allInfoLST[i][1] != allInfoLST[i][0]:
                            alleleLine = alleleLine + "AB" + "," 
                        else:
                            print("ERROR: " + allInfoLST[i])
                                
                
                PTorT43Out.write(alleleLine[:len(alleleLine)-1] + "\n")
            converCnt +=1
            
        rqtlOut.write(alleleInfoLine[:len(alleleInfoLine)-1]+"\n")
        
    finCnt +=1

#make parentID file  
for item in parentDict:
    outLine = "" + item +","+ parentDict[item] + "\n"
    parentOut.write(outLine)  
  
    
file.close()
BLASTFile.close()
parentOut.close()
rqtlOut.close() 
PTorT43Out.close()
