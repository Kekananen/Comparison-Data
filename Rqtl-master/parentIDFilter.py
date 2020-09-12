#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:40:10 2019
This script looks throught the ParentID.csv file and compiles a new HapMap 
file based on getting rid of monomorphic loci.

@author: Kathryn Kananen
"""
import sys

file = open(sys.argv[1], "r") #parentsID.csv
hapMapFile = open(sys.argv[2], "r") #SNPs_HapMap_DMz18-205v2.csv
out = open(sys.argv[3], "w") #SNPs_HapMap_DMz18-205v3.csv
PTorT43Out = open(sys.argv[4], "w") #SNPs_HAPMAP_DMz18-205v3_Parents.csv

GoodIDLST = [] #Holds all ID that are non monomorphic
parentDict = {}

for line in file:    
    if not line.startswith("ID"):
        allInfo = line.strip().split(",")
        #Where Types would be the parent allele
        PTType = allInfo[1]; T43Type = allInfo[3]; ID = allInfo[0]
        parentDict[ID] = PTType + "," + T43Type
        
        if not PTType == T43Type:
            GoodIDLST.append(ID)
            
cnt = 0
lineLST = []
for line in hapMapFile:
    
    if cnt > 5:
        allInfo = line.strip().split(",")
        ID = allInfo[0]
        
        if ID in GoodIDLST:
            out.write(line)
            lineLST.append(line)
    else:
      out.write(line) 
      lineLST.append(line)
      
    cnt +=1 
    
#This writes the PTorT43Out file which holds the converted information about
#PT in the PTxT43 lines.  For every pair the PT or T43 allele is represented 
#if present with a BB for PT a AA for T43, AB for a Hetrozygote and NN - for 
#missing data.  All other cases are removed.  This file also has been updated 
#to have the chromosome number and position data.  
cnt = 0
for line in lineLST: 
    allInfoLST = line.strip().split(",")
    if not allInfoLST[12] == "0":
    
        starredLines = "" 
        if cnt <= 4:
            for starredItem in allInfoLST:
                starredLines = starredLines + starredItem + ","        
            PTorT43Out.write(starredLines[:len(starredLines)-1]+"\n") 
     
        titleLine = ""
        if cnt == 5:
            for item in allInfoLST:
                titleLine = titleLine + item + "," 
    
            PTorT43Out.write(titleLine + "\n")
      
        alleleLine = ""
        if cnt > 5:
            #These are duplicated so a single donor allele given, say A, turns 
            #into be an AA in this case for each of the parents
            PT =(parentDict[allInfoLST[0]].split(",")[0] + 
                 parentDict[allInfoLST[0]].split(",")[0])
            T43 =(parentDict[allInfoLST[0]].split(",")[1] + 
                  parentDict[allInfoLST[0]].split(",")[1])
            for i in range(len(allInfoLST)):
                if i < 28:
                    alleleLine = alleleLine + allInfoLST[i] + ","
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
            
#            #removes odd case without PT in PTxT43 population
#            #Only look at those individuals in the PTxT43 population
#            noneParentLine = alleleLine.split(",")[:len(alleleLine.split(","))-11][:len(alleleLine.split(","))-10]           
#            del noneParentLine[len(noneParentLine)-2]
#            noneParentLine = ",".join(noneParentLine)
#
#            AAcnt = noneParentLine.split(",").count("AA")
#            BBcnt = noneParentLine.split(",").count("BB")
#            NAcnt = noneParentLine.split(",").count("NA")
#            ABcnt = noneParentLine.split(",").count("AB")
#            
#            #the raw number on the AAcnt side repesents the number of 
#            #that are being allowed for as a buffer of error
#            if not (len(noneParentLine.split(","))-30-NAcnt <= AAcnt):
#                if (ABcnt<=3 or BBcnt<=1):
#                    PTorT43Out.write(alleleLine[:len(alleleLine)-1]+"\n")
#                else:
#                    print(alleleLine[:len(alleleLine)-1])
        PTorT43Out.write(alleleLine[:len(alleleLine)-1]+"\n")
        cnt +=1 
        
file.close()
hapMapFile.close()
PTorT43Out.close()
out.close()
