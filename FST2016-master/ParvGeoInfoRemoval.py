"""""
This script is designed to remove information from the Genotypes file
The information is stored into a dictionary and then scans the bigInfoFile 
and seperates the the individauls by acession, check the confidence 
intervals, and then create a new list of those individuals
Created by Kathryn Kananen on 9/8/16
Finished by David E. Hufnagel on Sep 20, 2016
Updated on 11/10/16 by Kathryn Kananen
The script now only works for pure pops with modification to part 3
"""""
import sys

filteredRaces = open(sys.argv[1]) #genotypes_filtered.races
ZeaAllInfo = open(sys.argv[2])    #ZeaAllInfo.14col.noMML
out = open(sys.argv[3], 'w')      #pureParvPops.info or allParvPops.out


def SaveIntoDict(key, val, dictx):
    if not key in dictx:
        dictx[key] = [val,]
    else:
        dictx[key].append(val) 
        
        
###### Part 1 ###### 
# Holds the name of individuals in the population group such as North_Balsas.
# Should be used to get the race name of the individual.
filteredGenoDict = {}

# Uses a counter to get single information points out of filtered file by checking
# if the counter is at a certain point and grabbing information at the value 
# (index) of the counter.  Adds values to Dictionary.
count = 0
for line in filteredRaces:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        outRace = lineLst[2] # Individuals race
        outIndiv = lineLst[0] # Individuals name  
        filteredGenoDict[outIndiv] = outRace

###### Part 2 ######
# Holds information of each accession with high confidence individuals to look
# like >>> accession:[[name,race,conf],[name,race,conf],[name,race,conf]] where
# each list is an individual value.
finalDict = {}
for info in ZeaAllInfo:
    # Ignore the title line
    if not info.startswith("#"):
        allInfoLine = info.strip().split("\t")

        #Gets all necessary information from ZeaAllInfo file
        name = allInfoLine[0]
        tax = allInfoLine[1]
        accession = allInfoLine[2]
        hybrid = allInfoLine[3]
        race = filteredGenoDict[name]
        if not name in filteredGenoDict:
            race = filteredGenoDict[name]
        else:
            race = "NA"

#  		Excludes all non parv of mex individuals		
        if(tax == "Zea_mays_parviglumis"):
            allInfo = (name, race, hybrid)
            SaveIntoDict(accession, allInfo, finalDict)
print(finalDict)

###### Part 3 ######
# Checks to see if a group has a value that is a hybrid or not available and outputs 
# that group to a new file. Uses lists to cross reference each other.
out.write("#accession\tpopSize\tstatus\n")
for acc,inds in finalDict.items():
    isGood = True #starts true and is turned false if a nonpure individual is found in the pop
    for ind in inds:
        name = ind[0]; race = ind[1]; conf = ind[2]
        if not (conf == "mexHC" or conf == "parvHC"):
            isGood = False
    
    if isGood == True: #Remove this if statement for all parv population
        popSize = len(inds)
        newLine = "%s\t%s\t%s\n" % (acc,popSize,conf)
        out.write(newLine)



filteredRaces.close()
ZeaAllInfo.close()	
out.close()	