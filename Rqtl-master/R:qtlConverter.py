#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 10:36:50 2019
Converts the SNPs_HAPMAP_DMz18-205v2.csv to the read.cross file format needed
for rqtl analysis.  

The script takes in four files - two input and two output:
    1. SNPs_HAPMAP_DMz18-205v2.csv
    2. Metepec_pigment_intensity.csv
    3. The rqtl output file
    4. The file used to keep track of the order that holds the information of 
    the rows from the PTxT43 population.  These need to be repopulated from 
    the Metepec_pigment_intensity.csv as the order shifts.

The ouput file should look like this using my naming conventions below:
    
    Label   ID     ID     ID
            chrom  chrom  chrom
            pos    pos    pos
    row     AA     AB     AB
    row     AA     AA     AA
    row     NA     AB     NA
    
The row in this instance can be tracked through the orderFile, or the 4th file
given to the script, so the trait can still be tracked using the alleleIDs as
the markers in the rqtl output

If this is confusing, then run the script and look at the files to try to see 
how the variables connect to one another.

@author: Kathryn Kananen
"""
import sys

HAPMAPfile = open(sys.argv[1]) #SNPs_HAPMAP_DMz18-205v3.csv
phenoFile = open(sys.argv[2]) #Metepec_pigment_intensity.csv
out = open(sys.argv[3], "w") #Whatever you want your output file to be called
orderFile = open(sys.argv[4], "w") #File that holds order of families

"""""
Converts the alleles given as a tuple into a list for other methods that rely
on list based logic.
    Given: The alleles as a Tuple
    Output: The converted Alleles as a list
"""""
def AlleConvert(alleles):
    outLST = []
    for alleleLine in alleles:
        outLine = ""
        for allele in alleleLine:
            outLine = outLine + allele
        outLST.append(outLine + "\n")
        outLine = ""
        
    return outLST

#The alleleIDs were chosen as there was no other distinguishing marker in the 
#dataset that marked all of the SNPs that was 100% distintly unique.
IDs = "" #The alleleIDs 
rows = "" #The phenotype information in this case intensity
chroms = "," #Chromosome number
positions = "," #Location on the Chromosome

AlleleLST = [] #Holds the unconverted alleles from the HAPMAPfile

IDDictInfo = {} #Holds the position and chromosome for each SNP

cnt = 0
for line in HAPMAPfile:
    AlleleRows = ""
    #The info that is not title related
    if not line.startswith("*") and cnt > 5:
        allInfo = line.strip().split(",")   
        IDDictInfo[str(allInfo[0])] = allInfo[12] + "," + allInfo[11]
        
        ID = allInfo[0]
        IDs = IDs + ID + ","                 
        
        baseLST = [];
        for i in range(len(allInfo)):
            if i > 28:
                baseLST.append(allInfo[i] + ",")
        AlleleLST.append(baseLST)
        baseLST = []

    #Grabs the PTxT43, PT, and T43 row markers 
    if cnt == 4:
        rowsLST = line.strip().split(",")
        for i in range(len(rowsLST)):
            if rowsLST[i] != "*":    
                rows = rows + "," + rowsLST[i]
    cnt += 1

#Reverses the order of the columns and rows so they fit with the new labeling
RevAlleleLST = zip(*AlleleLST)
outAllelesLST = AlleConvert(RevAlleleLST) 

for ID in IDDictInfo:
    chroms = chroms + IDDictInfo[ID].split(",")[0] + ","
    positions = positions + IDDictInfo[ID].split(",")[1] + ","

#combine all of the information together to the output file
IDLine = "Average_Intensity,"+IDs[:len(IDs)-1]+ "\n"
chromLine = chroms[:len(chroms)-1] + "\n" 
posLine = positions[:len(positions)-1] + "\n"

out.write(IDLine)
out.write(chromLine)
out.write(posLine)

phenoRowsDict = {}
for line in phenoFile:
    if not line.startswith("PTxT43"):
        allInfo = line.strip().split(",")
        #Values that are not existant, skip them and assign them 
        #to be intensity of 0
        if not allInfo[1] == "#VALUE!":
            phenoRowsDict[allInfo[0]] = allInfo[1]
        else:
            phenoRowsDict[allInfo[0]] = "0"
        
outLine = ""
phenoLST = [] #List to be output in the orderFile
rowLST = rows[1:].split(",")
for row in rowLST:
    rowNum = row.split(" ")[-1:][0]
    #Only include rows that are PTxT43
    if rowNum.isdigit():
        if rowNum in list(phenoRowsDict.keys()):
            outLine = phenoRowsDict[rowNum] + ","
            #add the row number and the intensity to phenoLST
            phenoLST.append(rowNum + "_" + phenoRowsDict[rowNum])
            for i in range(len(rowLST)):
                compRowNum = rowLST[i].split(" ")[-1:][0]
                #If the comparison row is equal to the row number, then it is 
                #PTxT43 and then it can be added to the outLine
                if compRowNum == rowNum:
                    outLine = outLine + outAllelesLST[i]
            out.write(outLine[:len(outLine)-2]+"\n")
        else:
            print("ERROR")

#Write the output file for the orderFile
orderFile.write("ID_Trait\n")
for pheno in phenoLST:
    orderFile.write(pheno + "\n")
    
out.close()
phenoFile.close()
HAPMAPfile.close()
orderFile.close()
