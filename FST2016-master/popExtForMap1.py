# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 08:12:05 2016

This script was created for seperating out parv, mex, and hybrid populations
with the goal of mapping them later in R.

@author: Kathryn Kananen
"""
import sys

zeaAllInfo = open(sys.argv[1])
outMex = open(sys.argv[2], "w")
outParv = open(sys.argv[3], "w")
outHyb = open(sys.argv[4], "w")

# Adds headers to all files
header = "#name\ttax\ttype\tcoordinates\n"        
outParv.write(header)
outMex.write(header)
outHyb.write(header)
        




# scans through the ZeaAllInfo file
for line in zeaAllInfo:
    # removes the header
    if not line.startswith("#"):
        infoLine = line.split("\t")
        
        # Where the variables of interest are made
        name = infoLine[1]
        tax = infoLine[2]
        type = infoLine[3]
        coords = infoLine[5]+"\t"+infoLine[6]
        
        # For parv pops
        if type == "parvHC" or (type == "other" and tax == "Zea_mays_parviglumis"):
            desired = "%s\t%s\t%s\t%s\n" % (name, tax, type, coords)
            
            #writes line to file
            outParv.write(desired)
                # For parv pops
        if type == "mexHC" or (type == "other" and tax == "Zea_mays_parviglumis"):
            desired = "%s\t%s\t%s\t%s\n" % (name, tax, type, coords)
            
            #writes line to file
            outMex.write(desired)
        if type == "hybrid":
            desired = "%s\t%s\t%s\t%s\n" % (name, tax, type, coords)
            
            #writes line to file
            outHyb.write(desired)
            
outParv.close()
outMex.close()
outHyb.close()
zeaAllInfo.close()