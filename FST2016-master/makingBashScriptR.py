"""""
This script was created to make a bash script for the 
ZeaAllInfo.14col.allHybVeachPop.heirFstat file.
Created by Kathryn Kananen 6/4/15
"""""
import sys
from num2words import num2words

out = open(sys.argv[1], "w")

path = "~/Documents/Hufford_labs/David_data/chapter2/allPopFSTRun/part3FSTRun/"
pythonPath = path + "FSTpart2Run.R"
# dataPath = path + "ZeaAllInfo.14col.allHybVeachPop.heirFstat"

pastNum = []
# creates number range 106
for number in range(1,36):
	for nextNumber in range(1,36):
		if number != nextNumber:
			if nextNumber in pastNum:
			
				# formats first numbers
				spelledNum = num2words(number)
				repStep1 = spelledNum.replace("-", "_")
				finRep = repStep1.replace(" ", "_")

				# formats second numbers
				nextSpelledNum = num2words(nextNumber)
				nextRepStep1 = nextSpelledNum.replace("-", "_")
				nextFinRep = nextRepStep1.replace(" ", "_")		

				# outputs final line format
				newLine = finRep + "_" + nextFinRep
				finalLine2 = "Rscript --vanilla %s %sHybridFSTRun.txt %sHybridFSTRun.out\n" % (pythonPath, newLine, newLine)
				out.write(finalLine2)
				
	pastNum.append(number)


out.close()

# line format should go python, pathway, python script, inFile, outFile, num1, num2
# python ~/Users/User/Documents/Hufford_labs/David_data/chapter2/allPopFSTRun/HybVEachPop_heirFStat python BreakingFStatHybV.py ZeaAllInfo.14col.allHybVeachPop.heirFstat one_twoHybridFSTRun.txt 1 2 

# 
# path = "~/Documents/Hufford_labs/David_data/chapter2/allPopFSTRun/HybVEachPop_heirFStat/"
# pythonPath = path + "BreakingFStatHybV.py"
# dataPath = path + "ZeaAllInfo.14col.allHybVeachPop.heirFstat"
# 
# pastNum = []
# # creates number range 106
# for number in range(1,107):
# 	for nextNumber in range(1,107):
# 		if number != nextNumber:
# 			if nextNumber in pastNum:
# 			
# 				# formats first numbers
# 				spelledNum = num2words(number)
# 				repStep1 = spelledNum.replace("-", "_")
# 				finRep = repStep1.replace(" ", "_")
# 
# 				# formats second numbers
# 				nextSpelledNum = num2words(nextNumber)
# 				nextRepStep1 = nextSpelledNum.replace("-", "_")
# 				nextFinRep = nextRepStep1.replace(" ", "_")		
# 
# 				# outputs final line format
# 				newLine = finRep + "_" + nextFinRep
# 				finalLine2 = "python %s %s %sHybridFSTRun.txt %s %s\n" % (pythonPath, dataPath, newLine, str(number), str(nextNumber))
# 				out.write(finalLine2)
# 				
# 	pastNum.append(number)
# 
# 
# out.close()
# 
# # line format should go python, pathway, python script, inFile, outFile, num1, num2
# # python ~/Users/User/Documents/Hufford_labs/David_data/chapter2/allPopFSTRun/HybVEachPop_heirFStat python BreakingFStatHybV.py ZeaAllInfo.14col.allHybVeachPop.heirFstat one_twoHybridFSTRun.txt 1 2 