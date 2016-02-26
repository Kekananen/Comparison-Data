"""
This script is designed to compare differences between 
the Goebley and Hufnagel data by using a dictionary and 
printing to output.
created by Kathryn Kananen 1/29/2016

Creates a dictionary and runs through them in output 
"""
import sys

# Activate: python Comparison_Data.py _______.txt ________.txt Random_name.txt
Doebley2 = open(sys.argv[1]) 	#TPG050110.txt
Doebley1 = open(sys.argv[2]) 	#TPG_PM.txt
out = open(sys.argv[3], 'w') 

## Remove whitespace 
def RemoveSpaces(line):
    line = line.strip()
    line = line.split()
    line_list = []
    for word in line:
    	New_word = word.strip()
    	line_list.append(New_word)
    return(line_list)

#removes whitespace in TPG050110.txt and TPG_PM.txt
Doebley2.readline()
for line in Doebley2:
	Doebley2_file = RemoveSpaces(line)
	sys.exit()
Doebley1.readline()
for line in Doebley1:
	Doebley1_file = RemoveSpaces(line)
	sys.exit()

# Creates Dictionary for TPG050110.txt
DoebleyDict={}
for line in Doebley2:
	Doebley2_line = line[0]
	Doebley2_line2 = line[1]
	DoebleyDict[Doebley2_line]=Doebley2_line2
	print(DoebleyDict)


Doebley2.close()
Doebley1.close()
out.close()
