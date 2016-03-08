"""""
This script was designed for converting data into a usable SpaceMix format excluding the
R matrix needed to run it through R program.  Data should be made into a new file with all
missing data and both sample size and a random count of T,A,C,G alleles accounted for
Made by Kathryn Kananen 4/2/16
Updated by Kathryn Kananen 4/7/16
"""""
"""""
This script was designed for converting data into a usable SpaceMix format excluding the
R matrix needed to run it through R program.  Data should be made into a new file with all
missing data and both sample size and a random count of T,A,C,G alleles accounted for
Made by Kathryn Kananen 4/2/16
"""""
import sys, random

HufDoe_data = open(sys.argv[1]) #CombinedDD.txt
out = open(sys.argv[2], "w")


## makes column with name and genotype
for line in HufDoe_data:
	if not line.startswith("#"):
		samples = line.split("\t")
		samples_out = samples[0]
		geno_out = samples[12]
		data_out = (samples_out + "\t" + geno_out).replace("?", "N") #missing data replace
		out.write(data_out)
	
num_geno = len(geno_out.split(',')) #allele count
allele_choices = "AGTC" #creates a random allele list

## keeps alleles in pairs
formated_geno = geno_out.split(",") #formating locus to be together
for allele in formated_geno:
	geno = allele.strip()

""""" 
Issue: is printing 1 character not value as stated in geno previously, need to return 
one value for every pair (C_C) not (C). 
"""""
	
# scans every locus and puts in a List
LocusLST = [] 
selected_alleleLST = []
for locus in data_out:
	random_allele = random.choice(allele_choices)
	if locus == random_allele:
		selected_allele = random_allele
		selected_alleleLST.append(selected_allele) #adds random to List
		
"""""  
Issue: same as above as after it scans for one character I need it to place the pair it 
scanned into a list.  
"""""
HufDoe_data.close()
out.close()
