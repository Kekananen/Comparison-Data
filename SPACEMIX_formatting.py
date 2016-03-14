"""""
This script was designed for converting data into a usable SpaceMix format excluding the
R matrix needed to run it through R program.  Data should be made into a new file with all
missing data and both sample size and a random count of T,A,C,G alleles accounted for
Made by Kathryn Kananen 3/2/16
updated by Kathryn Kananen 3/7/16
updated by Kathryn Kananen 3/13/16
"""""
import sys, random

HufDoe_data = open(sys.argv[1]) #CombinedDD.txt
out = open(sys.argv[2], "w")
count_out = open(sys.argv[3], "w")


## makes column with name and genotype
genoLST = []
for line in HufDoe_data:
	if not line.startswith("#"):
		samples = line.split("\t")
		samples_out = samples[0]
		geno_out = samples[12]
		data_out = (samples_out + "\t" + geno_out).replace("?", "N") #missing data replace
		genoLST.append(data_out.strip())

# geno out is last line so fix in for loop below so reference the list made above	

geno_format1 = geno_out.split(",") # keeps alleles in pairs
num_geno = len(geno_out.split(',')) #allele count

allele_choices = "AGTC" #creates a random allele list

# Creates random allele instance LST for all values in line 
selected_alleles = []
for i in range(num_geno):
 	random_allele = random.choice(allele_choices) # random allele instance
	selected_alleles.append(random_allele)

#Creates list of counts of alleles using random alleles
locus_countLST = [] # for if they are found with relative alleles
allele_countLST = [] # allele counts for if they are present 
for indiv in genoLST:
	new_locus = indiv.split("\t")
	geno = new_locus[1]
	smallLST = []
	small_countLST = []
	cnt = 0 #counter 
	for allele in geno.split(","):
		first_allele = allele[0]
		last_allele = allele[2]
		missing_allele = "N"
		if first_allele == selected_alleles[cnt]:
			alleleC = 1
		elif first_allele != selected_alleles[cnt]:
			alleleC = 0
		if last_allele == selected_alleles[cnt]:
			alleleCL = 1
		elif last_allele != selected_alleles[cnt]:
			alleleCL = 0
			
		#runs for alleles if they are present
		if first_allele != missing_allele:
			count = 1
		elif first_allele == missing_allele:
			count = 0
		if last_allele != missing_allele:
			count2 = 1
		elif last_allele == missing_allele:
			count2 = 0
		else:
			print("Fix Me :'[")
			
		locus_count = alleleC + alleleCL
		allele_count = count + count2
		smallLST.append(locus_count)
		small_countLST.append(allele_count)
		
		cnt += 1
	locus_countLST.append(smallLST)
	allele_countLST.append(small_countLST)
out.write(repr(allele_countLST))
count_out.write(repr(locus_countLST))
		
HufDoe_data.close()
out.close()
count_out.close()
