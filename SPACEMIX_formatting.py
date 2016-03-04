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

## creates a random allele generator		
allele_choices = "AGTC"
random_allele = random.choice(allele_choices)
formated_geno = geno_out.split(",") #allows for geno to be used in pairs Ex:(G_G)

for allele in formated_geno:
	if random_allele == allele[0]:
		count = +1
		if allele[2] == random_allele:
			count =+1
			print(count)
		else:
			print(count)
			pass
	else:
		count = 0
		print(count)
				
		

## sets generator to sample a single line
# for allele in geno_out:
# 	print(allele)
	
	





HufDoe_data.close()
out.close()
"""""
# 1. Remove all extra data from new file (ignore not remove)
# 		- Import title line

# 2. Create a new file 
# 		- Make column start with Sample in first row then 
# 		name of individuals (good)
# 		- strip genotypes of _ and / 
# 		- Add genotypes into file
		
# 3. Replace all missing data with N (good) 
		
4. Randomly sample one allele at each locus (T for one row, 
   G for another row, A for a different row)
   		- List [ATGC] randomly select for one of four
   			> extract data with first iteration
   			> take in names (locus) as a list
   			> run through list of list save random seq for next line
   
5. Turn counted alleles into integers if they are present in
   that row (if no T present in row then return 0)


6. Make new file 
		- same header as the file created before
//// 4,5,6 all in same for loop	keep in mind to convert missing data to0 


later script
////8. List must be in R matrix form so num[1:30, 1:10000] 
////9. If using locations make a matrix of long, lat (waiting on David)
////10. Run through like normal
"""""
# for line in HufDoe_data:
# 	if line.startswith('#name'):
# 		geno = line[11:]
# 		geno_split = geno.split(":")[-1].strip().split(",")
# 	
# # 		split_geno = "Sample\t%s\ge
# #  		split_geno = "%s\Sample\t%s\n" % (geno_split)
# # 		out.write(split_geno)
# # 	else:
# # 		pass