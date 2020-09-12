"""""
Designed to transfer locus data to matix for SpaceMix analysis of Hufnagel data.  Will 
create three seperate files that contain the location of each population selected with 
corresponding matrices of a randomly selected allele and if the matrix was genotyped at 
the loci.  Populations that are from the same coordinates are used as one population type,
meaning they are sampled together if they are from the same geographical location.

Made by Kathryn Kananen 4/2/16
Updated by Kathryn Kananen 4/6/16
Cleaned up by David E. Hufnagel on 4/9/16
Updated by Kathryn Kananen 4/9/16
Added a locations matrix and grouped populations together
"""""
import sys, random


#if present, linearly sum lists 
def SaveIntoDictCounts(key, val, dictx):
    if not key in dictx:
        dictx[key] = val
    else:
        dictx[key] = [x+y for x,y in zip(dictx[key], val)]

#if present, check that new and old vals are the same
def SaveIntoDictLocs(key, val, dictx):
    if not key in dictx:
        dictx[key] = val
    else:
    	####### FIX ########
        if dictx[key] != val:
            print("ERROR WAS HERE!")
            sys.exit()


ZeaAllInfo = open(sys.argv[1]) #ZeaAllInfo.txt
AllelesPresent = open(sys.argv[2], "w") #AllelesPresent.txt
AllelesRandom = open(sys.argv[3], "w") #AlleleRandom.txt
LocationsOut = open(sys.argv[4], "w") #locationData.txt

# Make a list of randomly chosen letters the same length as the number of markers
genoLST = []
for line in ZeaAllInfo:
    if not line.startswith("#"):
        samples = line.split("\t")
        samples_out = samples[0]
        geno_out = samples[12].strip()
        tax = samples[1]
        country = samples[7]
        longitude = samples[5]
        latitude = samples[6]
        acc = samples[2]
        data_out= "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (samples_out, geno_out, tax, country, longitude, latitude, acc)
        genoLST.append(data_out.strip())

geno_format1 = geno_out.split(",") # keeps alleles in pairs
num_geno = len(geno_out.split(',')) #allele count (968)
allele_choices = "AGTC" #creates a random allele list

selected_alleles = []
for i in range(num_geno):
    random_allele = random.choice(allele_choices) # random allele instance
    selected_alleles.append(random_allele)

# Creates list of lists for alleles present and alleles chosen 
locus_countLST = [] # for if they are found with relative alleles
allele_countLST = [] # allele counts for if they are present 
locationsLST = []# locations of individuals or populations
accLST = []
for indiv in genoLST:
	new_locus = indiv.split("\t")
	name = new_locus[0]
	geno = new_locus[1]
	tax = new_locus[2]
	country = new_locus[3]
	lat = new_locus[4]
	long = new_locus[5]
	acc = new_locus[6]
	smallLST = []
	small_countLST = []
	#used_taxLST = [] #list of only individuals to be used
	cnt = 0

	if tax == "Zea_mays_parviglumis" or tax == "Zea_mays_mexicana" or (tax == "Zea_mays_mays" and country == "Mexico"):
		coord = (lat + "\t" + long)
		locationsLST.append(coord)
		
		# creates acc list 
		if not acc in accLST:
			accLST.append(acc) 
		
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
		
		smallLST.insert(0, acc)
		small_countLST.insert(0, acc)
		locus_countLST.append(smallLST)
		allele_countLST.append(small_countLST)

# creates dic 
alleleDict = {}
for line in allele_countLST:
	SaveIntoDictCounts(acc, line[1:], alleleDict)
	
locusDict = {}
for line in locus_countLST:
	SaveIntoDictCounts(acc, line[1:], locusDict)

coordDict = {}
for line in locationsLST:
	SaveIntoDictLocs(acc, line[1:], coordDict)
	
print(coordDict)
#output data in reasonable format
for i in range(len(locus_countLST)): #the lists have the same dimensions so they are parsed simultaneously
    outLineLoc = "\t".join(map(str, allele_countLST[i])) + "\n" #the map function converts numbers to strings in this case
    allLineLoc = "\t".join(map(str, locus_countLST[i])) + "\n"
    AllelesPresent.write(outLineLoc)
    AllelesRandom.write(allLineLoc)

for coords in locationsLST:
	new_line = coords + "\n"
	LocationsOut.write(new_line)

ZeaAllInfo.close()
AllelesPresent.close() 
AllelesRandom.close() 
LocationsOut.close()
