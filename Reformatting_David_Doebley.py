"""""
This script is designed to transfer two files in this case the ZeallInfo.txt file 
and TPG_PM file into one txt file for a spacemix analysis
created by Kathryn Kananen 2/18/2016
"""""

import sys
TPG_PM = open(sys.argv[1])
ZeaInfo = open(sys.argv[2])
out = open(sys.argv[3], 'w')


def RemoveSpaces(text):
    text = text.strip()
    now = text.split("\t")
    new = []
    for word in now:
        temp = word.strip().replace('/','_')
        new.append(temp)
    
    return(new)

#Output header info
out.write("#python %s\n" % (" ".join(sys.argv)))

ZeaInfo.readline()
titles = ZeaInfo.readline()
out.write(titles)

# makes Doebley data genotypes and name be transfered into new file
for line in TPG_PM:
	if not line.startswith('OTU'):
		Doebley_data = RemoveSpaces(line)
		name = Doebley_data[0]
		geno = ",".join(Doebley_data[2:])
		newline = "%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%s\n" % (name, geno)
		out.write(newline)

# makes Hufnagel data genotypes and name be transfered into new file 
for line in ZeaInfo:
	if not line.startswith('#'):
		out.write(line)

TPG_PM.close()
ZeaInfo.close()
out.close()

