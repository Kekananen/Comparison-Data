"""""
Made to run FST for ZeaAllInfo populations of relavance
created by Kathryn Kananen 4/16/15
"""""
import sys


zeaAllInfo14col = open(sys.argv[1])
out = open(sys.argv[2], "w")
inNumbers1 = sys.argv[3]
inNumbers2 = sys.argv[4]


for line in zeaAllInfo14col:
	line2 = line.split("\t")
	if line2[0] == "pop": #title line
		out.write(line)
 	if line2[0] == inNumbers1:
		out.write(line)
	if line2[0] == inNumbers2:
		out.write(line)
	else:
		pass

out.close()
zeaAllInfo14col.close()
