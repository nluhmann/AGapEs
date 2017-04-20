# in FPSAC, we can filter for families that have a higher mutltiplicity than 1, but we would still keep
# families that are not universal!



import sys

families = ["KIM","Nepal","Antiqua","Z1","Microtus","CO92","Pestoides","IP_31758","IP_32953","PB1","YPIII"]

families_filtered = sys.argv[1]
outfile = sys.argv[2]


file = open(families_filtered,"r")
first = True

out = open(outfile,"w")
for line in file:
	if line.startswith(">"):
		if not first:
			if len(counter) == 12:
				for elem in counter:
					out.write(elem)
		first = False
		counter = []
		counter.append(line)
	else:
		counter.append(line)
		

file.close()
