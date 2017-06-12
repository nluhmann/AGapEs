# in FPSAC, we can filter for families that have a higher mutltiplicity than 1, but we would still keep
# families that are not universal!



import sys

#families = ["KIM","Nepal","Antiqua","Z1","Microtus","CO92","Pestoides","IP_31758","IP_32953","PB1","YPIII"]



families_filtered = sys.argv[1]
species_file = sys.argv[2]
outfile = sys.argv[3]

families = []


species = open(species_file,"r")
for line in species:
	if not line.startswith("#"):
		families.append(line.rstrip("\n"))

species.close()


file = open(families_filtered,"r")
first = True

out = open(outfile,"w")
for line in file:
	if line.startswith(">"):
		if not first:
			if len(check_hash.keys()) == len(families) and len(counter) == len(families) +1:
				for elem in counter:
					out.write(elem)
				out.write("\n")
		first = False
		check_hash = {}
		#check_hash[line] = 1
		counter = []
		counter.append(line)
	elif not line == "\n":
		check_hash[line.split(".")[0]] = 1
		counter.append(line)

		
if len(check_hash.keys()) == len(families) and len(counter) == len(families) +1:
	for elem in counter:
		out.write(elem)
file.close()
