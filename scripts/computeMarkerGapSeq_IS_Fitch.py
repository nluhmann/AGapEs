import sys
import os


gap = sys.argv[1]
gap_coords = sys.argv[2]
alignment = sys.argv[3]
families = sys.argv[4]
assemblypath = sys.argv[5]
#ancestral_gaps = sys.argv[4:]






# 1) read gap coordinates to get connection between adjacencies and gapID (we do not really need the coordinates here!)

# 2) read sequences for families (family IDs are not doubled!)

# 3) read adjacencies: for each adjacency, get marker sequences and gap sequence -> this is the assembly
#		then compute mapping of reads to this assembly sequence, keep only mapped ones



#############################################################################################  	
### return sequence complement
baseComplement = {"A":"T","T":"A","C":"G","G":"C","N":"N","M":"M"}
def complement(sequence):
    result = []
    for n in range(len(sequence)):
        result.append(baseComplement[sequence[-(n+1)]])
    return ''.join(result)




# 1)
coords = open(gap_coords,"r")
coordsHash = {}
adjHash = {}
for line in coords:
	if line.startswith(">"):
		array = line.rstrip("\n").split(" ")
		adj = array[2]
		left = adj.split("-")[0]
		right = adj.split("-")[1]
		coordsHash[(left,right)] = array[0][1:]
		adjHash[array[0][1:]] = (left,right)
coords.close()


# 2)
fams = open(families,"r")
famsHash = {}
for line in fams:
	if line.startswith(">"):
		famID = line.rstrip("\n")[1:]
	else:
		seq = line.rstrip("\n")
		famsHash[int(famID)] = seq
fams.close()


# 3) read adjacencies and save assembly in a fasta file 
#		second file: save adjacency and gap information

# assemblypath = "./assemblies"
# if not os.path.exists(assemblypath):
#     os.makedirs(assemblypath)



#for gap in ancestral_gaps:
base = gap.split("/")[-1]
gap_number = base.split(".")[0]
gapID = gap_number.split("-")[0]

if "-" in gap_number:
	gap_version = gap_number.split("-")[1]
else:
	gap_version = ""

#	adj = adjHash[gapID]
#	left = adj[0]
#	right = adj[1]

left = adjHash[gapID][0]
right = adjHash[gapID][1]

if not os.path.isfile(alignment):
	#there is no gap, marker sequences directly adjacent
	gapSeq = ""	
else:
	with open(alignment, 'r') as content_file:
		gapSeq = content_file.read()
	gapSeq = gapSeq.rstrip("\n")


#get left marker sequence
if int(left) % 2 == 0:
	#just take sequence
	leftSeq = famsHash[int(left)/2]
else:
	#take sequence complement
	interimSeq = famsHash[(int(left)+1)/2]
	leftSeq = complement(interimSeq)
	#take also complement of gap sequence
	gapSeq = complement(gapSeq)
	
#get right marker sequence
if int(right) % 2 == 0:
	#take sequence complement
	interimSeq = famsHash[(int(right))/2]
	rightSeq = complement(interimSeq)
else:
	#just take sequence
	rightSeq = famsHash[(int(right)+1)/2]
	
assembly = leftSeq+gapSeq+rightSeq
gapStart = len(leftSeq) + 1
gapEnd = len(leftSeq)+len(gapSeq)

if not gap_version == "":	
	file = open(assemblypath+"/"+gapID+"-"+gap_version+".fasta_assembly","w")
else:
	file = open(assemblypath+"/"+gapID+".fasta_assembly","w")
file.write(">"+gapID+" "+left+" "+right+" "+str(len(assembly))+" "+str(gapStart)+"-"+str(gapEnd)+"\n"+assembly)
file.close()

		












