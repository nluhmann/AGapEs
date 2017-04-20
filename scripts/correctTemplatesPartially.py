import sys
import os


#INPUT: fasta_assembly
assFile = sys.argv[1]

#INPUT: partiallyGaps_coordinates / notFilledGaps_coordinates
gapcoords = sys.argv[2]

#INPUT: corresponding Fitch template fasta
fitch = sys.argv[3]

#OUTPUT: file
output = sys.argv[4]

#OUTPUT: index file     !!!!!!!!!! is not overwritten, but appended....delete before new computation!!!!!!!!!!!!
indexfile = sys.argv[5]

#INPUT: all partially filled gaps (.out)
gaps = sys.argv[6:]




#############################################################################################  	
### return sequence complement
baseComplement = {"A":"T","T":"A","C":"G","G":"C","N":"N","D":""}
def complement(sequence):
    result = []
    for n in range(len(sequence)):
        result.append(baseComplement[sequence[-(n+1)]])
    return ''.join(result)





#############################################################################################  
#compute edit Distance between two strings

def editDistance(first,second):
	if len(first) < len(second):
		return editDistance(second, first)
    
    # now we can assume that len(second) <= len(first)
	if len(second) == 0:
		return len(first)
	previous_row = range(len(second) + 1)

	for i, characterOne in enumerate(first):
		current_row = [i + 1]
		for j, characterTwo in enumerate(second):
			insertions = previous_row[j + 1] + 1 
			deletions = current_row[j] + 1       
			substitutions = previous_row[j] + (characterOne != characterTwo)
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	return previous_row[-1]



# 1) read fasta_assembly
file = open(assFile, "r")
assemblySeq = ""
for line in file:
	if line.startswith(">"):
		savedLine = line
		array = line.split(" ")
		gapID = array[0][1:]
		coords = array[4]
		gapStart = int(coords.split("-")[0])
		gapEnd = int(coords.split("-")[1])
	else:
		assemblySeq = assemblySeq + line.rstrip("\n")
file.close()

length = len(assemblySeq[gapStart-1:gapEnd])



#exit if there is no partially filled sequence
if "*" in gaps[0]:
	print "EXIT "+gaps[0]
	index = open(indexfile,"a")
	index.write("UNCOVERED"+"\t"+gapID+"\t"+"0"+"\t"+str(length)+"\n")
	index.close()
	sys.exit()




def sortkey(tup):
	if "/" in tup:
		basename = tup.split("/")[-1]
	else:
		basename = tup
	return int(basename.split("_")[-1].split("-")[0])

gapsSorted = sorted(gaps, key=sortkey)




# 1.1) read fitch template

if os.path.isfile(fitch):
	
	with open (fitch, "r") as file:
		template_fitch = file.read().replace('\n', '')
	file.close()
else:
	print "No template sequence available"
	template_fitch = ""

 

# 1.2) Read coordinates from file, save in list and cut assemblySeq accordingly


splitcoords = []
file = open(gapcoords,"r")
for line in file:
	if line.startswith(gapID+"\t"):	
		splitcoords.append(line.rstrip("\n"))
file.close()

lastStop = 0
assemblysplit = []
for splitcoord in splitcoords:
	array = splitcoord.split("\t")
	start = int(array[1].split("-")[0])
	stop = int(array[1].split("-")[1])
	if start == stop and stop == gapEnd:
		print "SONDERFALL"
		start = start - 1
	if array[2] == "LEFT":
			start = start
			stop = stop -1
	else:
		stop = stop -1

	first = assemblySeq[lastStop:start]
	second = assemblySeq[start:stop]
	assemblysplit.append(first)
	assemblysplit.append(second)
	
	lastStop = stop
	# if second == "":
# 		lastStop = stop + 1
# 	else:
# 		lastStop = stop

	
# if second == "":
# 	lastStop = lastStop - 1

	
	
assemblysplit.append(assemblySeq[lastStop:])


# file = open(gapcoords,"r")
# for line in file:
# 	if line.startswith(gapID+"\t"):	
# 		break
# file.close()
# 
# splitcoords = line.rstrip("\n").split("\t")
# assemblysplit = []
# lastStop = 0
# if gaps[0].split("_")[-1] == "0.out":
# 	firstbool = True
# else:
# 	firstbool = False
# 	
# for c in splitcoords[1:]:
# 	start = int(c.split("-")[0])
# 	stop = int(c.split("-")[1])
# 	if start == stop and stop == gapEnd:
# 		start = start - 1
# 	if firstbool:
# 		start = start - 2
# 		stop = stop - 1
# 		firstbool = False
# 
# 		
# 	else:		
# 		stop = stop -1
# 
# 	first = assemblySeq[lastStop:start]
# 	second = assemblySeq[start:stop]
# 	assemblysplit.append(first)
# 	assemblysplit.append(second)
# 	if second == "":
# 		lastStop = stop + 1
# 	else:
# 		lastStop = stop
# 
# 	
# if second == "":
# 	lastStop = lastStop - 1
# 
# 	
# 	
# assemblysplit.append(assemblySeq[lastStop:])

 
 
# 2) for each partially filled sequence, split assemblySeq accordingly, check if distance agrees and then substitute.
# 
newAssembly = assemblySeq
totaldistance = 0
coveredLength = 0
counter = 1
for gap in gapsSorted:
	part = open(gap,"r")
	#counter = (int(gap.split(".")[0].split("_")[-1])*2)+1
	for line in part:
		if line.startswith(">"):
			array = line.rstrip("\n").split(" ")
			distance = int(array[5])
			totaldistance = totaldistance + distance
			start = int(array[6])
			stop = int(array[7])
		else:
			readsequence = line.rstrip("\n")
	part.close()	

	# print editDistance(assemblysplit[counter],readsequence)
        if not readsequence == "":
            assemblysplit[counter] = readsequence
        else:
            print "NO GAPFILLING "+gap
 	coveredLength = coveredLength + len(readsequence)
 	counter = counter + 2
# 	print "READ"
# 	print readsequence
# 	print "OLD"
# 	print assemblysplit[counter]

# 3) rejoin assembly sequence to test distance (should equal to totaldistance now)
assSeq = "".join(x for x in assemblysplit if x is not '')
#dist = editDistance(assSeq, assemblySeq)
dist = totaldistance
if not dist == totaldistance:
	print "ERROR "+gapID+"\n"
	print assFile
	print assemblysplit
	print dist
	print totaldistance
	print "NEW"
	print assSeq
	print "OLD"
	print assemblySeq
	index = open(indexfile,"a")
	index.write("UNCOVERED"+"\t"+gapID+"\t"+"0"+"\t"+str(length)+"\n")
	index.close()
elif coveredLength == 0:
	index = open(indexfile,"a")
	index.write("UNCOVERED"+"\t"+gapID+"\t"+"0"+"\t"+str(length)+"\n")
	index.close()
else:
	print "CHECK "+gapID+"\n"
	#print assemblysplit
	
	gapSequence = assSeq[gapStart-1:gapEnd]
	distTemplate = editDistance(gapSequence, template_fitch)
	revGapSequence = complement(gapSequence)
	revDistTemplate = editDistance(revGapSequence, template_fitch)
	distTemplate = 0
	revDistTemplate = 1
	# 4) write fasta_ancestral file
	out = open(output,"w")
	#out.write(savedLine)
	if distTemplate < revDistTemplate:		
		out.write(gapSequence)
	else:
		out.write(revGapSequence)
	out.close()		
	
	index = open(indexfile,"a")
	index.write("PARTIALLY"+"\t"+gapID+"\t"+str(coveredLength)+"\t"+str(len(gapSequence))+"\n")
	index.close()
		
		
		
		
		
