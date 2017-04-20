#!/usr/bin/env python

import sys
import shortestPath_indels
import re
from operator import itemgetter
import time
import os


#INPUT: sam file with mappings to potential assembly
mappingFile = sys.argv[1]


#INPUT: corresponding Fitch template fasta
#fitch = sys.argv[2]

#INPUT: potential assembly fasta file (we need the start and end of each gap)
assFile = sys.argv[2]

#output file
output = sys.argv[3]

t0 = time.time()

#############################################################################################  	
### return sequence complement
baseComplement = {"A":"T","T":"A","C":"G","G":"C","N":"N","D":""}
def complement(sequence):
    result = []
    for n in range(len(sequence)):
        result.append(baseComplement[sequence[-(n+1)]])
    return ''.join(result)




#############################################################################################  
### get template sequence for gap from fitch reconstruction

# if os.path.isfile(fitch):
# 	
# 	with open (fitch, "r") as file:
# 		template = file.read().replace('\n', '')
# 	file.close()
# else:
# 	print "No template sequence available"
# 	template = ""

#############################################################################################  
### get start and stop for gap from assembly file

file = open(assFile, "r")
assemblySeq = ""
for line in file:
	if line.startswith(">"):
		savedLine = line
		array = line.split(" ")
		coords = array[4]
		gapStart = int(coords.split("-")[0])
		gapEnd = int(coords.split("-")[1])
	else:
		assemblySeq = assemblySeq + line.rstrip("\n")
file.close()

template = assemblySeq[gapStart -1:gapEnd]

#############################################################################################  
### test coverage for gap

def testCoverage(mappings, gapStart, gapEnd):
	length = gapEnd - gapStart	
	covArray=[0 for i in range(length)]
	breakpointArray = [[] for i in range(length)]
	for read in mappings:
		id = read[0]
		start = int(read[1])
		cigar = read[4]
		matchM = re.findall(r'(\d+)M', cigar)
		if len(matchM) > 1:
			l = int(matchM[0]) + int(matchM[1])
		else:
			l = int(matchM[0])		
		stop = start+l
		if (start <= gapStart) and (stop >= gapStart) and (stop <= gapEnd):
			overlap = stop - gapStart 
			for i in range (0, overlap):
				covArray[i] = covArray[i] + 1
				breakpointArray[i].append(id)
		elif start < gapEnd and start > gapStart and stop > gapEnd:
			overlap = gapEnd - start
			for i in range(1,overlap):
				covArray[-i] = covArray[-i] + 1
				breakpointArray[-i].append(id)
		elif start > gapStart and stop < gapEnd:
			diffStart = start - gapStart
			diffEnd = stop - gapStart
			for i in range (diffStart, diffEnd):
				covArray[i] = covArray[i] + 1
				breakpointArray[i].append(id)
		elif start < gapStart and stop > gapEnd:
			for i in range(len(covArray)):
				covArray[i] = covArray[i] + 1
				breakpointArray[i].append(id)

	#print "------------ TEST"
	#print "Coverage per gap base:"
	#print covArray
	
	if 0 in covArray:
		print "Yes, gap is not covered."
		pos = [i for i, e in enumerate(covArray) if e == 0]
		print "Gap: "+str(gapStart)+"-"+str(gapEnd)
		print "Not covered positions in gap: "
		print pos
		return False
	else:
		return True


	
	breakpoints = []
	for i in range(0,len(breakpointArray)-1):
		#look at sets of reads for neighbouring positions
		readsLeft = breakpointArray[i]
		readsRight = breakpointArray[i+1]
		if not (len(readsLeft) == 0 or len(readsRight) == 0):
			s1 = set(readsLeft)
			s2 = set(readsRight)
			inter = s1.intersection(s2)
			if len(inter) == 0:
				breakpoints.append(str(gapStart+i))
	print "Breakpoints: "+",".join(breakpoints)
	print "-------------"





#read cigar and return list of start and stop tuple of Insertions in the read
def readcigarInsertions(cigar):
	insertions = []
	pos =0
	while cigar:
		#suppose the cigar starts with a match
		matchM = re.findall(r'^(\d+)M', cigar)
		if not matchM == []:
			pos = pos + int(matchM[0]) 
			cigar = cigar[len(matchM[0])+1:]
			continue
		
		
		#suppose the cigar stars with an insertion	
		matchI = re.findall(r'^(\d+)I', cigar)
		if not matchI == []:
			ins = (pos, pos+int(matchI[0]) -1)
			insertions.append(ins)
			pos = pos + int(matchI[0]) 
			cigar = cigar[len(matchI[0])+1:]
			continue

		#if we ever reach this point, something is wrong
		print "error: insertions not parsed correctly: "+cigar
		break
	
	return insertions


def readcigarDeletions(cigar,seq):
	pos=0
	newSeq = ""
	while cigar:
		#suppose the cigar starts with a match
		matchM = re.findall(r'^(\d+)M', cigar)
		if not matchM == []:
			newSeq = newSeq + seq[pos:pos+int(matchM[0])]
			pos = pos + int(matchM[0]) 
			cigar = cigar[len(matchM[0])+1:]
			continue
		
		#suppose the cigar stars with an insertion	
		matchD = re.findall(r'^(\d+)D', cigar)
		if not matchD == []:
			newSeq = newSeq + int(matchD[0])*"D"
			cigar = cigar[len(matchD[0])+1:]
			continue

		#if we ever reach this point, something is wrong
		print "error: insertions not parsed correctly"
		break
		
		
	return newSeq
	




#############################################################################################  
### read mappingFile and save each mapping in array, then sort array

file = open(mappingFile, "r")
readMappings = []
for line in file:
	if not line == "\n":
			fields = line.split("\t")
			id = fields[0]
			readstart = int(fields[3])
			cigar = fields[5]
			seq = fields[9]
			listIns = []
			#if not seq == "*" and not "S" in cigar and not "H" in cigar:
			if not seq == "*":
				#find softclipping at beginning of sequence
				matchS = re.findall(r'^(\d+)S', cigar)
				if not matchS == []:
					#remove clipped part from sequence, we do not need to change start position
					#of the mapping
					seq = seq[int(matchS[0]):]
					cigar = cigar[len(matchS[0])+1:]
				#find softclipping at end of sequence
				matchS2 = re.findall(r'(\d+)S$', cigar)
				if not matchS2 == []:
					seq = seq[:-int(matchS2[0])]
					cigar = cigar[:-(len(matchS2[0])+1)]
					
				#suppose there are just insertions in the mapping:
				if "I" in cigar and not "D" in cigar:
					matchM = re.findall(r'(\d+)M', cigar)
					#the insertion in the read to the reference does not influence the length of the mapping in the
					# reference, just the length of the read
					stop = int(readstart) + sum(int(x) for x in matchM)
					listIns = readcigarInsertions(cigar)
				#suppose there are deletions in the mapping:
				elif "D" in cigar and not "I" in cigar:
					matchM = re.findall(r'(\d+)M', cigar)
					matchD = re.findall(r'(\d+)D', cigar)
					stop = int(readstart) + sum(int(x) for x in matchM) + sum(int(y) for y in matchD)
					seq = readcigarDeletions(cigar,seq)
					if not stop - readstart == len(seq):
						print "error: deletions not parsed correctly"
				#there are no InDels in the mapping:
				elif not "I" in cigar and not "D" in cigar:
					matchM = re.findall(r'(\d+)M', cigar)
					stop = int(readstart) + int(matchM[0])
				elif "I" in cigar and "D" in cigar:
					print "Not handled read with mapping: "+cigar
					continue
				
				
				if gapStart <= readstart <= gapEnd or gapStart <= stop <= gapEnd or readstart <= gapStart <= gapEnd <= stop:					
					readMappings.append([id,readstart,stop,seq,cigar,listIns])


testBool = testCoverage(readMappings, gapStart, gapEnd)
print "Coverage test: "+str(testBool)


if len(readMappings) > 0 and testBool:
#if len(readMappings) > 0:

	sortedMappings = sorted(readMappings, key=itemgetter(1,2))
	ret = shortestPath_indels.constructGraph(sortedMappings,gapStart,gapEnd,template)
	graph = ret[0]
	first = ret[1]
	last = ret[2]
	lastNodeSeq = ret[3]

	path = shortestPath_indels.shortestPath(graph,first,last)
	sequence = shortestPath_indels.reconstructSequence(path,graph)

	if not sequence == "":	
		cut_index = shortestPath_indels.freeEndGap(sequence,lastNodeSeq)
		
		if cut_index == 0:
			cut_sequence = sequence[0:cut_index]
		else:
			cut_sequence = sequence[0:cut_index+1]
	else:
		cut_sequence = sequence
		

	distance = shortestPath_indels.editDistance(cut_sequence,template)
	rev_sequence = complement(cut_sequence)
	rev_distance = shortestPath_indels.editDistance(rev_sequence,template)
	
	#print sequence with orientation that gives the smaller edit distance to the template
	if rev_distance < distance:
		distance = rev_distance
		cut_sequence = rev_sequence
	cut_sequence = re.sub('D+', '', cut_sequence)
	print "total distance: "+str(distance)
	if len(sequence) > 0:
		print "relative distance: "+str(float(distance)/len(sequence))
	else:
		print "sequence of length 0"
	print "read sequence:\n"
	print cut_sequence
	print "template sequence:\n"
	print template
	print "------------------------------------------------"
	out = open(output,"w")
	out.write(">"+savedLine.rstrip("\n")+" "+str(distance)+"\n")
	out.write(cut_sequence+"\n")
	out.close()
else:
	print "NO COVERING READS"
	print "------------------------------------------------"
#	out = open(output+".uncov","w")
	out = open(output,"w")
	out.write("UNCOV"+"\n")
	out.write(">"+savedLine.rstrip("\n"))
	out.write("Covered: "+str(testBool)+"\n")
	out.close()


print time.time() - t0, "seconds process time"



		
				

		

			
