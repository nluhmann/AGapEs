import sys


#ISgaps.info
info = sys.argv[1]

#grep from all bedgraphs
uncovered = sys.argv[2]

#relative, grep from fasta_assembly
gapCoords = sys.argv[3]

gapsBoth = sys.argv[4]


relHash = {}
file = open(info,"r")
for line in file:
	if "gap_" in line:
		gap = line.rstrip("\n")
		extantHash = {}
		ISHash = {}
	elif "IS" in line:
		
		all = line.rstrip("\n").split(" ")[2:-1]
		for elem in all:
			
			species = elem.split(":")[0]
			coords = elem.split(":")[2].split("-")
			start = int(coords[0])
			stop = int(coords[1])
			if species in ISHash:
				prevcoords = ISHash[species]
				if start < int(prevcoords[0]):
					prevcoords[0] = start
				if stop > int(prevcoords[1]):
					prevcoords[1] = stop
			else:
				ISHash[species] = coords
	
	
	elif "pestis" in line:
		species = line.split(":")[0].split(".")[0]
		coords = line.split(" ")[0].split(":")[1].split("-")
		start = coords[0]
		stop = coords[1]
		extantHash[species]=coords
	elif "CONN_COMP" in line: 
		startinterval = []
		stopinterval = []
		for spec in ISHash:
			extantCoords = extantHash[spec]
			IScoords = ISHash[spec]
			relStart = int(IScoords[0])-int(extantCoords[0])
			startinterval.append(relStart)
			relStop = int(IScoords[1])-int(extantCoords[0])
			stopinterval.append(relStop)
		
		minStart = min(startinterval)
		maxStart = max(startinterval)
		
		minEnd = min(stopinterval)
		maxEnd = max(stopinterval)
		
		
		#print gap+"\t"+str(minStart)+"\t"+str(maxStart)+"\t"+str(minEnd)+"\t"+str(maxEnd)
		relHash[gap] = [minStart,maxStart,minEnd,maxEnd]

file.close()


gapCoordsHash = {}
file = open(gapCoords,"r")

for line in file:
	gap = line.split(":")[0].split("/")[-1].split(".")[0]
	coords = line.rstrip("\n").split(" ")[-1].split("-")
	gapCoordsHash[gap]=coords

file.close()


gapsBothCovered = []
file = open(gapsBoth,"r")
for line in file:
	gap = line.rstrip("\n")
	gapsBothCovered.append(gap)
file.close()





def checkRange(coords,start,stop):
	coordsRangeStart = int(start) -20
	coordsRangeStop = int(stop) + 20
	
	if coordsRangeStart <= int(coords[0]) <= coordsRangeStop:
		return True
	elif coordsRangeStart <= int(coords[1]) <= coordsRangeStop:
		return True
	else:
		return False






#IS coords are relative to gap, while uncovered bedgraph coords are relative to whole fasta_assembly
#gapCoordsHash contains start and stop of gap relative to whole fasta_assembly


counter = 0
file = open(uncovered,"r")
for line in file:
	gap = line.split(":")[0].split(".")[0].split("/")[-1]
	gapID = gap.split("-")[0]
	if gapID in gapsBothCovered:
		start = line.split("\t")[1]
		stop = line.split("\t")[2]
		cov = line.rstrip("\n").split("\t")[3]
		gapCoords = gapCoordsHash[gap]
		if cov == "0" or cov == "1":
			if "-1" in gap:
				#print "uncovered in non-IS template"
		
				#check gap/marker border
				
				inRange = checkRange(gapCoords,start,stop)
				if inRange:
	
					print gap
					print cov
					print start
					print stop
					print gapCoords
					print "---------"
		
			elif "-2" in gap:
			
		
				#check gap/marker border
				gapCoords = gapCoordsHash[gap]
				inRange = checkRange(gapCoords,start,stop)
				if inRange:
	
					print gap
					print cov
					print start
					print stop
					print gapCoords
					print "---------"


		
		
		
			else:
				print "ERROR"
file.close














