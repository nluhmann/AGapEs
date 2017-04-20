#!/usr/bin/env python

import sys

#adjacency information files, based on reference occurences before linearization, considering all adjacencies

ingroup = ["Yersinia_pestis_Antiqua","Yersinia_pestis_KIM_10","Yersinia_pestis_Nepal516","Yersinia_pestis_CO92","Yersinia_pestis_Z176003"]
outgroup = ["Yersinia_pseudotuberculosis_IP_32953","Yersinia_pseudotuberculosis_YPIII","Yersinia_pseudotuberculosis_IP_31758","Yersinia_pseudotuberculosis_PB1","Yersinia_pestis_biovar_Microtus_str_91001","Yersinia_pestis_Pestoides_F","Yersinia_pestis_Angola"]


adjacencies = sys.argv[1]
is_coords = sys.argv[2]
gap_coords = sys.argv[3]
outfile = sys.argv[4]
conflictsoutfile = sys.argv[5]
nonConservedoutfile = sys.argv[6]
conservedoutfile = sys.argv[7]
ISgapsoutfile = sys.argv[8]


#key: genome, value: (start, stop, element)
coordsHash = {}
#read is_coords
coords = open(is_coords, "r")
for line in coords:
	fields = line[1:].split(":")
	genome = fields[0].split(".")[0]
	startc = int(fields[1].split("-")[0])
	stopc = int(fields[1].split(" ")[0].split("-")[1])
	element = fields[1].split(" ")[1].rstrip("\n")
	if genome in coordsHash:
		coordsHash[genome].append((startc,stopc,element))
	else:
		coordsHash[genome] = [(startc,stopc,element)]
coords.close()



coords = open(gap_coords,"r")
gapcoordsHash = {}
for line in coords:
	if line.startswith(">"):
		array = line.rstrip("\n").split(" ")
		adj = array[2]
		left = adj.split("-")[0]
		right = adj.split("-")[1]
		gapcoordsHash[(left,right)] = array[0][1:]
		gapcoordsHash[(right,left)] = array[0][1:]
coords.close()






adjacencyHash = {}



#read number and list of adj occ in ing and outg, compute min max gap length and difference
disc = open(adjacencies,"r")
for line in disc:
	fields = line.split("#")
	adjacencyID = fields[0].split("|")[0]
	firstMarker = fields[0].split(":")[1].split(" ")[0]
	secondMarker =	fields[0].split(":")[1].split(" ")[1]
	smallestIn = 50000000000000
	smallestOut = 50000000000000
	largestIn = 0
	largestOut = 0
	ingroupList = []
	outgroupList = []
	ingroupIS = []
	outgroupIS = []
	for field in fields[1:]:
			coords = field.split(":")[1][:-2]
			genome = field.split(":")[0].split(".")[0]
			start = int(coords.split("-")[0])
			stop = int(coords.split("-")[1])
			length = stop - start +1
			
			if genome in ingroup:
				ingroupList.append(field.rstrip("\n"))
				if length < smallestIn:
					smallestIn = length
				if length > largestIn:
					largestIn = length
			elif genome in outgroup:
				outgroupList.append(field.rstrip("\n"))
				if length < smallestOut:
					smallestOut = length
				if length > largestOut:
					largestOut = length
			else:
				print "BUUUUUG"
			
			if genome in coordsHash:
				liste = coordsHash[genome]
				for elem in liste:
					IS = False
					if elem[0] > start and elem[1] < stop:
						IS = True
					elif elem[1] > start and elem[0] < start and elem[1] < stop:
						IS = True
					elif elem[0] < stop and elem[1]	> stop and elem[0] > start:
						IS = True
					
					if IS:
						if genome in ingroup:
							ingroupIS.append((elem,genome))
						elif genome in outgroup:
							outgroupIS.append((elem,genome))
						else:
							x = 1
							#print "BUUUUUG"
			
	#difference = largest - smallest
	adjacencyHash[adjacencyID] = [firstMarker,secondMarker,len(ingroupList),ingroupList,len(outgroupList),outgroupList,smallestIn,largestIn,largestIn-smallestIn,
	smallestOut,largestOut,largestOut-smallestOut,len(ingroupIS),ingroupIS,len(outgroupIS),outgroupIS]

disc.close()


conComMarker = [[]]
conComIDs = [[]]

file = open(adjacencies,"r")
for line in file:
	fields = line.split("#")
	adjacencyID = fields[0].split("|")[0]
	firstMarker = int(fields[0].split(":")[1].split(" ")[0])
	secondMarker =	int(fields[0].split(":")[1].split(" ")[1])
	already = False
	
	for i in range(len(conComMarker)):
		if firstMarker in conComMarker[i] or secondMarker in conComMarker[i]:
			if not already:
				#add marker to already existing connected component
				if not i == 0:
					conComMarker[i].append(firstMarker)
					conComMarker[i].append(secondMarker)
					conComIDs[i].append(adjacencyID)
					already = True
				#take marker out of non-conflicting list, create new connected component
				else:
					
					if not len(conComMarker[0]) == 0:
						conComMarker.append([firstMarker,secondMarker])
						conComIDs.append([adjacencyID])
						already = True
						
						
						if firstMarker in conComMarker[i]:
							index = conComMarker[0].index(firstMarker)
							if index % 2 == 0:
								conComMarker[-1].append(conComMarker[0][index+1])
								conComMarker[-1].append(firstMarker)
								conComIDs[-1].append(conComIDs[0][index/2])
								conComMarker[0].pop(index+1)
								conComIDs[0].pop(index/2)
							else:
								conComMarker[-1].append(conComMarker[0][index-1])
								conComMarker[-1].append(firstMarker)
								conComIDs[-1].append(conComIDs[0][(index-1)/2])
								conComMarker[0].pop(index-1)
								conComIDs[0].pop((index-1)/2)
							conComMarker[0].remove(firstMarker)
							
							
						if secondMarker in conComMarker[i]:
							index = conComMarker[0].index(secondMarker)
							if index % 2 == 0:
								conComMarker[-1].append(conComMarker[0][index+1])
								conComMarker[-1].append(secondMarker)
								conComIDs[-1].append(conComIDs[0][index/2])
								conComMarker[0].pop(index+1)
								conComIDs[0].pop(index/2)
							else:
								conComMarker[-1].append(conComMarker[0][index-1])
								conComMarker[-1].append(secondMarker)
								conComIDs[-1].append(conComIDs[0][(index-1)/2])
								conComMarker[0].pop(index-1)
								conComIDs[0].pop((index-1)/2)
							conComMarker[0].remove(secondMarker)
							
						
					else:
						
						conComMarker[0].append(firstMarker)
						conComMarker[0].append(secondMarker)
						conComIDs[0].append(adjacencyID)
						already = True
			else:
				"there are overlapping connected components!"
		
	if not already:
		
		conComMarker[0].append(firstMarker)
		conComMarker[0].append(secondMarker)
		conComIDs[0].append(adjacencyID)
		already = True


#sanitycheck:
for i in range(1,(len(conComMarker))):
	if not set(conComMarker[i]).isdisjoint(set(conComMarker[0])):
		print set(conComMarker[i]).intersection(set(conComMarker[0]))



#merge overlapping connected components
for i in range(len(conComMarker)):
	for j in range(i+1,len(conComMarker)):
		if not set(conComMarker[i]).isdisjoint(set(conComMarker[j])):
			conComMarker[i].extend(conComMarker[j])
			conComMarker[j] = []
			conComIDs[i].extend(conComIDs[j])
			conComIDs[j] = []
			
for i in range(len(conComIDs)):
	for j in range(i+1,len(conComIDs)):
		if not set(conComIDs[i]).isdisjoint(set(conComIDs[j])):
			print "error"

file.close()

for i in range (len(conComMarker)):
	#no connected component
	if i == 0:
		for elem in conComIDs[0]:
			adjacencyHash[elem].append((0,0,0))
	else:
		marker = conComMarker[i]
		ids = conComIDs[i]
		vertices = len(set(marker))
		edges = len(marker)/2
		for id in ids:
			adjacencyHash[id].append((i,vertices,edges))



out = open(outfile,"w")
for elem in adjacencyHash:
	liste = adjacencyHash[elem]
	out.write(">"+elem+" "+liste[0]+" "+liste[1]+"\n")
	out.write(gapcoordsHash[(liste[0],liste[1])]+"\n")
	out.write("I_OCC "+str(liste[2])+"\n")
	for line in liste[3]:
		out.write(line+"\n")
	out.write("O_OCC "+str(liste[4])+"\n")
	for line in liste[5]:
		out.write(line+"\n")
	out.write("I_GAPS (smallest|largest|difference): "+str(liste[6])+" "+str(liste[7])+" "+str(liste[8])+"\n")
	out.write("O_GAPS (smallest|largest|difference): "+str(liste[9])+" "+str(liste[10])+" "+str(liste[11])+"\n")
	out.write("I_IS: " + str(liste[12])+" ")
	for elem in liste[13]:
		out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
	out.write("\n")
	out.write("O_IS: " + str(liste[14])+" ")
	for elem in liste[15]:
		out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
	out.write("\n")
	out.write("CONN_COMP (id|vertices|edges): "+str(liste[16][0])+" "+str(liste[16][1])+" "+str(liste[16][2])+"\n")
	out.write("\n")
out.close()	


sortedlist = sorted(adjacencyHash.items(), key=lambda x:x[1][16])
out = open(conflictsoutfile,"w")

for elem in sortedlist:
	liste = elem[1]
	if not liste[16][0] == 0:
		out.write(">"+elem[0]+" "+liste[0]+" "+liste[1]+"\n")
		out.write(gapcoordsHash[(liste[0],liste[1])]+"\n")
		out.write("I_OCC "+str(liste[2])+"\n")
		for line in liste[3]:
			out.write(line+"\n")
		out.write("O_OCC "+str(liste[4])+"\n")
		for line in liste[5]:
			out.write(line+"\n")
		out.write("I_GAPS (smallest|largest|difference): "+str(liste[6])+" "+str(liste[7])+" "+str(liste[8])+"\n")
		out.write("O_GAPS (smallest|largest|difference): "+str(liste[9])+" "+str(liste[10])+" "+str(liste[11])+"\n")
		out.write("I_IS: " + str(liste[12])+" ")
		for elem in liste[13]:
			out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("O_IS: " + str(liste[14])+" ")
		for elem in liste[15]:
	 		out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("CONN_COMP (id|vertices|edges): "+str(liste[16][0])+" "+str(liste[16][1])+" "+str(liste[16][2])+"\n")
		out.write("\n")
out.close()	


out = open(nonConservedoutfile,"w")

for elem in sortedlist:
	liste = elem[1]
	if liste[8] != 0 or liste[11] != 0:
		out.write(">"+elem[0]+" "+liste[0]+" "+liste[1]+"\n")
		out.write(gapcoordsHash[(liste[0],liste[1])]+"\n")
		out.write("I_OCC "+str(liste[2])+"\n")
		for line in liste[3]:
			out.write(line+"\n")
		out.write("O_OCC "+str(liste[4])+"\n")
		for line in liste[5]:
			out.write(line+"\n")
		out.write("I_GAPS (smallest|largest|difference): "+str(liste[6])+" "+str(liste[7])+" "+str(liste[8])+"\n")
		out.write("O_GAPS (smallest|largest|difference): "+str(liste[9])+" "+str(liste[10])+" "+str(liste[11])+"\n")
		out.write("I_IS: " + str(liste[12])+" ")
		for elem in liste[13]:
			out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("O_IS: " + str(liste[14])+" ")
		for elem in liste[15]:
	 		out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("CONN_COMP (id|vertices|edges): "+str(liste[16][0])+" "+str(liste[16][1])+" "+str(liste[16][2])+"\n")
		out.write("\n")
out.close()	


out = open(conservedoutfile,"w")

for elem in sortedlist:
	liste = elem[1]
	if liste[8] == 0 and liste[11] == 0:
		out.write(">"+elem[0]+" "+liste[0]+" "+liste[1]+"\n")
		out.write(gapcoordsHash[(liste[0],liste[1])]+"\n")
		out.write("I_OCC "+str(liste[2])+"\n")
		for line in liste[3]:
			out.write(line+"\n")
		out.write("O_OCC "+str(liste[4])+"\n")
		for line in liste[5]:
			out.write(line+"\n")
		out.write("I_GAPS (smallest|largest|difference): "+str(liste[6])+" "+str(liste[7])+" "+str(liste[8])+"\n")
		out.write("O_GAPS (smallest|largest|difference): "+str(liste[9])+" "+str(liste[10])+" "+str(liste[11])+"\n")
		out.write("I_IS: " + str(liste[12])+" ")
		for elem in liste[13]:
			out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("O_IS: " + str(liste[14])+" ")
		for elem in liste[15]:
	 		out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("CONN_COMP (id|vertices|edges): "+str(liste[16][0])+" "+str(liste[16][1])+" "+str(liste[16][2])+"\n")
		out.write("\n")
out.close()	

out = open(ISgapsoutfile,"w")

for elem in sortedlist:
	liste = elem[1]
	if len(liste[13]) != 0 or len(liste[15]) != 0:
		out.write(">"+elem[0]+" "+liste[0]+" "+liste[1]+"\n")
		out.write(gapcoordsHash[(liste[0],liste[1])]+"\n")
		out.write("I_OCC "+str(liste[2])+"\n")
		for line in liste[3]:
			out.write(line+"\n")
		out.write("O_OCC "+str(liste[4])+"\n")
		for line in liste[5]:
			out.write(line+"\n")
		out.write("I_GAPS (smallest|largest|difference): "+str(liste[6])+" "+str(liste[7])+" "+str(liste[8])+"\n")
		out.write("O_GAPS (smallest|largest|difference): "+str(liste[9])+" "+str(liste[10])+" "+str(liste[11])+"\n")
		out.write("I_IS: " + str(liste[12])+" ")
		for elem in liste[13]:
			out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("O_IS: " + str(liste[14])+" ")
		for elem in liste[15]:
	 		out.write(str(elem[1])+":"+str(elem[0][2])+":"+str(elem[0][0])+"-"+str(elem[0][1])+" ")
		out.write("\n")
		out.write("CONN_COMP (id|vertices|edges): "+str(liste[16][0])+" "+str(liste[16][1])+" "+str(liste[16][2])+"\n")
		out.write("\n")
out.close()		
