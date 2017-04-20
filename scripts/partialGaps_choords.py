import sys

gaps = sys.argv[1]

uncovered = sys.argv[2]

outfile = sys.argv[3]


gapHash = {}
coords = {}
file = open(gaps,"r")
for line in file:
	array = line.rstrip("\n").split(" ")
	gapID = array[0]
	gapHash[gapID] = [int(array[1]),int(array[2])]	
	coords[gapID] = []
file.close()



prev_gapID = ""
left = True
right = True
file = open(uncovered,"r")
for line in file:
	array = line.rstrip("\n").split("\t")
	gapID = array[0]
	if not gapID == prev_gapID and not prev_gapID == "":
		if left:
			coords[prev_gapID].append(gapHash[prev_gapID][0])
		if right:
			coords[prev_gapID].append(gapHash[prev_gapID][1])		
		left = True
		right = True
	prev_gapID = gapID		
	
	if gapID in gapHash:
		#we have to check if it also contains uncovered regions outside of gaps and determine filling mode.
		if int(array[2]) == gapHash[gapID][1]:
			right = False
			coords[gapID].append(int(array[1])+1)
		elif int(array[1])+1 == gapHash[gapID][0]:
			left = False
			coords[gapID].append(int(array[2]))
		elif int(array[1])+1 < gapHash[gapID][0] and int(array[2]) > gapHash[gapID][0]:
			#there is no left
			left = False
			coords[gapID].append(int(array[2]))
		elif int(array[1])+1 < gapHash[gapID][1] and int(array[2]) > gapHash[gapID][1]:
			#there is no right	
			right = False
			coords[gapID].append(int(array[1])+1)
		elif not int(array[1])+1 <= gapHash[gapID][0] and not int(array[2]) >= gapHash[gapID][1]:
			coords[gapID].append(int(array[1])+1)
			coords[gapID].append(int(array[2]))
if left:
	coords[gapID].append(gapHash[gapID][0])
if right:
	coords[gapID].append(gapHash[gapID][1])
	
		
file.close()




out = open(outfile,"w")
for key in gapHash:
	pos = sorted(coords[key])
	if not len(pos)%2 == 0:
		print "ERROR"
		print coords[key]
		print gapHash[key]
	else:
		
		for i in xrange(0,len(pos),2):
  			first = pos[i]
  			second = pos[i+1]
  			if first == gapHash[key][0]:
  				state = "LEFT" 			
  			elif second == gapHash[key][1]:
  				state = "RIGHT"
			else: 
				state = "INTERNAL"
			out.write(key+"\t"+str(first)+"-"+str(second)+"\t"+state+"\n")
out.close()
