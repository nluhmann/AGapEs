#!/usr/bin/env python


import heapq
import itertools
import os
import numpy
import copy






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

	
#############################################################################################  
#compute free-end-gap alignment between two strings, return first match position in first string

def freeEndGap(first,second):
	#define DP matrix, first string vertical, second string horizontal
	#Initialization: vertical 0, horizontal sum of insertions
	#Maximum: last row
	#Backtracking: until first column, then return index
	#matrix = [[0 for x in range(len(first))] for x in range(len(second))] 	
	if len(second) > len(first):
		return 0
	matrix = numpy.zeros((len(first)+1,len(second)+1))
	
	#Initialization
	for i in range(len(second)+1):
		matrix[0][i] = i
	
	#row-wise iteration
	for i in range(1,len(first)+1):
		for j in range(1,len(second)+1):
			ins = matrix[i-1][j] + 1
			dele = matrix[i][j-1] +1
			if first[i-1] == second[j-1]:
				sub = matrix[i-1][j-1]
			else:
				sub = matrix[i-1][j-1] + 1
			matrix[i,j] = min(ins,dele,sub)

	

	
	min_index = len(second)
	
	# backtracing
	j = min_index
	i = -1
	while not j == 0:
		value = matrix[i,j]
		if first == "":
			j = 0
			i = 0
			break	
		
		#first check diagonal
		if first[i] == second[j-1]:
			predecessor_value = matrix[i-1][j-1]
		else:
			predecessor_value = matrix[i-1][j-1] +1
		if value == predecessor_value:
			i = i - 1
			j = j -1
			continue
		#check indels
		if value == matrix[i-1][j] +1:
			i = i-1
			continue
		elif value == matrix[i][j-1] +1:
			j = j-1
			continue
		else:
			print "error in here"
			break

	index = len(first) +i 
	return index



#test method
#freeEndGap("ABCD","EFG")
				
	
    
#############################################################################################
### Dijkstra algorithm    
### graph: dict keys=vertices, value= (dict keys=neighbor vertices, values=edgeLength)
### output: distance matrix, where distances[v] is the distance from start to v
###			predecessor matrix, where predecessors[v] is the predecessor of v in the shortest path from start to v.

def Dijkstra(graph,start,end):
	heap = []                         
	entry_finder = {}               
	REMOVED = '<removed-task>'      
	counter = itertools.count()     

	#dict of distances
	distances = {}	
	#dict of predecessors
	predecessors = {}	
	#priority queue
	heap = []
	(heap,entry_finder,counter) = add_node(heap,entry_finder,counter,start,0)
	
	distances[start] = 0
	predecessors[start] = ""
	
	
	#initialization
	for vertex in graph:
		if vertex != start:
			distances[vertex] = float("inf")
			predecessors[vertex] = ""
			(heap,entry_finder,counter) = add_node(heap,entry_finder,counter,vertex,distances[vertex])
	
			
	while len(entry_finder) > 0:
		(currentVertex,heap,entry_finder) = pop_node(heap,entry_finder,REMOVED)
		if currentVertex == end:
			break
		for neighbor in graph[currentVertex]:
			if neighbor in entry_finder:
				newDist = distances[currentVertex] + graph[currentVertex][neighbor][0]
				if newDist < distances[neighbor]:
					distances[neighbor] = newDist
					predecessors[neighbor] = currentVertex
					#update priority queue
					entry_finder = remove_node(neighbor,entry_finder,REMOVED)
					(heap,entry_finder,counter) = add_node(heap,entry_finder,counter,neighbor,newDist)
	
			
	return (distances,predecessors)
	
#############################################################################################  
### helping functions to insert entry into queue with priority or update priority of an existing entry, 
### remove an entry from the queue and
### return entry with lowest priority that is not marked as removed

def add_node(heap,entry_finder,counter,node,priority=0):
    'Add a new node or update the priority of an existing node'
    if node in entry_finder:
        remove_node(node)
    count = next(counter)
    entry = [priority, count, node]
    entry_finder[node] = entry
    heapq.heappush(heap, entry)
    return (heap,entry_finder,counter)

def remove_node(node,entry_finder,REMOVED):
    'Mark an existing node as REMOVED.  Raise KeyError if not found.'
    entry = entry_finder.pop(node)
    entry[-1] = REMOVED
    return entry_finder

def pop_node(heap,entry_finder,REMOVED):
    'Remove and return the lowest priority task. Raise KeyError if empty.'
    while heap:
        priority, count, node = heapq.heappop(heap)
        if node is not REMOVED:
            del entry_finder[node]
            return (node,heap,entry_finder)
    raise KeyError('pop from an empty priority queue')
	



#############################################################################################  
### compute shortest path from start to end
### graph: dict keys=vertices, value= (dict keys=neighbor vertices, values=edgeLength)
### get path from end to start in predecessors
			
def shortestPath(graph,start,end):
	print "Find shortest path..."
	(distances,predecessors) = Dijkstra(graph,start,end)
	path = []
	if not predecessors[end] == "":
		while 1:
			path.append(end)
			if end == start: break
			end = predecessors[end]
		path.reverse()
	else:
		print "gap not covered?"
	return path

#############################################################################################  
### pretty print graph


#############################################################################################  
### construct graph from all mapped reads, assume list to be sorted by start position of the mapping
### return graph together with first and last node for shortest path
### readList: [(readID,start,stop,sequence,cigar)]

def constructGraph(readList, gapStart, gapEnd, template):
	print "constructGraph"
	print "Number of reads: "+str(len(readList))+" for template of length: "+str(len(template))
	graph = {}
	
	if not gapEnd-gapStart + 1 == len(template):
		print "error"
		print gapEnd - gapStart 
		print len(template)
	
	#take first read in the list as startRead, crop after gapStart
	startRead = copy.deepcopy(readList[0])
	startRead[0] = "start"
	overlap = startRead[2] - gapStart
	
	startRead[2] = gapStart 
	sequence = startRead[3] 
	if not overlap == 0:
		startRead[3] = sequence[0:-overlap]
	readList.insert(0, startRead)
	
	#readList[0] = startRead
	#take last read in the list as endRead, crop before gapEnd
	endRead = copy.deepcopy(readList[-1])
	endRead[0] = "end"
	overlap = gapEnd - endRead[1] + 1
	counter =0
	for ins in endRead[5]:
		if ins[0] <= overlap +1:
			overlap = overlap + ins[1]-ins[0]+1
			counter = counter + ins[1]-ins[0]+1
	endRead[1] = gapEnd + counter + 1
	sequence = endRead[3]
	endRead[3] = sequence[overlap:]
	#readList[-1] = endRead
	readList.append(endRead)
	
	#consider first read on pos i in readList, consider second read on pos j in readList
	for i in range(len(readList)):
		inserted = False
		for j in range(i+1,len(readList)):
			firstEnd = readList[i][2]
			secondStart = readList[j][1]
			secondEnd = readList[j][2]
			secondIndels = readList[j][5]
			
			#no need to add an edge if no overlap between reads
			if secondStart > firstEnd:
				break
				
			#overlap is needed to define suffix of second read, hence it should respect insertions
			overlap = firstEnd - secondStart 
			for ins in secondIndels:
				if ins[1] <= overlap:
					overlap = overlap + ins[1]-ins[0]+1
			
			#compute edit Distance between suffix of second read and template sequence
			secondSuffixStart = firstEnd - gapStart
			secondSuffixEnd = secondEnd - gapStart
			templateSeq = template[secondSuffixStart:secondSuffixEnd]
			secondSuffixSeq = readList[j][3][overlap:]
			
			distance = editDistance(secondSuffixSeq, templateSeq)
			
			
			#insert first read and directed edge to second read with distance into graph
			if readList[i][0] in graph:
				# if readList[j][0] in graph[readList[i][0]]:
# 					print readList[i][0]
# 					print readList[j][0]
# 					print "error, edge already inserted!"
# 				else:
					graph[readList[i][0]][readList[j][0]] = (distance,secondSuffixSeq)	
					inserted = True	
			else:
				graph[readList[i][0]] = {}
				graph[readList[i][0]][readList[j][0]] = (distance,secondSuffixSeq)
				inserted = True
		if not inserted:
			graph[readList[i][0]] = {}
	
	#insert last read into graph that has no neighboring edges
	graph[readList[-1][0]] = {}
	
	return (graph,readList[0][0],readList[-1][0], readList[-1][3])


#############################################################################################  
### translate shortest path into genome sequence based on extensions of overlapping read mappings

def reconstructSequence(path,graph):
	sequence = ""
	for i in range(len(path)-1):
		node = path[i]
		next = path[i+1]
		extension = graph[node][next][1]
		sequence = sequence + extension
		
	
	return sequence





#############################################################################################  	
### return sequence complement
baseComplement = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
def complement(sequence):
    result = []
    for n in range(len(sequence)):
        result.append(baseComplement[sequence[-(n+1)]])
    return ''.join(result)


	

# ############################
# #TEST
# 
# testtemplate = "AGCTAGCTAGCTAGCT"
# r1=["r1",152,156,"ATTAG",""]
# r2=["r2",154,161,"TAGTTAGC",""]
# r3=["r3",158,164,"TAGGTAG",""]
# r4=["r4",160,168,"GCTAGCTAG",""]
# r5=["r5",165,172,"CTAGCTAT",""]
# r6=["r6",168,174,"GCTATCG",""]
# testreads=[r1,r2,r3,r4,r5,r6]
# testgapStart = 154
# testgapStop = 171
# 
# testgraph = constructGraph(testreads,testgapStart,testgapStop,testtemplate)
# testpath = shortestPath(graph,"r1","r6")
# print path
# testsequence = reconstructSequence(path,graph)
# print sequence


# G = {'s':{'u':10, 'x':5}, 'u':{'v':1, 'x':2}, 'v':{'y':4}, 'x':{'u':3, 'v':9, 'y':2}, 'y':{'s':7, 'v':6}}
# path = shortestPath(G,"s","v")
	


	
	
	
