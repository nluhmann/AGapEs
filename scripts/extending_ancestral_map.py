import sys


## read ancestral map produced after scaffolding and augment it with gapFilling results
## each gap entry then contains the following information:
## gapID, adjacency, gap length in map, gap length covered by reads, consistent/conflicting/IS gap, dubious gap sequence flag OK/DUBIOUS
## gap sequences are marked as DUBIOUS if:
##		gap is completely not covered by reads (gap sequence is template sequence)
##		extant gap lengths are highly unconserved and differ more than IS length
##		regions in the gap sequence are only covered by single read
##		gap longer than X bp.
##		gap with high read sequence distance to template
##
##
##


#original ancestral map
ancestral_map_file = sys.argv[1]

#index file after full coverage gapFilling
index_file = sys.argv[2]

#index file after partial gapFilling consistent
partial_index_file = sys.argv[3]

#index file after partial gapFilling IS
partial_index_file_IS = sys.argv[4]

#adjacencies info file
info_file = sys.argv[5]





#read adjacencies info file and save status for each gap in hash
status_hash = {}
length_hash = {}

file = open(info_file,"r")
for line in file:
	if "gap" in line:
		ID = line.rstrip("\n")
		status = "CONSISTENT"
		length = []
	elif "GAPS" in line:
		arr = line.rstrip("\n").split(" ")
		length.append(arr[4])
	elif "IS" in line:
		arr = line.split(" ")
		if not arr[1] == "0":
			status = "IS"
	elif "CONN_COMP" in line:
		if not line.rstrip("\n")[-1] == "0":
			status = "CONFLICTING"
	
		status_hash[ID] = status
		length_hash[ID] = length
		
file.close()


ID_hash = {}
coverage_hash = {}
#read index file and save for all full covered gaps the length
file = open(index_file,"r")
for line in file:
	array = line.rstrip("\n").split("\t")
	ID_hash[array[1]] = array[0]
	spl=array[1].split("-")
	ID_hash[spl[1]+"-"+spl[0]] = array[0]
	if not array[2] == "FITCH" and not array[2] == "LEN_0":
		#length, covered length, dist to template
		coverage_hash[array[0]] = [array[4],array[4],array[3]]
	else:
		coverage_hash[array[0]] = ["0","0","-"]
file.close()


file = open(partial_index_file,"r")
for line in file:
	array = line.rstrip("\n").split("\t")
	if array[0] == "PARTIALLY":
		coverage_hash[array[1]] = [array[3],array[2],"-"]
		
	else:
		coverage_hash[array[1]] = [array[3],"0","-"]
file.close()

file = open(partial_index_file_IS,"r")
for line in file:
	array = line.rstrip("\n").split("\t")
	if array[0] == "PARTIALLY":
		coverage_hash[array[1]] = [array[3],array[2],"-"]
		
	else:
		coverage_hash[array[1]] = [array[3],"0","-"]
file.close()



#read map and augment
file = open(ancestral_map_file,"r")
out = open(ancestral_map_file+"_augmented","w")

for line in file:
	if "GAP" in line:
		array = line.rstrip("\n").split(" ")
		dubious = "OK"
		if len(array) == 3:
			#gap of length 0
			out.write(line)
		elif "MISSING_GAP" in line:
			out.write(line)
		else:
			gapID = ID_hash[array[1]]
			#check if dubious
			lengths = length_hash[gapID]
			cov_length = coverage_hash[gapID][1]
			if int(lengths[0]) > 3000:
				dubious = "DUB EXT_LENGTH_DIFF"
			elif int(array[2]) > 20000:
				dubious = "DUB LEN"
			elif not coverage_hash[gapID][2] == "-" and int(coverage_hash[gapID][2]) > 30:
				dubious = "DUB DIST"
			elif cov_length == "0":
				dubious = "DUB COV"
			status = status_hash[gapID]
			out.write("GAP"+" "+gapID+" "+array[1]+" "+array[3]+" "+array[4]+" "+array[2]+" "+cov_length+" "+status+" "+dubious+"\n")
	else:
		out.write(line)


file.close()
out.close()















