import sys


uncovered = sys.argv[1]


file = open(uncovered,"r")
old_gapID = ""
old_stop = ""
old_start = ""
cov_count = 0
covs = []
for line in file:
	array = line.rstrip("\n").split("\t")
	gapID = array[0]
	start = array[1]
	stop = array[2]
	cov = array[3]
	

	if gapID == old_gapID and start == old_stop:
		old_stop = stop
	elif not old_gapID == "":
		
		print old_gapID+"\t"+old_start+"\t"+old_stop+"\t"+str(covs[0])+"\t"+str(covs[-1])
		old_start = start
		old_stop = stop
		covs = []
		cov_count = 0	
	else:
		old_start = start
		old_stop = stop
	old_gapID = gapID
	if cov == "1":
		cov_count = int(stop) - int(start)
		covs.append(cov_count)
		cov_count = 0	
	else:
		covs.append(0)
	
covs.append(cov_count)
print old_gapID+"\t"+old_start+"\t"+old_stop+"\t"+str(covs[0])+"\t"+str(covs[-1])		
file.close()