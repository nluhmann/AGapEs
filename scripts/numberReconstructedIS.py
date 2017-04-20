import sys
import math

IS_info = sys.argv[1]

#gaps extracted from augmented map
IS_gaps = sys.argv[2]


lengths = {}
annotatedspecs_gap = {}

outgroup = ["Yersinia_pseudotuberculosis_IP_31758","Yersinia_pseudotuberculosis_IP_32953","Yersinia_pseudotuberculosis_PB1","Yersinia_pseudotuberculosis_YPIII","Yersinia_pestis_biovar_Microtus_str_91001","Yersinia_pestis_Pestoides_F"]
ingroup = ["Yersinia_pestis_CO92","Yersinia_pestis_Z176003","Yersinia_pestis_Antiqua","Yersinia_pestis_Nepal516","Yersinia_pestis_KIM_10"]

co = ["Yersinia_pestis_CO92","Yersinia_pestis_Z176003","Yersinia_pestis_Antiqua"]
kim = ["Yersinia_pestis_Nepal516","Yersinia_pestis_KIM_10"]


dollo_ancestral = 0
dollo_ids = set()

file = open(IS_info, "r")
for line in file:
	if "gap_" in line:
		gapID = line.rstrip("\n")
		all_extant = {}
		annotated_specs = set()
		annotatedspecs_gap[gapID] = []
	elif "_IS" in line:
		#Is annotations
		array = line.rstrip("\n").split(" ")
		annos = array[2:]
		if len(annos) == 0:
			print "ERROR"				
		for elem in annos:
				if not elem == "":
					annotated_specs.add(elem.split(":")[0])
					coords = elem.split(":")[-1]
					start = int(coords.split("-")[0])
					stop = int(coords.split("-")[1])
					length = abs(stop - start)
					annotatedspecs_gap[gapID].append(length)
	elif "Yersinia" in line:
		spec = line.split(":")[0][:-4]
		coords = line.rstrip("\n").split(":")[1][:-2]
		start = int(coords.split("-")[0])
		stop = int(coords.split("-")[1])
		length = stop - start
		all_extant[spec]=length
	
	elif line == "\n":
		#compute average length!
		avg_IS_list = []
		avg_noIS_list = []
		for spec in all_extant:
			if spec in annotated_specs:
				avg_IS_list.append(all_extant[spec])
			else:
				avg_noIS_list.append(all_extant[spec])

		avg_IS = sum(avg_IS_list) / float(len(avg_IS_list))
		if not len(avg_noIS_list) == 0:
			avg_noIS = sum(avg_noIS_list) / float(len(avg_noIS_list))
		else:
			avg_noIS = "-"
		lengths[gapID] = [avg_IS,avg_noIS]
		
		ingroup_bool = False
		outgroup_bool = False
		for elem in annotated_specs:
			if elem in outgroup:
				outgroup_bool = True
			elif elem in ingroup:
				ingroup_bool = True
		if ingroup_bool and outgroup_bool:
			dollo_ancestral += 1	
			dollo_ids.add(gapID)
		elif ingroup_bool:
			co_bool = False
			kim_bool = False
			for elem in annotated_specs:
				if elem in co:
					co_bool = True
				elif elem in kim:
					kim_bool = True
			if co_bool and kim_bool:
				dollo_ancestral += 1
				dollo_ids.add(gapID)
		
		
file.close()

counter_noIS = 0
counter_IS = 0
IS_ids = set()
limitset = set()


file =open(IS_gaps, "r")
for line in file:
	array = line.split(" ")
	id = array[1]
	length = int(array[5])
	length_diff_IS = abs(lengths[id][0]-length)
	
	if not lengths[id][1] == "-":
		length_diff_noIS = abs(lengths[id][1]-length)
	else:
		length_diff_noIS = 10000000000
	if length_diff_IS >= length_diff_noIS:
		counter_noIS += 1
	else:
		counter_IS += 1
		IS_ids.add(id)
		for elem in annotatedspecs_gap[id]:
			if elem < 400:
				limitset.add(id)
file.close()

print "no IS: "+str(counter_noIS)
print "IS: "+str(counter_IS)
print "Limit: "+str(len(limitset))
print "Dollo: "+str(dollo_ancestral)

diff = IS_ids.symmetric_difference(dollo_ids)
print IS_ids
print diff
