import sys


# read adjacency Info file, for each gap containing an IS, do:
# parse genomes containing an IS, then read extant.fasta for the respective gap and sort them into the two categories



adjacency_info_file = sys.argv[1]
alignments = sys.argv[2]


file = open(adjacency_info_file,"r")
IS_genomes = []
for line in file:
	if line == "\n":
		#collected all the infos, so go for it
		if IS_genomes:
			extant_genomes = open(alignments+gap+".fasta","r")
			extant_IS_list = []
			extant_NOIS_list = []
			
			first = True
			for line in extant_genomes:
				if line.startswith(">"):
					
					if not first:
						if species in IS_genomes and not "pseudotuberculosis" in species:
							extant_IS_list.append(header+sequence+"\n")
						elif not "pseudotuberculosis" in species:
							extant_NOIS_list.append(header+sequence+"\n")
						#else:
							#print species
					first = False
					header = line
					species = line.split(":")[0][1:].split(".")[0]
					sequence = ""
				else:
					sequence = sequence + line.rstrip("\n")
			if species in IS_genomes and not "pseudotuberculosis" in species:
				extant_IS_list.append(header+sequence+"\n")
			elif not "pseudotuberculosis" in species:
				extant_NOIS_list.append(header+sequence+"\n")
			extant_genomes.close()

			if not len(extant_IS_list) == 0:
				extant_IS = open(alignments+gap+"-2"+".fasta_extant","w")
				extant_IS.write(''.join(extant_IS_list))
				extant_IS.close()

			if not len(extant_NOIS_list) == 0:
				extant_NOIS = open(alignments+gap+"-1"+".fasta_extant","w")
				extant_NOIS.write(''.join(extant_NOIS_list))
				extant_NOIS.close()
			
			
	
	
		IS_genomes = []
	elif "gap_" in line:
		gap = line.rstrip("\n")
	
	elif "I_IS" in line or "O_IS" in line:
		line = line.rstrip("\n")
		if not len(line) == 8:
			array = line.split(" ")
			for i in range(2,len(array)-1):
				IS_genomes.append(array[i].split(":")[0])


file.close
