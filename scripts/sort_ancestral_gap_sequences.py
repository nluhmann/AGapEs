import sys
import os

#gaps_coordinates_and_length
gaps = sys.argv[1]

#finishing folder
fi = sys.argv[2]

#assembly folder
ass = sys.argv[3]

out = fi+"alignments_final/"

fitch = fi+"alignments/"

ass_cons = ass+"consistent/"

ass_conf = ass+"conflicts/"

ass_is = ass+"IS_gaps/"

ass_is_IDs = ass+"IS_gaps/gaps_filled"

#############################################################################################  	
### return sequence complement
baseComplement = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
def complement(sequence):
    result = []
    for n in range(len(sequence)):
        result.append(baseComplement[sequence[-(n+1)]])
    return ''.join(result)


#1 read all IS filled gapIDs, find gap seq in ass_is
IS_gaps = []
file = open(ass_is_IDs, "r")
for line in file:
	if not line == "\n":
		IS_gaps.append(line.rstrip("\n"))
file.close()


def check_gaps(firstpath, secondpath, fitchpath, gapid, IS_gaps):
	sequence = ""
	if os.path.isfile(firstpath) and os.path.isfile(secondpath):
		#read and check if they are the same
		seq1 = readSeq(firstpath)
		seq2 = readSeq(secondpath)
		if seq1 == seq2 or gapid in IS_gaps:
			sequence = seq1
			template = readSeq(fitchpath)
			dist = editDistance(sequence, template)
			if dist > 0:
					rev = complement(sequence)
					dist_rev = editDistance(rev, template)
					if dist_rev < dist:
						sequence = rev
						dist = dist_rev
			#print "distance "+str(dist)
			
		else:
			#take fitch template
			sequence = readSeq(fitchpath)
	else:
		if os.path.isfile(firstpath):
			sequence = readSeq(firstpath)
			template = readSeq(fitchpath)
			dist = editDistance(sequence, template)
			if dist > 0:
					rev = complement(sequence)
					dist_rev = editDistance(rev, template)
					if dist_rev < dist:
						sequence = rev
						dist = dist_rev
			#print "distance "+str(dist)
		else:
			sequence = readSeq(secondpath)
			template = readSeq(fitchpath)
			dist = editDistance(sequence, template)
			if dist > 0:
					rev = complement(sequence)
					dist_rev = editDistance(rev, template)
					if dist_rev < dist:
						sequence = rev
						dist = dist_rev
			#print "distance "+str(dist)	
	return sequence
	

def readSeq(path):
	seq = ""
	file = open(path, "r")
	for line in file:
		if not line.startswith(">"):
			seq = seq + line.rstrip("\n")
	file.close()
	return seq	
	
	
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



file = open(gaps, "r")
for line in file:
	if line.startswith(">"):
		gapid = line.split(" ")[0][1:]
		interval = line.rstrip("\n").split(" ")[-1]
		sequence = ""
		if not interval == "0-0":
			#check if we filled the sequence, otherwise take Fitch template
			if os.path.isfile(ass_cons+gapid+".out"):
				#print "assembly consistent"
				sequence = readSeq(ass_cons+gapid+".out")
				template = readSeq(fitch+gapid+".fasta_ancestral")
				if not sequence == "UNCOV":
				    dist = editDistance(sequence, template)
				    if dist > 0:
					rev = complement(sequence)
					dist_rev = editDistance(rev, template)
					if dist_rev < dist:
						sequence = rev
						dist = dist_rev
				    #print "distance "+str(dist)
				else:
					sequence = template
				
			elif os.path.isfile(ass_conf+gapid+"-1.out") or os.path.isfile(ass_conf+gapid+"-2.out") or os.path.isfile(ass_conf+gapid+".out"):
				#print "assembly conflicts"
				if os.path.isfile(ass_conf+gapid+".out"):
					sequence = readSeq(ass_conf+gapid+".out")
					template = readSeq(fitch+gapid+".fasta_ancestral")
					if not sequence == "UNCOV":
						dist = editDistance(sequence, template)
						if dist > 0:
							rev = complement(sequence)
							dist_rev = editDistance(rev, template)
							if dist_rev < dist:
								sequence = rev
								dist = dist_rev
						#print "distance "+str(dist)
					else:
						sequence = template
				else:
					# if only one version filled: take it
					# if both versions filled: if both filling the same: take it
					#							else: take fitch
					sequence = check_gaps(ass_conf+gapid+"-1.out",ass_conf+gapid+"-2.out",fitch+gapid+".fasta_ancestral", gapid, IS_gaps)
				
			elif os.path.isfile(ass_is+gapid+"-1.out") or os.path.isfile(ass_is+gapid+"-2.out"):
				#print "assembly is"
				# if only one version filled: take it
				# if both versions filled: if both filling the same: take it
				#							else: take fitch
				sequence = check_gaps(ass_is+gapid+"-1.out",ass_is+gapid+"-2.out",fitch+gapid+".fasta_ancestral",gapid, IS_gaps)
				
			elif os.path.isfile(fitch+gapid+".fasta_ancestral"):
				#print "fitch"
				sequence = readSeq(fitch+gapid+".fasta_ancestral")
				
			else:
				print "ERROR"
				print gapid
		outfile = open(out+gapid+".fasta_ancestral","w")
		outfile.write(sequence)
		outfile.close()

file.close()




