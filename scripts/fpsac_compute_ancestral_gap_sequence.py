# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing the ancestral sequence of a gap  using Fitch algorithm

# argument 1 = Input  multiple alignment file
# argument 2 = Input  phylogenetic tree with marked ancestor
# argument 3 = Output ancestral gap FASTA file

import sys,tree,string

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

# Auxiliary procedure that takes a phylogenetic tree with a specified ancestor
# and a set of discrete caracters and produces the ancestral character
# according to the Fitch procedure.
# characters format: dictionary, keys are current species and values
# are any state (taken as a string)

def Fitch(species_tree,characters):
	########################################################
	leaves = tree.getLeaves(species_tree,tree.getRoot(species_tree))
	extant_species = []
	for l in leaves:
	    spe = tree.getName(species_tree,l)
	    extant_species.append(spe)
	ancestral_species = []
	for node in tree.getNodes(species_tree):
	    if not tree.isLeaf(species_tree,node):
		leaves = tree.getLeaves(species_tree,node)
		species = []
		for l in leaves:
		    species.append(tree.getName(species_tree,l))
		species.sort()
		ancestral_species.append([tree.getName(species_tree,node),species])
	########################################################
	postorder = tree.getNodes(species_tree)
	index = {}
	for p in range(len(postorder)):
	    index[postorder[p]] = p
	isapostorder = False
	while not isapostorder:
	    isapostorder = True
	    for p in range(len(postorder)):
		if not tree.isRoot(species_tree,postorder[p]):
		    i = index[tree.getParent(species_tree,postorder[p])]
		    if p > i:
			isapostorder = False
			index[postorder[p]] = i
			index[postorder[i]] = p
			temp = postorder[p]
			postorder[p] = postorder[i]
			postorder[i] = temp
	########################################################
	value = {}
	for current in postorder:
		if tree.isLeaf(species_tree,current):
			spe = tree.getName(species_tree,current)
			if characters.has_key(spe):
				value[current] = [characters[spe]]
			else:
				value[current] = []
		else:
			value0 = value[tree.getChildren(species_tree,current)[0]]
			value1 = value[tree.getChildren(species_tree,current)[1]]
			inter = []
			for v in value0:
				if v in value1:
					inter.append(v)
			if inter == []:
				value[current] = value0 + value1
			else:
				value[current] = inter
	postorder.reverse()
	for current in postorder:
		if len(value[current]) == 1 or tree.isRoot(species_tree,current):
			value[current] = value[current][0]
		else:
			value[current] = value[tree.getParent(species_tree,current)]
	return value[tree.getAncestor(species_tree)]
#enddef

gap_alignment_file=open(sys.argv[1],"r").readlines()
gap_lines={}
for l in gap_alignment_file:
	if l[0]==">":
		species=l.split('.')[0].split('>')[1]
	else:
		gap_lines[species] = l.split()[0]
		nbcolumns=len(gap_lines[species]) 

species_tree=tree.readTree(open(sys.argv[2],"r").readline())
sequence=""
for c in range(nbcolumns):
	character={}
	for s in gap_lines.keys():
		character[s]=gap_lines[s][c]
	value=Fitch(species_tree,character)
	if value in ["A","C","G","T"]:
		sequence = sequence + value
output_file=open(sys.argv[3],"w")
output_file.write(sequence)
