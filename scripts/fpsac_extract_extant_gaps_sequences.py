# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Generating a FASTA file for an extant gap

# argument 1 = Input  extant genomes FASTA file
# argument 2 = Input  gaps coordinates and lengths
# argument 3 = Output gaps FASTA file prefix

import sys
from operator import itemgetter

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

base_complement = {"A":"T","T":"A","C":"G","G":"C","R":"","M":""}
def complement(sequence):
    result = []
    for n in range(len(sequence)):
        result.append(base_complement[sequence[-(n+1)]])
    return ''.join(result)

extant_sequences_file=open(sys.argv[1],"r").readlines()
sequences={}
for l in extant_sequences_file:
    if l[0]==">":
        chr=l.split('>')[1].rstrip('\n')
    else:
        sequences[chr]=l.rstrip('\n')

gaps_file=open(sys.argv[2],"r").readlines()
for l in gaps_file:
    if l[0]==">":
        gap=l.split(' ')[0].split('>')[1]
        interval=l.split(' ')[len(l.split(' '))-1].rstrip('\n')
        if interval!="0-0" and interval!="-1--1":
            interval_min=int(interval.split('-')[0])
            interval_max=int(interval.split('-')[1])
            output_file=open(sys.argv[3]+gap+".fasta","w")
            generate_fasta_file=True
        else:
            generate_fasta_file=False
    elif len(l)>1 and generate_fasta_file:
        chr=l.split(':')[0]
        start=int(l.split(':')[1].split('-')[0])-1
        end=int(l.split(':')[1].split('-')[1].split(' ')[0])-1
        lg=end-start+1
        sign=l.split(' ')[1].rstrip('\n')
        if start<=end:
            output_file.write(">"+l.rstrip('\n')+" "+gap+"\n")
            if sign=="+":
                output_file.write(sequences[chr][start:end+1]+"\n")
            else:
                output_file.write(complement(sequences[chr][start:end+1])+"\n")
