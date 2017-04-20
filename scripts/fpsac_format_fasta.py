# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Formatting FASTA sequences to have the sequence on one line and
# creating a separate file recording the length of each sequences

# argument 1: Input  FASTA file
# argument 2: Output formatted FASTA file
# argument 3: Output file recording the length of each sequence

import sys,string

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

contigs_file=open(sys.argv[1],"r").readlines()
output_file1=open(sys.argv[2],"w")
output_file2=open(sys.argv[3],"w")

contig_count=0
for line in contigs_file:
    if line[0]=='>':
        if contig_count>0:
            output_file1.write(">"+contig_name)
            output_file2.write(">"+contig_name+" length "+str(len(contig_seq))+'\n')
            output_file1.write("\n"+contig_seq+"\n")
        contig_count+=1
        contig_name=line.split('>')[1].rstrip('\n').rstrip(' ')
        contig_seq=""
    else:
        contig_seq=contig_seq+line.strip('\n')
output_file1.write(">"+contig_name)
output_file2.write(">"+contig_name+" length "+str(len(contig_seq))+'\n')
output_file1.write("\n"+contig_seq+"\n")
                                                                                    
