from Bio import SeqIO
from Bio.Seq import Seq


fasta_input = open('/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/stx_subtyping/stx2a_whole_genes.fa', 'r')

seq_dict = {}

for each in SeqIO.parse(fasta_input, 'fasta'):
	name = each.name
	sequence = str(each.seq)
	if sequence in seq_dict:
		seq_dict[sequence].append(name)
	else:
		seq_dict[sequence] = []
		seq_dict[sequence].append(name)

#print len(seq_dict)

#print seq_dict

for each in seq_dict:
	
	
	print
	print len(seq_dict[each]),
	print each,
	for every in seq_dict[each]:
		print every, '\t',
	
	#print seq_dict[each]
	#print len(seq_dict[each])
	#print '>', each

print
