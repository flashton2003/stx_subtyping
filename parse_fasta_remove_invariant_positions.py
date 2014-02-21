from Bio import SeqIO
from Bio.Seq import Seq
'''

inhandle = '/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/stx_subtyping/stx2c_whole_genes.fa'
fasta_input = open(inhandle, 'r')
fasta_output = open('/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/stx_subtyping/stx2c_variant_positions.fa', 'w')

seq_dict = {}

for each in SeqIO.parse(fasta_input, 'fasta'):
	i = 0
	for every in each.seq:
		i += 1
		if i in seq_dict:
			seq_dict[i].append(every)
		else:
			seq_dict[i] = []
			seq_dict[i].append(every)

#print len(seq_dict)
#print seq_dict
keep_list = []

for each in seq_dict:
	if len(set(seq_dict[each])) > 1:
		keep_list.append(each)

#print keep_list
print len(keep_list)

fasta_input.close()

fasta_input = open(inhandle, 'r')

for each in SeqIO.parse(fasta_input, 'fasta'):
	seq_to_keep = ''
	i = 0
	for every in each.seq:
		i += 1
		if i in keep_list:
			seq_to_keep += every
	fasta_output.write('>%s\n' % each.description)
	fasta_output.write('%s\n' % seq_to_keep)
	

fasta_output.close()
fasta_input.close()

'''

fasta_input = open('/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/stx_subtyping/stx2c_variant_positions.fa', 'r')

letter_to_number_dict = {'A':1, 'G': 2, 'C': 3, 'T': 4, 'Y':5, 'V':6, 'D':7, 'R':8, 'K':9}

for each in SeqIO.parse(fasta_input, 'fasta'):
	print
	print each.description, '\t',
	for every in each:
		print letter_to_number_dict[every],

fasta_input.close()
'''