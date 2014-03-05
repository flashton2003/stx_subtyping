from __future__ import division

#################################### imports ####################################

import os, sys, pysam
from os import system	
from datetime import datetime


#################################### functions ####################################


def usage():
	print '\nthis script maps specific reads against reference and parses output\n'
	print 'Usage: <path/to/fastq_read1> <path/to/fastq_read2> <path/to/contigs.fa> <output_root>'
	sys.exit()

def get_opts():
	if len(sys.argv) != 5:
		usage()
	else:
		return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def get_s_name(fastq_read1):
	s_name = fastq_read1.split('/')[-1].split('.')[0]
	return s_name

def make_output_dirs(output_root):
	if not os.path.isdir('%s/tmp/bwa' % output_root):
		system ('mkdir %s/tmp/bwa' % output_root)
	if not os.path.isdir('%s/results' % output_root):
		system('mkdir %s/results' % output_root)

	bwa_output_dir = '%s/tmp/bwa' % output_root
	results_dir =  '%s/results' % output_root
	return bwa_output_dir, results_dir

def map_to_stx(ref, fastq_read1, fastq_read2, bwa_output_dir, s_name):
	print '######  Running bwa...' + str(datetime.time(datetime.now())).split('.')[0]
	# check that is bwa isn't found, it returns a sensible error, when it is being run by shell script
	#system ('/usr/local/bwa/bwa aln %s %s > %s/%s.1.sai' % (ref, fastq_read1, bwa_output_dir, s_name))
	#system ('/usr/local/bwa/bwa aln %s %s > %s/%s.2.sai' % (ref, fastq_read2, bwa_output_dir, s_name))
	#system ('/usr/local/bwa/bwa sampe %s %s/%s.1.sai %s/%s.2.sai %s %s > %s/%s.sam' % (ref, bwa_output_dir, s_name, bwa_output_dir, s_name, fastq_read1, fastq_read2, bwa_output_dir, s_name))

	system ('bwa mem %s %s %s > %s/%s.sam' % (ref, fastq_read1, fastq_read2, bwa_output_dir, s_name))

def sam_to_bam(bwa_output_dir, s_name):
	print '######  Running sam_to_bam...' + str(datetime.time(datetime.now())).split('.')[0]
	system('samtools view -h -F 4 -bS -o %s/%s.unique.bam %s/%s.sam' % (bwa_output_dir, s_name, bwa_output_dir, s_name))
	system('samtools sort %s/%s.unique.bam %s/%s.unique.sorted' % (bwa_output_dir, s_name, bwa_output_dir, s_name))
	system('samtools index %s/%s.unique.sorted.bam' % (bwa_output_dir, s_name))
	if os.path.exists('%s/%s.unique.sorted.bam' % (bwa_output_dir, s_name)):
		system('rm -rf %s/%s.sam' % (bwa_output_dir, s_name))
		system('rm -rf %s/%s.unique.bam' % (bwa_output_dir, s_name))

def make_pile_up(bamfile, reference):
	bamfile = pysam.Samfile(bamfile, "rb")
	mapped_bases = {}
	for pileupcolumn in bamfile.pileup(reference):
		#print pileupcolumn
		for pileupread in pileupcolumn.pileups:
			#print pileupread
			#print '\tbase in base %s = %s' % (pileupcolumn.pos, pileupread.alignment.seq[pileupread.qpos])
			pos = pileupcolumn.pos
			base = pileupread.alignment.seq[pileupread.qpos].split()
			if pos in mapped_bases:
				mapped_bases[pos].append(base[0])
			else:
				mapped_bases[pos] = base
	#for each in mapped_bases:
	#	print each, mapped_bases[each]
	return mapped_bases

def parse_pileup(mapped_bases, reference):
	ref_to_uniq_pos_dict = {'stx2a_gi|14892|emb|X07865.1|':{1174:'G', 1179:'C', 1191:'G'}, 'stx2c_gi|15718404|dbj|AB071845.1|':{1054:'A', 1173:'A', 1191:'A'}, 'stx2d_gi|30313370|gb|AY095209.1|':{1037:'C', 1173:'A', 1178:'T'}, 'stx2b_gi|49089|emb|X65949.1|':{1060:'C', 1083:'A', 1158:'C'}, 'stx2e_gi|8346567|emb|AJ249351.2|':{966:'C', 1108:'G', 1192:'G'}, 'stx2f_gi|254939478|dbj|AB472687.1|':{1015:'C', 964:'G', 897:'G'}, 'stx2g_gi|30909079|gb|AY286000.1|':{914:'C', 938:'A', 1079:'C'}, 'stx1a_gi|147832|gb|L04539.1|':{468:'A', 528:'A', 712:'G'}, 'stx1c_gi|535088|emb|Z36901.1|':{741:'T', 905:'G', 922:'T'}, 'stx1d_gi|28192582|gb|AY170851.1|':{478:'C', 501:'G', 508:'T'}}
	uniq_pos = ref_to_uniq_pos_dict[reference]
	i = 0
	for each in mapped_bases:
		#print each, mapped_bases[each]
		if each in uniq_pos:
			
			depth = len(mapped_bases[each])
			uniq_mapping = mapped_bases[each].count(uniq_pos[each])
			if depth > 10:
				#print uniq_mapping / depth
				if uniq_mapping / depth >= 0.9:
					i += 1
					#print each, mapped_bases[each]

	if i == 3:
		return reference
		#print '%s is present in sample' % reference

			#print each, mapped_bases[each]


def run_pileup_analysis(bwa_output_dir, s_name, results_dir):
	bamfile = '%s/%s.unique.sorted.bam' % (bwa_output_dir, s_name)
	references = ['stx2a_gi|14892|emb|X07865.1|' ,'stx2c_gi|15718404|dbj|AB071845.1|', 'stx2d_gi|30313370|gb|AY095209.1|', 'stx2b_gi|49089|emb|X65949.1|', 'stx2e_gi|8346567|emb|AJ249351.2|', 'stx2f_gi|254939478|dbj|AB472687.1|', 'stx2g_gi|30909079|gb|AY286000.1|', 'stx1a_gi|147832|gb|L04539.1|', 'stx1c_gi|535088|emb|Z36901.1|', 'stx1d_gi|28192582|gb|AY170851.1|']
	#references = ['stx2a_gi|14892|emb|X07865.1|']
	outhandle = open('%s/%s.MapSNP.txt' % (results_dir, s_name), 'w')
	outhandle.write('%s\t' % (s_name))
	for reference in references:
		#print reference
		mapped_bases = make_pile_up(bamfile, reference)
		matched_stx = parse_pileup(mapped_bases, reference)
		if matched_stx is not None:
			matched_stx = matched_stx.split('_')[0]
			outhandle.write('%s\t' % matched_stx)
	outhandle.write('\n')
	outhandle.close()



#################################### main ####################################

fastq_read1, fastq_read2, output_root, ref = get_opts()
s_name = get_s_name(fastq_read1)
bwa_output_dir, results_dir = make_output_dirs(output_root)

map_to_stx(ref, fastq_read1, fastq_read2, bwa_output_dir, s_name)
sam_to_bam(bwa_output_dir, s_name)

run_pileup_analysis(bwa_output_dir, s_name, results_dir)







