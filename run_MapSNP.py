#################################### imports ####################################

import os, sys
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

def make_bwa_dir(output_root):
	if not os.path.isdir('%s/tmp/bwa' % output_root):
		system ('mkdir %s/tmp/bwa' % output_root)
	bwa_output_dir = '%s/tmp/bwa' % output_root
	return bwa_output_dir

def map_to_stx(ref, fastq_read1, fastq_read2, bwa_output_dir, s_name):
	print '######  Running bwa...' + str(datetime.time(datetime.now())).split('.')[0]
	# check that is bwa isn't found, it returns a sensible error, when it is being run by shell script

	#system ('/usr/local/bwa/bwa aln %s %s > %s/%s.1.sai' % (ref, fastq_read1, bwa_output_dir, s_name))
	#system ('/usr/local/bwa/bwa aln %s %s > %s/%s.2.sai' % (ref, fastq_read2, bwa_output_dir, s_name))
	#system ('/usr/local/bwa/bwa sampe %s %s/%s.1.sai %s/%s.2.sai %s %s > %s/%s.sam' % (ref, bwa_output_dir, s_name, bwa_output_dir, s_name, fastq_read1, fastq_read2, bwa_output_dir, s_name))

	system ('/usr/local/bwa-0.7.4/bwa mem -v 1 %s %s %s > %s/%s.sam' % (ref, fastq_read1, fastq_read2, bwa_output_dir, s_name))

def sam_to_bam(bwa_output_dir, s_name):
	print '######  Running sam_to_bam...' + str(datetime.time(datetime.now())).split('.')[0]
	system('samtools view -h -q 30 -bS -o %s/%s.unique.bam %s/%s.sam' % (bwa_output_dir, s_name, bwa_output_dir, s_name))
	system('samtools sort %s/%s.unique.bam %s/%s.unique.sorted' % (bwa_output_dir, s_name, bwa_output_dir, s_name))
	system('samtools index %s/%s.unique.sorted.bam' % (bwa_output_dir, s_name))
	if os.path.exists('%s/%s.unique.sorted.bam' % (bwa_output_dir, s_name)):
		system('rm -rf %s/%s.sam' % (bwa_output_dir, s_name))
		system('rm -rf %s/%s.unique.bam' % (bwa_output_dir, s_name))


#################################### main ####################################

fastq_read1, fastq_read2, output_root, ref = get_opts()
s_name = get_s_name(fastq_read1)
bwa_output_dir = make_bwa_dir(output_root)

#map_to_stx(ref, fastq_read1, fastq_read2, bwa_output_dir, s_name)
#sam_to_bam(bwa_output_dir, s_name)









