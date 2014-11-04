#################################### imports ####################################

import sys, os
from os import system
from datetime import datetime

### think about snakemaking this pipeline.

#################################### functions ####################################

def usage():
	print '\nPipeline to subtype stx genes from assembly fasta and fastq. This script will subtype a single strain, it is expected that for high throughput useage a command line loop or shell script will be used to call this script. Put about naming conventions of input files. Put about expected naming convention i.e. sample_name.<R>1.fastq\n'
	print 'Usage: <path/to/fastq_read1> <path/to/fastq_read2> <path/to/contigs.fa> <output_root - default to pwd>'
	sys.exit()

def get_opts():
	if len(sys.argv) != 5:
		usage()
	else:
		return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def make_output_dir(output_root):
	if not os.path.isdir(output_root):
		system ('mkdir %s' % output_root)
	if not os.path.isdir('%s/tmp' % output_root):
		system ('mkdir %s/tmp' % output_root)

def get_scripts_dir():
	return os.path.dirname(os.path.realpath(__file__))

def get_ref_fasta():
	scripts_dir = os.path.dirname(os.path.realpath(__file__))
	ref = '%s/refs/stx_genes_with_flanking_regions.fa' % scripts_dir
	return ref
	# insert test to check for correct blast, bwa indexes

def run_MapSNP(fastq_read1, fastq_read2, output_root):
	scripts_dir = get_scripts_dir()
	ref = get_ref_fasta()
	print '###  START:' + str(datetime.time(datetime.now())).split('.')[0]
	print '###  Running MapSNP workflow...' + str(datetime.time(datetime.now())).split('.')[0]
	system ('python %s/run_MapSNP.py %s %s %s %s' % (scripts_dir, fastq_read1, fastq_read2, output_root, ref))
	

def run_AssBLAST(contigs, output_root):
	scripts_dir = get_scripts_dir()
	ref = get_ref_fasta()
	print '### START:' + str(datetime.time(datetime.now())).split('.')[0]
	print '###  Running AssBLAST workflow...' + str(datetime.time(datetime.now())).split('.')[0]
	system ('python %s/run_AssBLAST.py %s %s %s' % (scripts_dir, contigs, output_root, ref))



#################################### main ####################################



fastq_read1, fastq_read2, contigs, output_root = get_opts()
make_output_dir(output_root)
#run_MapSNP(fastq_read1, fastq_read2, output_root)
run_AssBLAST(contigs, output_root)
