
#################################### imports ####################################

from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os, sys, operator
from os import system
from datetime import datetime

#################################### functions ####################################

def usage():
	print '\nthis script blasts assembled contigs against a reference and parses output\n'
	print 'Usage: <path/to/fastq_read1> <path/to/fastq_read2> <path/to/contigs.fa> <output_root>'
	sys.exit()

def get_opts():
	if len(sys.argv) != 4:
		usage()
	else:
		return sys.argv[1], sys.argv[2], sys.argv[3]

def get_s_name(contigs):
	s_name = contigs.split('/')[-1].split('.')[0]
	return s_name

def make_output_dirs(output_root):
	if not os.path.isdir('%s/tmp/blast' % output_root):
		system ('mkdir %s/tmp/blast' % output_root)
	if not os.path.isdir('%s/results' % output_root):
		system('mkdir %s/results' % output_root)

	bwa_output_dir = '%s/tmp/blast' % output_root
	results_dir =  '%s/results' % output_root
	return bwa_output_dir, results_dir

def blast(query, database, out):
	print '######  BLASTing...' + str(datetime.time(datetime.now())).split('.')[0]
	blastx_cline = NcbiblastxCommandline(cmd = 'blastn', query = query , db = database, evalue = 1e-20, outfmt = 5, out = out)
	stdout, stderr = blastx_cline()

def parse_blast_output(infile, results_dir):
	print '######  Parsing BLAST results...' + str(datetime.time(datetime.now())).split('.')[0]
	blast_xml = open(infile, 'r')
	blast_records = NCBIXML.parse(blast_xml)
	s_name = (infile.split('/')[-1]).split('_')[0]
	outhandle = open('%s/%s.AssBLAST.txt' % (results_dir, s_name), 'w')
	i = 0
	results_dict = {}
	for record in blast_records:
		
		if record.alignments:
			#for alignment in record.alignments:
			alignment = record.alignments[0]
			hsp = alignment.hsps[0]
			#for hsp in alignment.hsps:
				#if i < 1:
				#	i += 1
			db_match = str(alignment.title.split('; ')[-1])
			db_match = db_match.replace('_', '|')
			db_match = db_match.replace(' ', '|').split('|')[5]

			if db_match in results_dict:
				results_dict[db_match] += hsp.positives
			else:
				results_dict[db_match] = hsp.positives
			
			#print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (s_name, record.query, len(hsp.match), hsp.positives, db_match, hsp.expect, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end)
	
	outhandle.write('%s\t' % s_name)
	#print results_dict
	sorted_results = sorted(results_dict.iteritems(), key = operator.itemgetter(1))
	#print sorted_results
	sorted_results.reverse()
	#print sorted_results

	for each in sorted_results:
		print each[0], each[1]
		outhandle.write('%s:%s\t' % (each[0], each[1])) 
	outhandle.write('\n')





#################################### main ####################################

contigs, output_root, ref = get_opts()
s_name = get_s_name(contigs)
blast_output_dir, results_dir = make_output_dirs(output_root)
blast_xml = '%s/%s.xml' % (blast_output_dir, s_name)

blast(contigs, ref, blast_xml)
parse_blast_output(blast_xml, results_dir)

















