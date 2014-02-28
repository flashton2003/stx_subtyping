from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os

######################################## inputs etc ########################################

# Single BLAST
'''
blast_input = '/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/velvet_contigs/H121420197_contigs.fa'
blast_db = '/Users/philip/vtec_proj/stx_analysis_16.01.14/refs/stx2a.fa'
blast_db_name = blast_db.split('/')[-1].split('.')[0]
xml_output = '/Users/philip/vtec_proj/stx_analysis_16.01.14/blast_output/H121420197_vs_stx2a.xml'
'''

# multiple BLASTs

contigs_directory = '/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/velvet_contigs'
blast_db = '/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/ref/stx_genes_with_flanking_regions.fa'
blast_db_name = blast_db.split('/')[-1].split('.')[0]
output_root = '/Volumes/NGS2_DataRAID/projects/tim/vtec_grant/stx_subtyping/allele_analysis'

to_do_list = ['H120300214', 'H103620472', 'H123680355', 'H124380400', 'H103120149', 'H124040469']

not_to_do_list = ['E64497', 'H092220118', 'H093100566', 'H093360372', 'H094080746', 'H094680415', 'H102440973', 'H102980326', 'H103460518', 'H103620472', 'H110320350', 'H113480397', 'H113580158', 'H121320380', 'H121320381', 'H121320842', 'H121380674', 'H121620536', 'H121700348', 'H121800662', 'H121800664', 'H121800665', 'H121860384', 'H122000289', 'H122260441', 'H122340673', 'H122360382', 'H122360386', 'H122360387', 'H122380383', 'H122420237', 'H122420245', 'H122460427', 'H122480279', 'H122540612', 'H122580782', 'H122600406', 'H122620619', 'H122660484', 'H122700503', 'H122700504', 'H122720636', 'H122760538', 'H122860602', 'H122880421', 'H122900485', 'H122920159', 'H122920160', 'H122920161', 'H122920162', 'H122920499', 'H122940311', 'H122960360', 'H122960363', 'H122980188', 'H122980190', 'H122980191', 'H122980192', 'H122980845', 'H123000125', 'H123000129', 'H123000130', 'H123020474', 'H123020477', 'H123020478', 'H123020480', 'H123040572', 'H123060304', 'H123060819', 'H123120508', 'H123120509', 'H123120511', 'H123120512', 'H123120513', 'H123140453', 'H123140834', 'H123180608', 'H123240465', 'H123280435', 'H123300001', 'H123380175', 'H123540539', 'H123560642', 'H123560644', 'H123560647', 'H123560650', 'H123580709', 'H123600594', 'H123660383', 'H123680355', 'H123680356', 'H123680361', 'H123740296', 'H123820321', 'H123840423', 'H123880472', 'H123880634', 'H123960529', 'H124180619', 'H124180620', 'H124180621', 'H124180625', 'H124240276', 'H124240277', 'H124240278', 'H124280635', 'H124280637', 'H124280639', 'H124280643', 'H124280644', 'H124280645', 'H124280647', 'H124280649', 'H124340142', 'H124340144', 'H124340145', 'H124340146', 'H124340152', 'H124340153', 'H124340154', 'H124340157', 'H124340160', 'H124340163', 'H124340168', 'H124340172', 'H124340173', 'H124340174', 'H124340180', 'H124340184', 'H124380366', 'H124380371', 'H124380378', 'H124380379', 'H124380382', 'H124380394', 'H124380398', 'H124380399', 'H124380400', 'H124380409', 'H124380411', 'H124380414', 'H124380417', 'H124380420', 'H124380422', 'H124380425', 'H124380426', 'H124380427', 'H124380434', 'H124480895', 'H124680285', 'H124680287', 'H125160541', 'H125280570', 'H125281263', 'H130300235', 'H130600599']
#to_do_list = ['H094080746']

######################################## functions  ########################################


def blast(query, database, out):
	blastx_cline = NcbiblastxCommandline(cmd = 'blastn', query = query , db = database, evalue = 1e-20, outfmt = 5, out = out)
	stdout, stderr = blastx_cline()

def parse_blast_output(infile):
	output_handle = open(infile, 'r')
	blast_records = NCBIXML.parse(output_handle)
	ins_cont_dict = {}
	s_name = (infile.split('/')[-1]).split('_')[0]
	for record in blast_records:
		i = 0
		if record.alignments:
			for alignment in record.alignments:
				for hsp in alignment.hsps:
					#if i < 1:
					#	i += 1
					db_match = str(alignment.title.split('; ')[-1])
					if db_match == 'gnl|BL_ORD_ID|1 Stx2c_gi|15718404|dbj|AB071845.1|':
						if len(hsp.match) > 1000:
						#print record.query
							#print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (s_name, record.query, len(hsp.match), hsp.positives, db_match, hsp.expect, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end)
							
							print '>%s_2c' % s_name
							print hsp.query
							'''
							correct_orientation = hsp.query[100:]
							#print correct_orientation
							if correct_orientation.startswith('ATG'):
								correct_orientation = correct_orientation[:-100]
								print correct_orientation
							else:
								reverse_orientation = hsp.query[100:]

								reverse_orientation = reverse_orientation[:-100]
								reverse_orientation = Seq(reverse_orientation, IUPAC.unambiguous_dna)
								#print reverse_orientation.reverse_complement()
								if reverse_orientation.reverse_complement().startswith('ATG'):
									
									print reverse_orientation.reverse_complement()
							'''

def run_blast(contigs_directory, blast_db, output_root):
	for each_dir in os.listdir(contigs_directory):
		s_name = each_dir.split('_')[0]
		if s_name in to_do_list:
			#print each_dir
			inhandle = '%s/%s' % (contigs_directory, each_dir)
			#print 's_name is', s_name[0]
			blast_results = '%s/%s_vs_%s.xml' % (output_root, s_name, blast_db_name)
			blast(inhandle, blast_db, blast_results)


def run_parse_output(output_root):
	for each in os.listdir(output_root):
		if each.split('_')[0] in to_do_list:
			if each.endswith('%s.xml' % blast_db_name):
				inhandle = '%s/%s' % (output_root, each)
				#print every
				parse_blast_output(inhandle)

######################################## main ########################################



#blast(blast_input, blast_db, xml_output)
#parse_blast_output(xml_output)

#run_blast(contigs_directory, blast_db, output_root)
run_parse_output(output_root)











