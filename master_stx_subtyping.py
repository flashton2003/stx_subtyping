from os import system
import os, glob, re, sys

###################################### variables and outdirs ######################################


#root_dir2 = root_dir.split('/')[:-1]
#root_dir2 = '/'.join(str(d) for d in root_dir2)
home = os.path.expanduser('~')

tool = 'stx_subtype'
scripts_dir = '%s/scripts/%s_scripts' % (home, tool)
logs_dir = '%s/logs' % (home)

stx_pipeline_scripts_dir = os.path.dirname(os.path.realpath(__file__))

if not os.path.exists(scripts_dir):
	os.makedirs(scripts_dir)

if not os.path.exists(logs_dir):
	os.makedirs(logs_dir)



###################################### functions ######################################

def clean_up_script_and_log_dirs(scripts_dir):
	system ('rm -rf %s/*' % scripts_dir)
	system('rm -rf %s/log.%s*' % (logs_dir, tool))

def get_opts():
	if len(sys.argv) != 5:
		usage()
	else:
		return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def usage():
	print "\nScript to run stx_subtyping and on a directory of paired fastqs and a directory of contigs\nWrites and qsubs a bunch of shell scripts to ~/scripts/velvet_scripts\nThe path to the contigs and the fastqs has to be the absolute path, the list can be relative.\n"
	print 'Usage: <list> <fastq_dir - can be .gz> <assembly_dir - flat dir of _contigs.fa files> <out dir, where do you want your assemblies>'
	sys.exit()

def make_output_dir(outdir):
	if not os.path.exists(outdir):
		os.makedirs(outdir)

def make_to_do_list(list_file):
	to_do_list = []
	with open(list_file, 'r') as fo:
		for sample in fo:
			sample = sample.strip()
			to_do_list.append(sample)
	return to_do_list



def write_stx_subtyping_scripts(list_file, scripts_dir, fastq_dir, assembly_dir, outdir):
	read1_re = re.compile(r'1.fastq')
	read2_re = re.compile(r'2.fastq')
	to_do_list = make_to_do_list(list_file)
	for each in os.listdir(fastq_dir):
		if each.split('.')[0] in to_do_list:
			R1 = ''
			R2 = ''
			sample_name = each.split('.')[0]
			#if sample_name in to_do_list:
			#print sample_name
			for each in glob.glob('%s/%s*' % (fastq_dir, sample_name)):
				m = read1_re.search(each)
				if m:
					R1 = each
					read_name_1 = R1.split('/')[-1].split('.')
					read_name_1 = '.'.join(read_name_1[0:2])
					for each in glob.glob('%s/%s*' % (fastq_dir, sample_name)):
						n = read2_re.search(each)
						if n:
							R2 = each
							contigs = '%s/%s_contigs.fa' % (assembly_dir, sample_name)
							try:
								os.path.isfile(contigs)
							except:
								print 'there is no corresponding contigs file \'%s\', naming convention is samplename_contigs.fa' % (contigs)				
							outhandle = open('%s/%s_stx_subtype.sh' % (scripts_dir, sample_name), 'w')
							command = '#! /bin/bash\n#$ -o %s/log.stx_subtype.stdout\n#$ -e %s/log.stx_subtype.stderr\n#$ -m e\n#$ -wd %s\n#$ -N stx_subtype_%s\n. /etc/profile.d/modules.sh\nmodule load bwa/0.7.5a\nmodule load blast+/2.2.27\nmodule load python/2.7.5\nmodule load biopython/python2.7/1.61\nmodule load samtools/0.1.19\nmodule load pysam/python2.7/0.7.5\n\npython %s/stx_subtyping.py %s %s %s %s\n' % (logs_dir, logs_dir, home, sample_name, stx_pipeline_scripts_dir, R1, R2, contigs, outdir)
							outhandle.write(command)
							outhandle.close()

def execute_stx_subtyping_scripts(scripts_dir):
	for each in os.listdir(scripts_dir):
		system ('qsub %s/%s' % (scripts_dir, each))


###################################### main ######################################

clean_up_script_and_log_dirs(scripts_dir)

list_file, fastq_dir, assembly_dir, outdir = get_opts()

make_output_dir(outdir)

write_stx_subtyping_scripts(list_file, scripts_dir, fastq_dir, assembly_dir, outdir)

execute_stx_subtyping_scripts(scripts_dir)
















