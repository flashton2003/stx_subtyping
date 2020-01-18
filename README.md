### stx_subtyping: a tool for subtyping shiga toxin genes from fastqs or fastas.

**Installation**

git clone https://github.com/flashton2003/stx_subtyping

**Requirements**

Tested with Biopython v1.6.3
Test with pysam v0.6
Tested with bwa v0.7.5
Need to have bwa installed and available in your PATH (i.e. when you type `bwa` from any location, it runs bwa)

Tested with python v2.7.6
Tested with Mac OS X 10.9


```
conda create -n stx_subtyping python=2.7 bwa biopython pysam samtools=1.1
```

**Usage**

python /path/to/script/stx_subtyping.py -h
python /path/to/script/stx_subtyping.py dual_stx_subtyping -h

python /path/to/script/stx_subtyping.py dual_stx_subtyping /path/to/fastq1.fastq /path/to/fastq2.fastq /path/to/contigs.fa /path/to/whereyouwanttheoutput

There are 3 running options. The best approach is dual_stx_subtyping, but if only fastqs or fastas are available, the other options can be used.

dual\_stx\_subtyping - Takes a pair of fastqs and a set of contigs for a strain, outputs info into output_dir

assblast_only - Use this when you only have contigs

mapsnp_only - Use this when you only have fastqs, no assemblies

If in doubt, try -h
