
#!/bin/bash
#$ -N ITSxpress
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -l h='symbiosis'
#$ -l mem_free=15G



for infile in *_R1_001.fastq
do
base="$(basename ${infile} _R1_001.fastq)"
itsxpress --fastq ${infile} \
--fastq2 ${base}_R2_001.fastq \
--log ${base}_logfile.txt \
--outfile ${base}_filtered.fastq \
--region ITS1 --taxa Fungi --threads 2
done