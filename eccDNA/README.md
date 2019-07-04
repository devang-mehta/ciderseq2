# cider

Usage: CIDER-main.py [OPTIONS] LIST_FILE GENOME_FILE

  Script CIDER

Options:
  --blastn_mode TEXT          Blastn mode : local or remote, local recommended
  
  --path_2_ncbi_blastdb TEXT  Input local nt DB
  
  --blast_threads INTEGER     Number of cores for Blast
  
  --gap_window INTEGER        Number of gap between two hits
  
  --help                      Show this message and exit.
  
  
#Arguments of CIDER-main.py

#LIST_FILE : tab separate file with CIDER-deconcat output of replicates (same samples)

input/lima_output.lbc14--lbc14-ccs-099-cider.fasta  input/lima_output.lbc14--lbc14-ccs-099-cider.stat

input/lima_output.lbc22--lbc22-ccs-099-cider.fasta  input/lima_output.lbc22--lbc22-ccs-099-cider.stat

#GENOME_FILE : path to genome file, like 

PATH-2-genome.fasta

#Example

#Take input directory : must contain a taxonomy file - provided can be update with Bio-MustCore perl

ncbi.tax

#Create output directory - this is where the program write

$ mkdir output

$ ./CIDER-main.py input/WT.list /media/vol2/scratch/lcornet/cassava/CIDER/script/input/Athaliana_167_TAIR9.fasta

#Prerequiste

ncbi blast in PATH

#Local ncbi nt database - no recommended to use remote

update_blastdb.pl (furnished with NCBI-BLAST)

  
