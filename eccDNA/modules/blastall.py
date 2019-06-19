#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""blast part of CIDER pashing pipeline"""

import os
import subprocess
import pickle
from Bio.Blast.Applications import NcbiblastnCommandline

#launch makeblastdb on a fasta file
def make_blastdb(subject_file, database):
	subject_short = subject_file.replace(".fasta", "")
	blastdb_cmd = 'makeblastdb -in {0}.fasta -dbtype nucl -out {1} -parse_seqids'.format(subject_short, database)

	DB_process = subprocess.Popen(blastdb_cmd,
	                              shell=True)
	DB_process.wait()
	
#launch blastn with Biopython
def blastn(database, query_file, blast_output, thread): 
	blastn_cmd = NcbiblastnCommandline(query=query_file, db=database, evalue=0.001, outfmt=6, out=blast_output, num_threads=thread)
	stdout, stderr = blastn_cmd() #no output to terminal

#launch blast on local version on nt db
def blastn_ncbi(database, query_file, blast_output, thread, mode):
	if (mode == 'local'):
	
		blastn_cmd = 'blastn -db {0} -query {1} -out {2} -outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" -evalue=1e-3 -num_alignments=10 -num_threads={3}'.format(database, query_file, blast_output, thread)
		
		DB_process = subprocess.Popen(blastn_cmd,
									  shell=True)
		DB_process.wait()
		
	elif(mode == 'remote'):
		
		blastn_cmd = 'blastn -db nt -query {0} -out {1} -outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" -evalue=1e-3 -num_alignments=10 -num_threads={2} -remote'.format(query_file, blast_output, thread)
		
		DB_process = subprocess.Popen(blastn_cmd,
									  shell=True)
		DB_process.wait()

#Get NCBI blast taxids output and use BIOmustCORE to get linage
def get_lineage(blast_output, taxdir):
	#parse blast results for taxids
	my_out = open('output/list.taxid', "w")
	
	infile = open(blast_output)
	
	seen_taxid = {}
	for line in infile:
		blast_record = line.replace("\n", "")
		split_liste = blast_record.split("\t")
		taxid = split_liste[12]
		
		if (taxid not in seen_taxid):
			print('pass')
			my_out.write(taxid + "\n")
			seen_taxid[taxid] = 1

			
	#run fetch tax
	fetch_cmd = 'fetch-tax.pl output/list.taxid  --taxdir={0} --item-type=taxid'.format(taxdir)

	DB_process = subprocess.Popen(fetch_cmd,
	                              shell=True)
	DB_process.wait()
	
	#parse fetch-tax output for dict creation
	lineage_of = {}
	
	infetch = 'output/list.tax'
	for line in infetch:
		tax_record = line.replace("\n", "")
		split_liste = tax_record.split("\t")
		taxid = split_liste[0]
		lineage = split_liste[3]
		lineage_of[taxid] = lineage
	
	#print list of ids
	with open('output/lineage.pickle', 'wb') as pickleout:
		pickle.dump(lineage_of, pickleout, protocol=pickle.HIGHEST_PROTOCOL)	
	
	
	
	
	   
#Test part of the modules, call as main
if __name__ == "__main__": 
	#Test of makeblastdb def
	#subject_file = '/home/luc/Documents/Python-course/CIDER-cassava/debug/dbfile.fasta'
	#database = '/home/luc/Documents/Python-course/CIDER-cassava/debug/test-debug/dbfile'
	#make_blastdb(subject_short, database)
	#Test of local blastn
	#query_file = '/home/luc/Documents/Python-course/CIDER-cassava/debug/queryfile.fasta'
	#query_short = query_file.replace(".fasta", "")
	#blast_output = '/home/luc/Documents/Python-course/CIDER-cassava/debug/test-debug/queryfile.blastn'
	#thread = 1
	#blastn(database, query_file, blast_output, thread)
	get_lineage('/media/vol2/scratch/lcornet/cassava/CIDER/script/output/replicatefile1.fasta_on_NCBI.blast6TAXIDS', '/media/vol2/scratch/lcornet/synteny/genomes/fasta/64-cyanos/I-1.5/Orthogroups/test/taxdir/')
	
