#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""ciderseq-phasing.py: Main script of the cider phasing algorithm."""
__author__ = "Luc Cornet, Syed Shan-e-Ali Zaidi"
__copyright__ = "Copyright 2019, University of Liège"
__version__ = "1.0.0"
__maintainer__ = "Luc Cornet"
__email__ = "luc.cornet@uliege.be"
__status__ = "Production"


import click
import os
import glob

#modules import
from modules.subfasta import *
from modules.blastall import *
from modules.replicatefinder import *
from modules.deconcatparser import *
from modules.originassess import *
from modules.order import *
from modules.shortenfasta import *

@click.command()
###ARGUMENT#####
#Input tab separate list of deconcat files input/WTT.list
@click.argument('list_file', type=click.Path(exists=True,readable=True))
#Genome host file path : input/manihot_esculenta_07Oct2017_aWswf_tme3.fasta
@click.argument('genome_file', type=click.Path(exists=True,readable=True))
###OPTIONS#####
@click.option('--ncbiblast', default='no', help='NCBI blast, yes or no - first separate sequences between virus and host')
@click.option('--blastn_mode', default='local', help='Blastn mode : local or remote, local recommended')
@click.option('--path_2_ncbi_blastdb', default='/media/vol1/databases/nt-20170809/nt', help='Input local nt DB')
@click.option('--blast_threads', default=40, help='Number of cores for Blast')
@click.option('--gap_window', default=150, help='Number of gap between two hits')


def main(list_file, genome_file, ncbiblast, path_2_ncbi_blastdb, blast_threads, gap_window, blastn_mode):
	"""Script CIDER"""
	##take list of replicates samples
	seen_files = {}
	infile = open(list_file)
	for line in infile:
		list_record = line.replace("\n", "")
		split_list = list_record.split("\t")
		fasta = split_list[0]
		stat = split_list[1]
		seen_files[fasta] = stat
		
	##Check Origin of sequence based on NCBI blast	
	if (ncbiblast == 'yes'): #not do by default
		for fasta in sorted(seen_files):
			#get base name
			basename = fasta.replace("input/", "")
			basename = basename.replace(".fasta", "")
			#blast on ncbi : no possibility of remote, only on local copy of ncbi nt db : database, query_file, blast_output, thread
			database = path_2_ncbi_blastdb
			blast_output = "output/" + basename + '_on_NCBI.blast6TAXIDS'
			blastn_ncbi(database, fasta, blast_output, blast_threads, blastn_mode)
			#read blast output and determine origin of sequence : Virus or not
			ncbi_assess(blast_output)
			#produce final fasta with reduce set of non virus sequence : dict of sequence, original fasta file, loop or id, basename in output/ dir
			#original fasta file reduce to host and virus
			virus_fasta('output/origin.pickle', fasta, basename)
		
	
	##reduce fasta based on deconcat outstat and order hits
	loop = 0
	for fasta in sorted(seen_files):
		loop += 1
		#get cider output
		stat = seen_files[fasta]
		#modify name to get non-virus sequence : olny do the modification if Virus disctionction done by NCBI blast
		if (ncbiblast == 'yes'):
			fasta = fasta.replace("input/", "output/")
			fasta = fasta.replace(".fasta", "-nonvirus.fasta")
		#produce list of > 0 round of deconcat
		deconcat_parser(stat)
		#filter fasta file (file1.fasta) based on deconcat list : : dict of sequence, original fasta file, loop or id, basename
		reduce_fasta('output/deconcat.pickle', fasta, loop, 'rfile')
		#reorganise sequence
		rfile_name = 'output/rfile' + str(loop) + '.fasta'
		rfile_blast = 'output/rfile' + str(loop) + '.fasta_on_genomedb.blastn6'
		#makeblastdb of genome : file, db
		make_blastdb(genome_file, 'output/genomedb')
		#blast verified origin fasta file on host genome : db, query file, blast6 output, threads
		blastn('output/genomedb', rfile_name, rfile_blast, blast_threads)
		#try to find et re-organise partial hits
		reorder_sequence(rfile_blast, gap_window)
		#produce fasta with new coordinate
		shorten_fasta(rfile_name, loop)
		#blast reorder fasta file on host genome : db, query file, blast6 output, threads
		reorder_name = 'output/reorder' + str(loop) + '.fasta'
		reorder_blast = 'output/reorder' + str(loop) + '.fasta_on_genomedb.blastn6'
		blastn('output/genomedb', reorder_name, reorder_blast, blast_threads)
		#check how many sequence can be detected into one hit after reordering
		assess_fasta(reorder_name, reorder_blast, loop)
		
	##find replicates among assess files
	assess_lists = glob.glob("output/assess*.fasta")
	#loop in list 2by2
	assess_length = len(assess_lists)
	iterid = 0
	while assess_length > 1: #at least two files to compare
		#iterid = 0 
		for assess_main in assess_lists:
			iterid += 1
			
			print(iterid)
			print(assess_lists)
		
			#collect length of sequence
			fasta2_dico(assess_main)
			#make blastdb of second file (sub) : file, db
			assess_sub = assess_lists[1]
			db = assess_sub.replace(".fasta", "")
			make_blastdb(assess_sub, db)
			#launch blast with file1 in query : db, query file, blast6 output, threads
			db_name = db.replace("output/", "")
			blast_out = assess_main + '_on_' + db_name + '.blastn6'
			blastn(db, assess_main, blast_out , 1)
			#turn replicate finder
			replicate_finder(blast_out)
			#reduce fasta from replicate finder list : dict of sequence, original fasta file, loop or id, basename in output/ dir
			reduce_fasta('output/replicatefinder.pickle', assess_main, iterid, 'replicatefile')
			#Prepare next loop
			#First : delete the two files investigate 2by2
			assess_lists.remove(assess_main)
			assess_lists.remove(assess_sub)
			#add the resulting file: replicatefile to assess list
			rep_file = 'output/replicatefile' + str(iterid) + '.fasta'
			assess_lists.append(rep_file)
			#break the for loop
			break
		#re-evaluate length of list of assess-replicate files
		assess_length = len(assess_lists)
	
	
if __name__ == '__main__':
    main()
