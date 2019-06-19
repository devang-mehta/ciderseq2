#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""assessment of origin"""

import pickle

#Parsing blast fmt6 result
def origin_assess(blast6):
	
	#load dicts of subfasta with pickles
	length_seq = {}
	seq_seq = {}
	with open('output/fasta2dico.pickle', 'rb') as pickelin: #reading binary	
		length_seq = pickle.load(pickelin)
		seq_seq = pickle.load(pickelin)
		
	#define dico for results
	blast_out = {}
	#read blast 
	count = 0
	result_handle = open(blast6)
	assess_id = {}	
	for line in result_handle:
		count += 1
		blast_record = line.replace("\n", "")
		split_liste = blast_record.split("\t")
		query = split_liste[0]
		subject = split_liste[1]
		identity = split_liste[2]
		length = split_liste[3]
		mismatch = split_liste[4]
		gap = split_liste[5]
		qstart = split_liste[6]
		qend = split_liste[7]
		sstart = split_liste[8]
		send = split_liste[9]
		evalue = split_liste[10]
		bitscore = split_liste[11]	
		
		#compute percent length
		query_length = length_seq[query]
		percent_len = (float(length)/float(query_length)) * 100
		
		if (float(percent_len) >= 20) and (float(identity) >= 95) and (float(evalue) <= 1e-3):
			if not query in assess_id:
				assess_id[query] = 1
			else: 
				assess_id[query] += 1
				
	#delete query if too much hits
	temp_dic = {}
	for key in assess_id :
		num = assess_id[key]
		if (int(num) < 50):
			temp_dic[key] = 1	
	assess_id = temp_dic
	#print(assess_id)
			
	#print list of ids
	with open('output/originassess.pickle', 'wb') as pickleout:
		pickle.dump(assess_id, pickleout, protocol=pickle.HIGHEST_PROTOCOL)

#parse NCBI blast results along with linage file produce by fetch tax: grep on virus or not		
def ncbi_assess(blast_output):
	
	#load lineage database
	lineage_of = {}
	infetch = open('input/ncbi.tax')
	for line in infetch:
		tax_record = line.replace("\n", "")
		split_liste = tax_record.split("\t")
		taxid = split_liste[0]
		lineage = split_liste[3]
		lineage_of[taxid] = lineage

			
	#read blast output
	host_of = {}
	virus_of = {}
	ids = 'NA'
	
	infile = open(blast_output)
	for line in infile:
		blast_record = line.replace("\n", "")
		split_liste = blast_record.split("\t")
		query = split_liste[0]
		subject = split_liste[1]
		identity = split_liste[2]
		length = split_liste[3]
		mismatch = split_liste[4]
		gap = split_liste[5]
		qstart = split_liste[6]
		qend = split_liste[7]
		sstart = split_liste[8]
		send = split_liste[9]
		evalue = split_liste[10]
		bitscore = split_liste[11]
		taxid = split_liste[12]
		
		lineage = "NA"
		if taxid in lineage_of:
			lineage = lineage_of[taxid]

		if (query != ids): #only first hit
			ids = query
			
			if ("Virus" not in lineage):
				host_of[query] = 1
			elif ("Virus" in lineage):
				virus_of[query]= 1
	
	
	#print host_in in pickle
	with open('output/origin.pickle', 'wb') as pickleout:
		pickle.dump(host_of, pickleout, protocol=pickle.HIGHEST_PROTOCOL)	
		pickle.dump(virus_of, pickleout, protocol=pickle.HIGHEST_PROTOCOL)		
				
			
	
	
			
				
if __name__ == "__main__": 
	blast6 = '/home/luc/Documents/Python-course/CIDER-cassava/output/replicatefile1_on_genomedb.blastn6'
	ncbi_assess('/media/vol2/scratch/lcornet/cassava/CIDER/script/final/output/lima_output.lbc14--lbc14-ccs-099-cider_on_NCBI.blast6TAXIDS')
