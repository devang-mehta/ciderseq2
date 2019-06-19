#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""replicate finder"""

import pickle

#Parsing blast fmt6 result
def replicate_finder(blast6):
	
	#load dicts of subfasta with pickles
	length_seq = {}
	seq_seq = {}
	with open('output/fasta2dico.pickle', 'rb') as pickelin: #reading binary	
		length_seq = pickle.load(pickelin)
		seq_seq = pickle.load(pickelin)
	
	#parse results of blast to find replicate
	result_handle = open(blast6)
	replicate_id = {}	
	for line in result_handle:
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
		
		query_length = length_seq[query]
		
		if (float(identity) >= 99):# and (int(length) >= int(query_length)) :
			#replicate found
			replicate_id[query] = 1
			
	#print list of ids
	with open('output/replicatefinder.pickle', 'wb') as pickleout:
		pickle.dump(replicate_id, pickleout, protocol=pickle.HIGHEST_PROTOCOL)
				
		
if __name__ == "__main__": 
	pickel_file = '/home/luc/Documents/Python-course/CIDER-cassava/debug/test-debug/subfasta.pickle'
	blast6 = '/home/luc/Documents/Python-course/CIDER-cassava/debug/queryfile.blastn'
	pick_idsfile = '/home/luc/Documents/Python-course/CIDER-cassava/debug/test-debug/ids.pickle'
	replicate_finder(pickel_file, blast6, pick_idsfile )
