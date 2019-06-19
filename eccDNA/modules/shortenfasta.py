#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""sub fasta"""

import pickle

def shorten_fasta(fasta_input, loop):
	
	#load coordinate
	listchains = []
	with open('output/reorder.pickle', 'rb') as pickelin: #reading binary
		listchains = pickle.load(pickelin)
		
	#print(listchains)
	
	#load sequence from fasta
	seq_seq = {}
	sequence = str()
	
	infile = open(fasta_input)
	defline = "NA"
	for line in infile:
		fasta_record = line.replace("\n", "")
		if '>' in fasta_record:
			defline = defline.replace(">", "")
			seq_seq[defline] = sequence
			#next defline
			defline = fasta_record
			sequence = ""
		else :
			#concat sequence
			sequence = sequence + str(fasta_record)
	#load last line too
	defline = defline.replace(">", "")
	seq_seq[defline] = sequence
			
	
	#print(seq_seq)
	#reorder and fasta print
	fileout = 'output/reorder' + str(loop) + '.fasta'
	my_out = open(fileout, "w")	
	
	for coo in listchains:
		split_liste = coo.split("\t")
		length = len(split_liste)
		
		query = str(split_liste[0])
		
		part1 = split_liste[1]
		chunks_1 =  part1.split("-")
		min1 = int(chunks_1[0])-1 #-1 count from 0 in python and not in blast report
		max1 = int(chunks_1[1])-1
		
		part2 = split_liste[2]
		chunks_2 =  part2.split("-")
		min2 = int(chunks_2[0])-1
		max2 = int(chunks_2[1])-1
			
		
		seq = seq_seq[query]
		p1 = seq[min1:max1]
		p2 = seq[min2:max2]	
		final_seq = str(p1) + str(p2)
		
		#case when third hit
		if (length >= 4):
			part3 = split_liste[3]
			chunks_3 =  part3.split("-")
			min3 = int(chunks_3[0])-1
			max3 = int(chunks_3[1])-1
			p3 = seq[min3:max3]
			final_seq = str(p1) + str(p2) + str(p3)
		
		my_out.write('>' + query + "\n")
		#fasta: 60 charatcter length
		i = len(final_seq)
		start = 0
		end = 60
		while start < i:
			if end <= i:
				substring = final_seq[start:end]
				my_out.write(substring + "\n")
				start = start + 60
				end = end + 60
			elif end > i:
				substring = final_seq[start:i]
				my_out.write(substring + "\n")
				start = start + 60
				end = end + 60

def assess_fasta(fasta_input, blast6, loop):	
		
	#load sequence from fasta
	seq_seq = {}
	sequence = str()
	
	infile = open(fasta_input)
	defline = "NA"
	for line in infile:
		fasta_record = line.replace("\n", "")
		if '>' in fasta_record:
			defline = defline.replace(">", "")
			seq_seq[defline] = sequence
			#next defline
			defline = fasta_record
			sequence = ""
		else :
			#concat sequence
			sequence = sequence + str(fasta_record)
	#load last line too
	defline = defline.replace(">", "")
	seq_seq[defline] = sequence	
	
	#load blast output, 
	result_handle = open(blast6)
	current_query = 'NA'
	
	assess_list = {}
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
		
		if (query != current_query):
			query_sequence = seq_seq[query]
			query_length = len(query_sequence)
			percentage = float(length)/float(query_length) * 100
			if (percentage >= 95):
				assess_list[query] = 1
			current_query = query
	
	#print assess fasta file
	fileout = 'output/assess' + str(loop) + '.fasta'
	my_out = open(fileout, "w")	
	
	for query in assess_list:
		my_out.write('>' + query + "\n")
		#fasta: 60 charatcter length
		final_seq = seq_seq[query]
		i = len(final_seq)
		start = 0
		end = 60
		while start < i:
			if end <= i:
				substring = final_seq[start:end]
				my_out.write(substring + "\n")
				start = start + 60
				end = end + 60
			elif end > i:
				substring = final_seq[start:i]
				my_out.write(substring + "\n")
				start = start + 60
				end = end + 60		
		
		
	
				
				
if __name__ == "__main__": 
	#fasta2_dico parser 
	fasta_input = '/media/vol2/scratch/lcornet/cassava/CIDER/script/output/rfile2.fasta'
	loop = '2'
	shorten_fasta(fasta_input, loop)
	#assess_fasta('/home/luc/Documents/Python-course/CIDER-cassava/output/reorder1.fasta', '/home/luc/Documents/Python-course/CIDER-cassava/output/reorder1.fasta_on_genomedb.blastn6', '1')

	
	
