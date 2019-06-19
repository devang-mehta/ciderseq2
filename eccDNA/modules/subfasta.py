#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""sub fasta"""

import pickle

#filtering fasta based on positive ids
def fasta2_dico(fasta_input):
	#parse blast of input file to get sequence length : fasta reader
	#dico declare
	length_seq = {}
	seq_seq = {}
	
	infile = open(fasta_input)
	defline = ""
	sequence = ""
	count = 0
	for line in infile:
		fasta_record = line.replace("\n", "")
		if '>' in fasta_record:
			if (count > 0) :
				length = len(sequence)
				defline = defline.replace(">", "")
				length_seq[defline] = length
				seq_seq[defline] = sequence
				#next defline
				defline = fasta_record
				sequence = ""
			
			if (count == 0):
				count += 1
				defline = fasta_record

		else :
			#concat sequence
			sequence = sequence + str(fasta_record)
			
	length = len(sequence)
	defline = defline.replace(">", "")
	length_seq[defline] = length
	seq_seq[defline] = sequence
	
	#write with pickles
	#load dicts in tuples for writing with pickles
	my_dicts = [length_seq, seq_seq]
	pick_file = 'output/fasta2dico.pickle';
	with open(pick_file, 'wb') as pickleout:
		pickle.dump(length_seq, pickleout, protocol=pickle.HIGHEST_PROTOCOL)
		pickle.dump(seq_seq, pickleout, protocol=pickle.HIGHEST_PROTOCOL)
		
def reduce_fasta (ids_file, fasta_input, loop, name):
	
	#load list if replicate id of sequence
	deconcatfilters = []
	with open(ids_file, 'rb') as pickelin: #reading binary
		deconcatfilters = pickle.load(pickelin)
		
	#read fasta file
	infile = open(fasta_input)
	
	#outfile = 'output/file' + str(loop) + '.fasta'
	my_out = open('output/' + str(name) + str(loop) + '.fasta', "w")
	tag = 0
	for line in infile:
		fasta_record = line.replace("\n", "")
		if '>' in fasta_record:
			defline = fasta_record.replace(">", "")
			if defline in deconcatfilters:
				tag = 1
			else: 
				tag = 0
				
			#print if found in dico
			if (tag == 1):
				my_out.write('>' + defline + "\n")
		else:
			if (tag == 1):
				my_out.write(fasta_record + "\n")
				

		
def virus_fasta (ids_file, fasta_input, basename):
	
	print('ENTER')
	
	#load list if replicate id of sequence
	host_of = {}
	virus_of = {}
	with open(ids_file, 'rb') as pickelin: #reading binary
		host_of = pickle.load(pickelin)
		virus_of = pickle.load(pickelin)
		
	#read fasta file
	infile = open(fasta_input)
	
	#outfile = 'output/file' + str(loop) + '.fasta'
	my_virus = open('output/' + str(basename) + '-virus.fasta', "w")
	my_host = open('output/' + str(basename) + '-nonvirus.fasta', "w")
	
	hosttag = 0
	virustag = 0
	
	print(host_of)
	
	for line in infile:
		fasta_record = line.replace("\n", "")
		if '>' in fasta_record:
			#prepare next iteration
			hosttag = 0
			virustag = 0
			
			#defline
			split_liste = fasta_record.split(" ")
			defline = split_liste[0]
			defline = defline.replace(">", "")
			
			#checi if next sequence is virus or host
			if defline in host_of:
				hosttag = 1
			elif defline in virus_of:
				virustag = 0
				
			#print if found in dico
			if (hosttag == 1):
				my_host.write('>' + defline + "\n")
			elif (virustag == 1):
				my_virus.write('>' + defline + "\n")
				
		else:
			if (hosttag == 1):
				my_host.write(fasta_record + "\n")
			elif (virustag == 1):
				my_virus.write(fasta_record + "\n")
				
	
		
	
if __name__ == "__main__": 
	#fasta2_dico parser 
	#fasta_input = '/home/luc/Documents/Python-course/CIDER-cassava/debug/test-debug/deconcat-rep1.fa'
	#fasta2_dico(fasta_input)
	#reduce_fasta part
	virus_fasta('output/origin.pickle', 'input/lima_output.lbc14--lbc14-ccs-099-cider.fasta', 'lima_output.lbc14--lbc14-ccs-099-cider')
	
