#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""cdhit parser"""


def shorten_fasta(cdhit):
	
	#load sequence from fasta
	cluster = 0
	sequence = 0
	
	infile = open(cdhit)
	tag = 0
	
	file1 = 0
	file2 = 0
	
	FILE1 = 0
	FILE2 = 0
	
	for line in infile:
		cdhit_record = line.replace("\n", "")
		if 'Cluster' in cdhit_record:
			#reset
			tag = 0
			
			if (file1 > 0) and (file2 > 0):
				cluster += 1
				sequence = file1 + file2 + sequence
				FILE1 += file1
				FILE2 += file2
				
			file1 = 0
			file2 = 0
			#next defline
			tag += 1
		else:
			#count sequence origin
			if (tag == 1):
				if 'file1' in cdhit_record:
					file1 += 1
					
				elif 'file2' in cdhit_record:
					file2 += 1
					
	print(cluster)
	print(sequence)
	print(FILE1)
	print(FILE2)
				
			
if __name__ == "__main__":
	shorten_fasta('cdhit.fasta.clstr')

