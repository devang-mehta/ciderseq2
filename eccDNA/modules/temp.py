#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""reorder CIDER hits"""

import pickle

#Parsing blast fmt6 result
def reorder_sequence(blast6):
	#read blast 
	count = 0
	result_handle = open(blast6)
	
	parse_of = {}
		
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
		
		query_id = query + '-' + str(count)
		
		if (float(identity) >= 99):
			parse_of[query_id] = dict()
			parse_of[query_id]['scaffold'] = subject
			parse_of[query_id]['coord_query'] = [qstart,qend]
			parse_of[query_id]['coord_subject'] = [sstart,send]
			
			
	#find sequence separate in two hits : double for loop
	listchains = []
	for main_query in sorted(parse_of):
		#declare query coordinate for Master hit
		split_liste = main_query.split("-")
		main = split_liste[0]
		main_scaffold = parse_of[main_query]['scaffold']
		main_coords = parse_of[main_query]['coord_query']
		min_main = main_coords[0]
		max_main = main_coords[1]
		#take sequence at the beginning (0:5)
		if (0 <= int(min_main) <= 5):
			for sub_query in parse_of:
				#declare query coordinate for sub hit
				split_liste = sub_query.split("-")
				sub = split_liste[0]
				sub_scaffold = parse_of[sub_query]['scaffold']
				sub_coords = parse_of[sub_query]['coord_query']
				min_sub = sub_coords[0]
				max_sub = sub_coords[1]
				#deaclare all subject coordinates: Master and sub hit
				main_subjects = parse_of[main_query]['coord_subject']
				min_Smain = main_subjects[0]
				max_Smain = main_subjects[1]
				sub_subjects = parse_of[sub_query]['coord_subject']
				min_Ssub = sub_subjects[0]
				max_Ssub = sub_subjects[1]
				#proceed of same query is on same scaffold
				if (sub == main) and (main_scaffold == sub_scaffold):			
				#Now look if on the same scaffolds, there is a hit following the coordinates (with window)
				
			
		
#Based on coordinate, find partial hits and reconstruct
def find_part(window, sub, min_main, max_main, min_sub, max_sub, min_Smain, max_Smain, min_Ssub, max_Ssub):
	
	#declare the chain to return, modified later if found complementary partial hits
	chain = "Unfound"
	#Check if on subject coordinate, a seconde is present
	window_value = float(window)/2
	
	num = int(min_Smain) - int(window_value)
	maxi = int(min_Smain) + int(window_value)
	
	#find a second hit of same query who begin +- at the end of main one (-10 : main)
	if (num <= int(min_sub) <= maxi):
		#re-organise sequence based on subject coordinate : minimam min value eq first sequence
		
		#check for subject orientation
		if (min_Smain > max_Smain) and (min_Ssub > max_Ssub):
			#return the block, whole reverse orientation
			tmp = min_Smain
			min_Smain = max_Smain
			max_Smain = tmp
			tmp = min_Ssub
			min_Ssub = max_Ssub
			max_Ssub = min_Ssub
						
		#Order hits
		if (min_Smain < min_Ssub):
			#case1 
			#check adjacent location
			maxi1 = int(min_Ssub) + int(window_value)
			num1 = int(min_Ssub) - int(window_value)
			#check if the hits on subject are in the same window range as on query coordinates
			if (num1 <= int(max_Smain) <= maxi1):
				end = int(min_sub) - 1
				chain = (str(sub) + "\t" + str(min_sub) + "-" + str(max_sub) + "\t" + str(min_main) + "-" + str(end) + "\n")
		elif (min_Ssub < min_Smain):
			#case2
			#check adjacent location
			maxi2  = int(max_Ssub) + int(window_value)
			num2 = int(max_Ssub) - int(window_value)
			if (num2 <= int(min_Smain) <= maxi2):
				end = int(min_sub) - 1
				chain = (str(sub) + "\t" + str(min_sub) + "-" + str(max_sub) + "\t" + str(min_main) + "-" + str(end) + "\n")	
				
	return(chain)
			
		
		
			
def get_master(master, parse_of):
	
	
	



if __name__ == "__main__": 
	blast6 = '/media/vol2/scratch/lcornet/cassava/CIDER/script/output/rfile2.fasta_on_genomedb.blastn6'
	test = '/media/vol2/scratch/lcornet/cassava/CIDER/script/test/test.blastn'
	reorder_sequence(blast6)








