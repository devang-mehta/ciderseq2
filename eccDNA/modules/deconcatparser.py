#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""deconcat parser"""

import pickle

#Parsing deconcat stat
#Print if of sequence with > 0 run of deconcat
def deconcat_parser(deconcat_stat):
	
	deconcatfilters = []
	#read stat file
	count = 0
	result_handle = open(deconcat_stat)	
	for line in result_handle:
		deconcat_record = line.replace("\n", "")
		split_liste = deconcat_record.split("\t")
		query = split_liste[0]
		cycle = split_liste[1]
		
		if (int(cycle) > 0):
			deconcatfilters.append(query)
	
	pick_file = 'output/deconcat.pickle';
	with open(pick_file, 'wb') as pickleout:
		pickle.dump(deconcatfilters, pickleout, protocol=pickle.HIGHEST_PROTOCOL)
		
if __name__ == "__main__":
	deconcat_stat = '/home/luc/Documents/Python-course/CIDER-cassava/debug/test-debug/deconcat.stat'
	deconcat_parser(deconcat_stat)
