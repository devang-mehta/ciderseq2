"""phase.py

this module is part of the ciderseq distribution.

"""
#read resultfiles
import sys,os
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna

def phase(settings,genome,deconfile,annofile,logger):
	#check if genome exists
	if genome=='' or not genome in settings['phasegenomes']:
		print('\nphase genome missing or does not exists in config-file: '+genome+'\n')
		sys.exit(1)
	#check if files to process exists
	if not os.path.isfile(deconfile) or not os.path.isfile(annofile):
		print('\nresult file or sequence file missing or does not exists : '+deconfile+"/"+annofile+'\n')
		sys.exit(1)
	basefilename=os.path.splitext(os.path.basename(deconfile))[0]
	#cleanup
	for fformat in settings['outputformat']:
		output_file=settings['outputdir']+"/"+basefilename+"."+fformat
		#remove file
		if os.path.isfile(output_file):
			os.remove(output_file)

	#read sequence file
	#purpose is that at and of script we write a .gb annotation file
	seqrecord={}
	for record in SeqIO.parse(deconfile, "fasta"):
		seqrecord[record.id]=record
	
	#read json result file
	results={}
	with open(annofile) as result_file:
		results=json.load(result_file)
	
	#loop through result
	for seq in results:
		logger.debug("phase > "+seq)
		#check if all proteins are in result and proteins are correctly stranded.
		complete=1;
		reverted=0;
		for k in sorted(settings['phasegenomes'][genome]['proteins']):
			if not k in results[seq]['proteins']: #check if protein from config exists in result
				complete=0 #incomplete - key missing
			elif not results[seq]['proteins'][k]['strand'] in [-1,1]:
				complete=0 #incomplete - contradicting strands
			else:
				#check strand and set reverse flag
				if int(results[seq]['proteins'][k]['strand']) != int(settings['phasegenomes'][genome]['proteins'][k]['strand']):
					reverted+=1 #increase reverse, at the end reverse must have the length of proteins = reverse all 
		
		if complete==1 and (reverted==0 or reverted==len(settings['phasegenomes'][genome]['proteins'])):
			#all proteins are in result and the in strand 
			#set sequence length
			moleculelength=len(seqrecord[seq])
			
			if reverted!=0:
				#the annotations are wrong in strand, we have to turn everything around
				for k in sorted(settings['phasegenomes'][genome]['proteins']): #do it for every protein
					#variables for the new nucleotide positions
					revstart=0
					revend=0
					
					if int(results[seq]['proteins'][k]['maxnucpos']) < moleculelength:
						#the piece is within the sequence
						revstart = moleculelength - int(results[seq]['proteins'][k]['maxnucpos']) + 1 #verified!!
						revend = moleculelength - int(results[seq]['proteins'][k]['minnucpos']) + 1 #verified!!
					else:
						#the piece overhangs sequence end
						revstart=moleculelength - (int(results[seq]['proteins'][k]['maxnucpos'])-moleculelength) + 1 #verified!!
						revend=moleculelength + (moleculelength - int(results[seq]['proteins'][k]['minnucpos'])) + 1 #verified!!
					
					results[seq]['proteins'][k]['minnucpos']=revstart
					results[seq]['proteins'][k]['maxnucpos']=revend
	
				#sequence is reverted, create reverse complement
				seqrecord[seq].seq=seqrecord[seq].seq.reverse_complement()
				#print "this one"
			
			#sequence cutting position
			#set current position of the phaseto protein
			phasegene_original_startpos=int(results[seq]['proteins'][settings['phasegenomes'][genome]['phaseto']]['minnucpos'])
			
			#correct the phasestart with adjustments if phaseto protein starts not at position 1
			phaseadjustment=(3*(int(results[seq]['proteins'][settings['phasegenomes'][genome]['phaseto']]['minprotpos'])-1)) 
			
			#sequence cut/concat position
			seq_cut_pos = (phasegene_original_startpos-int(settings['phasegenomes'][genome]['offset']))-phaseadjustment
	
			#cut and concat the sequence
			seqrecord[seq].seq = seqrecord[seq].seq[seq_cut_pos:]+seqrecord[seq].seq[0:seq_cut_pos]
	
			#break_point_pos after cut and splice needed for genes located before seq_cut_pos
			break_point_pos = moleculelength - seq_cut_pos
	
			for k in sorted(settings['phasegenomes'][genome]['proteins']):
				newstart=0
				newend=0
				if int(results[seq]['proteins'][k]['minnucpos']) < seq_cut_pos: 
					#everything that is located before new 1 pos
					newstart = int(results[seq]['proteins'][k]['minnucpos']) + break_point_pos 
					newend = int(results[seq]['proteins'][k]['maxnucpos']) + break_point_pos 
				else:
					nop=0
					newstart=int(results[seq]['proteins'][k]['minnucpos']) - seq_cut_pos - 1 #verified!! needed to genbank ouput 0-based/1-based transition
					newend=int(results[seq]['proteins'][k]['maxnucpos']) - seq_cut_pos 
				
				logger.debug(str(k)+":"+str(results[seq]['proteins'][k]['minnucpos'])+":"+str(results[seq]['proteins'][k]['maxnucpos'])+":"+str(newstart)+":"+str(newend))
				
				#add the feature to the sequence
				seqrecord[seq].features.append(SeqFeature(FeatureLocation(newstart,newend)
					,type='gene'
					,strand=settings['phasegenomes'][genome]['proteins'][k]['strand']
					,qualifiers={"locus_tag":k}))
			
		else:
			#remove the record from seqrecord
			if seq in seqrecord:
				logger.debug( "discard key : "+seq+ " (complete="+str(complete)+", reverted="+str(reverted)+")")
				del seqrecord[seq]	
			else:
				logger.debug( "key alreday discarded : "+seq )
			
	#print result
	for fformat in settings['outputformat']:
		output_file=settings['outputdir']+"/"+basefilename+"."+fformat
		#write file
		with open(output_file, "at") as output_handle:
			for seqs in seqrecord:
				seqrecord[seqs].name=seqrecord[seqs].name[-15:] #this is a bit a hack.... 
				seqrecord[seqs].description=seqrecord[seqs].description[-15:]#this is a bit a hack.... 
				seqrecord[seqs].id=seqrecord[seqs].id[-15:]#this is a bit a hack.... 
				seqrecord[seqs].seq=Seq(str(seqrecord[seqs].seq),generic_dna)
				SeqIO.write(seqrecord[seqs],output_handle,fformat)

	return settings['outputdir']+"/"+basefilename
