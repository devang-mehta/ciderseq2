"""annotate.py

this module is part of the ciderseq distribution.

"""
import os,sys
import uuid
import tempfile
import json
from Bio import SeqIO
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML


def annotate(settings,inputfile,logger):
	#init blast
	if not settings['blastinit']=='':
		exec(settings['blastinit'],globals())
	#create temporary input/output directory
	tmpdir=tempfile.mkdtemp()
	#this will hold the results	
	rdata={}
	#read file and annotate
	for record in SeqIO.parse(inputfile, 'fasta'):
		logger.debug("annotate > "+str(record.id))
		#create new result entry
		rdata[record.id]={}
		rdata[record.id]['length']=len(record.seq)
		rdata[record.id]['proteins']={}
		#print(record.id,'-',len(record.seq))
	
		#define blast base file
		blastfile=str(tmpdir)+"/"+str(uuid.uuid4())
	
		#write blast input file
		with open(blastfile+'.in', "wt") as output_handle:
			SeqIO.write(record,output_handle,'fasta')
		output_handle.close()
		#build tblastn command
		cline=NcbitblastnCommandline(
				    query = settings['tblastndb']
				, subject = blastfile+".in" 
				,  evalue = settings['evalue']
				,     out = blastfile+".out"
				,    task = settings['blastexe']
				,  outfmt = 5) #xml
		#debug cline command
		logger.debug(cline)
		#execute command
		stout, stderr = cline()
		#dictionary for results
		aligns={}
		#read result
		for b_record in NCBIXML.parse(open(blastfile+'.out', 'r')):
			#read all alignments
			for alignment in b_record.alignments:
				#read all hsp 
				for hsp in alignment.hsps:
					#query-frame and subject-frame
					qframe,sframe = hsp.frame
					#does align key exists in aligns?
					if b_record.query in aligns.keys():
						#append alignment
						aligns[b_record.query].append([sframe,hsp.sbjct_start,hsp.sbjct_end,hsp.query_start,hsp.query_end,hsp.expect])
					else:
						#create key
						aligns[b_record.query]={}
						#all alignment
						aligns[b_record.query]=[[sframe,hsp.sbjct_start,hsp.sbjct_end,hsp.query_start,hsp.query_end,hsp.expect]]
	
		#remove input/output files
		os.remove(blastfile+'.in')
		os.remove(blastfile+'.out')
					
		#evaluate results
		#loop through all keys
		for key in sorted(aligns):
			rdata[record.id]['proteins'][key]={}
			#print(key,len(aligns[key]),aligns[key])
			rdata[record.id]['proteins'][key]['hsps']=aligns[key]
			#def and set result variables
			minprotpos=len(record.seq)+99
			maxprotpos=1
			minnucpos=-1
			maxnucpos=-1
			strand=0
			#loop through aligns
			for i in range(0,len(aligns[key])):
				#extract elements
				f2_frame,f2_start,f2_end,q2_start,q2_end,evalue=aligns[key][i]
				
				#check same frame direction
				if f2_frame < 0 and strand==0:
					#reverse strand and strand not yet set
					strand=-1
				elif f2_frame > 0 and strand==0:
					#forward strand and strand not yet set
					strand=1
				elif f2_frame > 0 and strand>0:
					#ok - forward strand
					nop=0
				elif f2_frame < 0 and strand<0:
					#ok - reverse strand
					nop=0
				else:
					#strand contradicting, break here, the sequence will not be further processed
					#print("strand error")
					strand=-99
					break
				#reverse, we have to switch start and end 
				if f2_frame < 0:
					#switch nuc start-end
					f2_start,f2_end=f2_end,f2_start
				
				#def smallest protein position
				if q2_start < minprotpos:
					minprotpos=q2_start
					#with that, set smalles nucleotide position
					minnucpos=f2_start
				#def largest protein position
				if q2_end > maxprotpos:
					maxprotpos=q2_end
					#with that, set largest nucleotide position
					maxnucpos=f2_end
	
			#strand and protein start-end position are ready
			rdata[record.id]['proteins'][key]['strand']=strand
			rdata[record.id]['proteins'][key]['minprotpos']=minprotpos
			rdata[record.id]['proteins'][key]['maxprotpos']=maxprotpos
	
			#evaluate 
			if strand < 0:
				#reverse minnucposition in result is current maxnucpos
				rdata[record.id]['proteins'][key]['minnucpos']=maxnucpos
				#check if we run over breakpoint
				if minnucpos > maxnucpos:
					#normal piece, just turn minnucpos and maxnucpos
					#print("ok reverse")
					#print(strand,maxnucpos,minnucpos,minprotpos,maxprotpos)		
					rdata[record.id]['proteins'][key]['case']='ok reverse'
					rdata[record.id]['proteins'][key]['maxnucpos']=minnucpos
				else:
					#over breakpoint
					#print("over breakpoint reverse")
					#print(strand,minnucpos,maxnucpos,minprotpos,maxprotpos)		
					#print(strand,maxnucpos,(len(record.seq)+minnucpos),minprotpos,maxprotpos)
					rdata[record.id]['proteins'][key]['case']='over breakpoint reverse'
					rdata[record.id]['proteins'][key]['maxnucpos']=(len(record.seq)+minnucpos)
			else:
				rdata[record.id]['proteins'][key]['minnucpos']=minnucpos
				#check if we run over breakpoint
				if minnucpos > maxnucpos:
					#over breakpoint
					#print("over breakpoint forward")
					#print(strand,minnucpos,maxnucpos,minprotpos,maxprotpos)		
					#print(strand,minnucpos,(len(record.seq)+maxnucpos),minprotpos,maxprotpos)
					rdata[record.id]['proteins'][key]['case']='over breakpoint forward'
					rdata[record.id]['proteins'][key]['maxnucpos']=(len(record.seq)+maxnucpos)
				else:
					#normal piece
					#print("ok forward")
					#print(strand,minnucpos,maxnucpos,minprotpos,maxprotpos)		
					rdata[record.id]['proteins'][key]['case']='ok forward'
					rdata[record.id]['proteins'][key]['maxnucpos']=maxnucpos
			
	#write result
	outputfile=settings['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]
	#print(json.dumps(rdata))
	with open(outputfile+".json",'w') as fout:
		json.dump(rdata,fout)
	fout.close()

	os.rmdir(tmpdir)
	#return genome
	return outputfile
				
