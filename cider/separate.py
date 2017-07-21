"""separate.py

this module is part of the ciderseq distribution.

"""
import os
import uuid
import tempfile
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

def separate(settings,seq,logger):
	logger.debug("separate > "+str(seq.id))

	#define return variable as nohit
	genome='nohit'
	
	#init blast
	if not settings['blastinit']=='':
		exec(settings['blastinit'],globals())

	#create temporary input/output directory
	tmpdir=tempfile.mkdtemp()

	#define blast base file
	blastfile=str(tmpdir)+"/"+str(uuid.uuid4())
	
	#write blast input file
	with open(blastfile+'.in', "wt") as output_handle:
		SeqIO.write(seq,output_handle,'fasta')
	output_handle.close()
	
	#create blast command
	cline=NcbiblastnCommandline(
			   query = blastfile+".in" 
			,     db = settings['blastndb']
			, evalue = float(settings['evalue'])
			,    out = blastfile+".out"
			,    cmd = settings['blastexe']
			,   task = 'blastn'
			, outfmt = 5) #xml output
	logger.debug(cline)
	#execute command
	stout, stderr = cline()
	#read xml blast result 
	for b_record in NCBIXML.parse(open(blastfile+'.out', 'r')):
		#do we have a blast hit
		if len(b_record.alignments) > 0:
			genome=b_record.alignments[0].title.split(" ")[1]

	#delete temporary input/output files
	os.remove(blastfile+'.in')
	os.remove(blastfile+'.out')
	#remove temporary directory
	os.rmdir(tmpdir)
	#open outputfile
	outputfile=settings['outputdir']+"/"+str(seq.id)+"."+genome+".separate"
	fout = open(outputfile+".fa",'wt')
	#write sequence
	SeqIO.write(seq,fout,'fasta')
	#close file
	fout.close()
	#return genome
	logger.debug("separate = "+genome)
	return genome,outputfile
	
		    