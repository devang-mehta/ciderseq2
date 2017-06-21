#!/usr/bin/env python
"""ciderseq.py: Main script of the ciderseq algorithm."""
__author__ = "Matthias Hirsch-Hoffmann, Devang Mehta"
__copyright__ = "Copyright 2017, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Production"
import sys,os
import click
import logging
import tempfile
import json
from Bio import SeqIO

@click.command()
#which arguments and options are needed? 
#arguments=mandatory, option=possible
#*************** ARGUMENETS *****************************
#config-file
@click.argument('configfile', type=click.Path(exists=True,readable=True))
#input-file
@click.argument('inputfile', type=click.Path(exists=True,readable=True))
#*************** OPTIONS *****************************
#input-filetype
@click.option('--format', default='fasta', type=click.Choice(['fasta','fastq','tab','gb']), help='Input-file format (default=FASTA).')
#options to turn processing steps off
@click.option('--no-separation', is_flag=True, help='will not perform separation step.')
@click.option('--no-alignment', is_flag=True, help='will not perform alignment step.')
@click.option('--no-deconcatenation', is_flag=True, help='will not perform de-concatenation step.')
@click.option('--no-annotation', is_flag=True, help='will not perform annotation step.')
@click.option('--no-phasing', is_flag=True, help='will not perform phasing step.')

#plot graph
#https://matplotlib.org/examples/statistics/histogram_demo_histtypes.html

def main(configfile,inputfile,format,no_separation,no_alignment,no_deconcatenation,no_annotation,no_phasing):
	print(configfile)
	print(os.path.dirname(inputfile))
	
	#read config file
	settings={}
	with open(configfile) as config_file:
		try:
			settings=json.load(config_file)	
		except ValueError:
			print("\nError in config file.\n")
			sys.exit(1)

	#validate config file!!!!!!


	#add input-filepath to outputdir settings
	settings['outputdir']=os.path.dirname(inputfile)+"/"+settings['outputdir']
	#define output directories for steps
	steps=['separate','align','deconcat','annotate','phase']
	#tempfiles-dictionary
	tempfiles={}
	for step in steps:
		settings[step]['outputdir']=os.path.dirname(inputfile)+"/"+settings[step]['outputdir']
		#check and create if necessary
		if not os.path.isdir(settings[step]['outputdir']):
			try:
				os.mkdir(settings[step]['outputdir'])
			except:
				print("\nCould not create directory:"+settings[step]['outputdir']+".\n")
				sys.exit(1)
		#create step-file-list
		tempfiles[step]=[] #list for return files
		
#define logger!!!!########################################################################
	#check loglevel
	numeric_level = getattr(logging, settings['loglevel'].upper(), None)
	if not isinstance(numeric_level, int):
		print('\nInvalid log level in config file: ' + settings['loglevel']+"\n")
		sys.exit(1)
	#check outputdirectory
	if not os.path.isdir(settings['outputdir']):
		try:
			os.mkdir(settings['outputdir'])
		except:
			print("\nCould not create directory:"+settings['outputdir']+".\n")
			sys.exit(1)
	#define log-filename
	logfilename=settings['outputdir']+"/"+os.path.basename(inputfile)+".log"
	#open logfile
	logging.basicConfig(format='%(asctime)s %(message)s',filename=logfilename,level=numeric_level)
	logger = logging.getLogger()
#define logger end!!!!########################################################################
	#process input file, sequence by sequence
	for record in SeqIO.parse(inputfile, format):
		#update record-id, replace / to _
		record.id=record.id.replace('/','_')
		#log
		logger.debug("ciderseq > "+record.id)
		#.....
		genome='' #for now....
	
		if not no_separation:
			#perform separation step
			from cider.separate import separate
			genome,sepfile=separate(settings['separate'],record,logger)
			tempfiles['separate'].append(sepfile)
		if not genome=='nohit':
			if not no_alignment: 
				#perorm alignment step
				from cider.align import align
				record,alignfile = align(settings['align'],record,genome,logger)
				tempfiles['align'].append(alignfile)
			deconfile=''
			if not no_deconcatenation:
				#perform deconcaatenation step
				from cider.deconcatenate import deconcatenate
				deconfile=deconcatenate(settings['deconcat'],record,logger)
				tempfiles['deconcat'].append(deconfile)
			else:
				deconfile=inputfile
			annofile=''
			if not no_annotation:
				#perform annotation step
				from cider.annotate import annotate
				annofile=annotate(settings['annotate'],deconfile+".fa",logger)
				tempfiles['annotate'].append(annofile)
			else:
				annofile=inputfile
				
			if not no_phasing:
				#perform phasing step
				from cider.phase import phase
				phasefile=phase(settings['phase'],genome,deconfile+".fa",annofile,logger)
				tempfiles['phase'].append(phasefile)
				
	#join output-files 
	for step in steps:
		outfilename=settings[step]['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]
		if step=='separate' or step=='align':
			#align,align,deconcat - fasta file
			file_summary(step,tempfiles[step],outfilename+".fa","fasta",'')
		elif step=='deconcat':
			#align,align,deconcat - fasta file
			file_summary(step,tempfiles[step],outfilename+".fa","fasta",'.fa')
			if settings['deconcat']['statistics']==1:
				file_summary('stat',tempfiles[step],outfilename+".stat","",'.stat')
		elif step=='annotate':
			#annotate - json file
			file_summary(step,tempfiles[step],outfilename+".json","json",'')
		elif step=='phase':
			for fformat in settings[step]['outputformat']:
				file_summary(step,tempfiles[step],outfilename+"."+fformat,fformat,"."+fformat)

#will summarize all detail files
def file_summary(step,filelist,outfilename,fformat,ending):
	#open output file handler
	fout = open(outfilename,'wt')
	annotation=[] #list for annotations to append and written at the end
	for f in filelist:
		#read/write fasta file
		if step=='separate' or step=='align' or step=='deconcat' or step=='phase':
			#separate, align, deconcat, phase_fa = fasta file
			for record in SeqIO.parse(f+ending, fformat):
				SeqIO.write(record,fout,fformat)
		#read json file
		elif step=='annotate':
			#annotate - json file
			with open(f+ending) as ffile:
				annotation.append(json.load(ffile))
		elif step=='stat':
			#concat simple ascii txt files
			with open(f+ending) as ffile:
				for line in ffile:
					fout.write(line)
		#remove file
		os.remove(f+ending)
	#dump annotation
	if step=='annotate':
		#dump json output
		json.dump(annotation,fout)
	#close outputfile-handler
	fout.close()

if __name__ == '__main__':
	
	main()
