"""align.py

this module is part of the ciderseq distribution.

"""
import os
import uuid
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC

def align(settings,seq,genome,logger):
	logger.debug("align > "+str(seq.id))
	#check genome in config
	if not genome in settings['targets']:
		#genome does not exists in settings file
		status=0
		print('\ngenome does not exists in config file : '+genome+'\n')
		sys.exit(1)
	elif not os.path.isfile(settings['targets'][genome]):
		#no sequence file for genome
		print('\ngenome sequence file in config file does not exists : '+settings['targets'][genome]+'\n')
		sys.exit(1)
	else:
		#init muscle
		if not settings['muscleinit']=='':
			exec(settings['muscleinit'],globals())
		
		for record in SeqIO.parse(settings['targets'][genome], "fasta"):
			#read target genome file
			baseid = str(record.id)
			nbaseid = 'n_'+baseid
			myseq = str(record.seq)
			
			#calculate delimiter or create nseq = sequence with half in front and at end
			coord= int(len(myseq)/2)
			nseq=myseq[coord:]+myseq+myseq[0:coord]
			#create temporary input/output directory
			tmpdir=tempfile.mkdtemp()
			# define file
			musclefile = tmpdir+"/"+str(uuid.uuid4())
			#define in 
			in_file = musclefile+".in"
			# define out file
			out_file = musclefile+".out"

			# write in_file fasta file
			fout=open(in_file,'wt')
			fout.write('>'+str(seq.id)+'\n') #submitted record id
			fout.write(str(seq.seq)+'\n') #submitted record sequence
			fout.write('>3_'+baseid+'\n') #target genome id
			fout.write(myseq+myseq+myseq+'\n') #target genome sequence * 3
			fout.write('>'+nbaseid+'\n') #ntarget genome id
			fout.write(nseq+'\n') #ntarget sequence 
			fout.close()
			#execute muscle
			muscle_exe = settings['muscleexe']
			muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
			logger.debug(muscle_cline)
			stdout, stderr = muscle_cline()	
			#read result
			resseq={} #fill result into dict to have direct access to seqrecord
			#read out_file file
			for resrecord in SeqIO.parse(out_file, "fasta"):
				resseq[str(resrecord.id)]=str(resrecord.seq)

			#delete temporary input/output files
			os.remove(in_file)
			os.remove(out_file)
			#remove temporary directory
			os.rmdir(tmpdir)
		
			#perform sliding window
			startpos=0
			endpos=len(resseq[seq.id])
			for i in range(0,(len(resseq[seq.id])-9)):
				if resseq[seq.id][i:(i+settings['windowsize'])].count('-') <= 1:
					startpos = i
					break
					
			for i in reversed(range(0,(len(resseq[seq.id])-(settings['windowsize']-1)))):
				if resseq[seq.id][i:(i+settings['windowsize'])].count('-') <= 1:
					endpos = (i+settings['windowsize'])
					break
			
			result_seq_r=SeqRecord(Seq(resseq[seq.id][startpos:endpos].replace('-',''),IUPAC.unambiguous_dna)
					,id=seq.id,description='')

			#open outputfile
			outputfile=settings['outputdir']+"/"+str(seq.id)+"."+genome+".align"
			fout = open(outputfile+".fa",'wt')
			#write sequence
			SeqIO.write(seq,fout,'fasta')
			#close file
			fout.close()

			return result_seq_r,outputfile
	