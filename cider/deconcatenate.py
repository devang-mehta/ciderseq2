"""deconcatenate.py

this module is part of the ciderseq distribution.

"""
#this script deconcats only one single sequence
import os
import sys
import uuid
import tempfile
import logging
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
#################################################################################################
#consensus score
def _conscore(consensus):
	l = consensus.split('#') #split consensus by #
	#print l
	m=[]
	for i in range(0,len(l)): #loop through list
		if len(l[i])>0: #eliminate all no match
			m.append(l[i]) #keep consenses
	
	#print m
	t=0
	for i in range(0,len(m)): #calculate total length of matches
		t+=len(m[i])
	
	score=0	
	if float(len(m)) > 0:
		score=float(t)/float(len(m)) #calculate score of total length devided by quantity of pieces.	
	return score	
#################################################################################################
#uncirle sequence	
def _uncircle_seq(settings,seq,logger):
	
	divisor = int(round(float(len(str(seq.seq)))/float(settings['fragmentsize']))) #compute the number or pieces
	
	if divisor == 1:
		#pieces too small, interrupt process
		return 0, 0, -1, '','',''
	
	#divide sequence		
	coord= int(len(str(seq.seq))/divisor)
	
	if coord < settings['fragmentsize']:  #check coord size, sould be > fragment size or
		coord = settings['fragmentsize'] #set coord = fragment size
	
	seqstart={}
	seqend={}
	seqcon={}
	seqscore={}
	
	maxscore=0
	maxindex=0
	maxrevidx=0
	
	for fwd_rev_flag in range(0,2): #0=fwd, 1=rev
		#print fwd_rev_flag
		for i in range(1,divisor): #i is counter
			ncoord=coord*i
			startseq=str(seq.seq)[0:ncoord]
			#print "*"*60
			#print startseq
			if fwd_rev_flag == 1:
				tmp=Seq(startseq)
				startseq=str(tmp.reverse_complement())	
				#print startseq
			
			#print "*"*60
			endseq=str(seq.seq)[ncoord:]
			#print ">0:"+str(ncoord)

			myseqs=[]
			myseqs.append(SeqRecord(Seq(startseq,IUPAC.unambiguous_dna)
					,id='start_'+str(seq.id),name='',description=''))
			myseqs.append(SeqRecord(Seq(endseq,IUPAC.unambiguous_dna)
					,id='end_'+str(seq.id),name='',description=''))

			#create id for in filename
			filename="0-"+str(ncoord)+"-"+str(uuid.uuid4());
			in_file = settings['tempdirdir']+"/"+filename+".in"
			# write in_file fasta file
			with open(in_file, "wt") as output_handle:
				SeqIO.write(myseqs,output_handle,'fasta')
			output_handle.close()
			#define output filename			
			out_file = settings['tempdirdir']+"/"+filename+".out"
			#define muscle exe			
			muscle_exe = settings['muscleexe']
			#create muscle command
			muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
			#log
			#logger.debug(muscle_cline)
			#execute muscle command
			stdout, stderr = muscle_cline()	
			resseq={} #fill result into hash
			#read out_file file
			for resrecord in SeqIO.parse(out_file, "fasta"):
				resseq[str(resrecord.id)]=str(resrecord.seq)
				#print str(resrecord.id)[0:20]+"\t"+str(resrecord.seq)
	
			#create consensus
			resseq['consensus']=''
			#keep longest peaces
			
			#sys.stdout.write("consensus...........\t")
			for ii in range(0,len(resseq['start_'+str(seq.id)])):
				#print resseq['start_'+str(seq.id)][i:(i+1)],resseq['end_'+str(seq.id)][i:(i+1)]
				
				if resseq['start_'+str(seq.id)][ii:(ii+1)]==resseq['end_'+str(seq.id)][ii:(ii+1)]:
					#sys.stdout.write("*")
					resseq['consensus']+="*"
				else:
					#sys.stdout.write(" ")
					resseq['consensus']+="#"
		
			#remove files
			os.remove(in_file)
			os.remove(out_file)
			
			seqscore[str(i)+'-'+str(fwd_rev_flag)]=_conscore(resseq['consensus'])
			seqstart[str(i)+'-'+str(fwd_rev_flag)]=resseq['start_'+str(seq.id)]
			seqcon[str(i)+'-'+str(fwd_rev_flag)]=resseq['consensus']
			seqend[str(i)+'-'+str(fwd_rev_flag)]=resseq['end_'+str(seq.id)]
			
			if seqscore[str(i)+'-'+str(fwd_rev_flag)] > maxscore:
				logger.debug(str(i)+":"+str(fwd_rev_flag)+":"+str(seqscore[str(i)+'-'+str(fwd_rev_flag)]))
				maxscore=seqscore[str(i)+'-'+str(fwd_rev_flag)]
				maxindex=i
				maxrevidx=fwd_rev_flag
		
	return (maxindex,maxrevidx,seqscore[str(maxindex)+'-'+str(maxrevidx)],seqstart[str(maxindex)+'-'+str(maxrevidx)],seqcon[str(maxindex)+'-'+str(maxrevidx)],seqend[str(maxindex)+'-'+str(maxrevidx)])
#################################################################################################
def _find_start_end(sequence):
	#find position where we have 9 out of 10 bases
	startpos=0
	endpos=len(sequence)
	#from start
	for i in range(0,(len(sequence)-9)):
		if sequence[i:(i+10)].count('-') <= 1:
			startpos = i
			break
	#from end
	for i in reversed(range(0,(len(sequence)-9))):
		if sequence[i:(i+10)].count('-') <= 1:
			endpos = (i+10)
			break
	#return startpos and endpos
	return (startpos,endpos)	
#################################################################################################
def _prep_new_seq(sstart,send,revcomp,logger):
	#routine to return new sequence
	startpos_start=len(sstart)
	endpos_start= 0

	startpos_end=len(send)
	endpos_end=0
	
	startpos_start,endpos_start = _find_start_end(sstart)
	startpos_end,endpos_end = _find_start_end(send)
	
	logger.debug('find_start_end:'+str(startpos_start)+":"+str(endpos_start)+":"+str(startpos_end)+":"+str(endpos_end))
	mycase=''
	if startpos_start < startpos_end:
		#startseq overlaps endseq
		if endpos_start > endpos_end:
			#endseq is completely in startseq
			sstart=sstart.replace('-','')
			if revcomp==1:
				#case 1d
				logger.debug("case 1d")
				mycase='1d'
				#reverse complement piece, revert it
				tmp=Seq(sstart)
				sstart=str(tmp.reverse_complement()) #create reverse complement
			else:
				#case 1b
				logger.debug("case 1b")
				mycase='1b'
			return [sstart],mycase
		else:
			if revcomp==0:
				#endseq overlaps end of startseq
				#continue with both sequences
				logger.debug("case 2")
				mycase='2'
				return [sstart.replace('-',''),send.replace('-','')],mycase
			else:
				#the start piece is reverse complement
				logger.debug("case 4")
				mycase='4'
				#cut start piece
				nsstart=sstart[startpos_start:startpos_end] #cutted start-sequence
				#cut end piece
				nsend=send[startpos_end:endpos_end] #complete end-sequence
				#revert nsstart
				nsstart=nsstart.replace('-','')
				tmp=Seq(nsstart)
				nsstart=str(tmp.reverse_complement()) #create reverse complement
				#join
				return [(nsstart+nsend).replace('-','')],mycase
	else:
		#endseq overlaps start of startseq
		if endpos_end > endpos_start:
			#startseq is completely in endseq
			#if reverse complement the sstart piece is eliminated, so no further things to do
			if revcomp==1:
				logger.debug("case 1c")
				mycase='1c'
			else:
				logger.debug("case 1a")
				mycase='1a'
			return [send.replace('-','')],mycase
		else:
			if revcomp==0:
				#startseq overlaps end of endseq
				#we have to cut and splice
				logger.debug("case 3")
				mycase='3'
				nsstart=sstart[startpos_start:endpos_start] #complete startseq
				nsend = send[startpos_end:startpos_start]
				return [(nsstart+nsend).replace('-','')],mycase
			else:
				#reverse complement thing
				logger.debug("case 5")
				mycase='5'
				nsstart=sstart[startpos_start:endpos_start] #complete startseq
				#revert nsstart
				nsstart=nsstart.replace('-','')
				tmp=Seq(nsstart)
				nsstart=str(tmp.reverse_complement()) #create reverse complement
				#cut end piece
				nsend=send[startpos_end:startpos_start] #cutted end-sequence
				return [(nsstart+nsend).replace('-','')],mycase
#################################################################################################
def _process_seq(settings,seq,circles,deconcatcases,fout,foutstats,logger):
	logger.debug("process : "+str(seq.id)+" at fragment size "+str(settings['fragmentsize'])+" - round " +str(circles))
	#print seq
	idx, revcomp, circlescore, sstart,scon,send = _uncircle_seq(settings,seq,logger)
	logger.debug("uncircle : "+str(idx)+":"+str(revcomp)+":"+str(circlescore))
	if circlescore > 20:
		logger.debug('recircle')
		nseqs,deconcase = _prep_new_seq(sstart,send,revcomp,logger) #returns array of sequences
		deconcatcases[deconcase]+=1
		logger.debug(deconcatcases)
		counter=0
		for s in nseqs:
			if len(nseqs)>1:
				seq.id=str(seq.id)+str(counter)
			
			#create new sequence record
			newseq=SeqRecord(Seq(s,IUPAC.unambiguous_dna)
					,id=str(seq.id),name='',description='')
			
			_process_seq(settings,newseq,(circles+1),deconcatcases,fout,foutstats,logger)
			counter+=1
	else:
		#write result
		logger.debug("write result : "+ str(seq.id) +" : "+str(circlescore))
		#write seq
		SeqIO.write(seq,fout,'fasta')
		#write statistics
		if settings['statistics']==1:
			foutstats.write(seq.id+"\t"+str(circles))
			foutstats.write("\t"+str(circlescore))
			for c in sorted(deconcatcases):
				foutstats.write("\t"+str(deconcatcases[c]))
			foutstats.write("\n")
#################################################################################################
def deconcatenate(settings,seq,logger):
	logger.debug("deconcat > "+str(seq.id))
	basefilename=settings['outputdir']+"/"+str(seq.id)+".deconcat."+str(settings['fragmentsize'])
	#open outputfile
	outputfile=basefilename+".fa"
	fout = open(outputfile,'wt')
	#open statistic file
	if settings['statistics']==1:
		statfile=basefilename+".stat"
		foutstats = open(statfile,'wt')
	#define stat directory
	deconcatcases={'1a':0,'1b':0,'1c':0,'1d':0,'2':0,'3':0,'4':0,'5':0}
	#load muscle environment if necessary
	if not settings['muscleinit']=='':
		exec(settings['muscleinit'],globals())
	#define tempdir for muscle alignment output
	settings['tempdirdir']=tempfile.mkdtemp()
	#call deconcat first circle
	_process_seq(settings,seq,0,deconcatcases,fout,foutstats,logger)
	#close
	fout.close()
	if settings['statistics']==1:
		foutstats.close() 
	#remove tempdir
	os.rmdir(settings['tempdirdir'])
	return basefilename
#################################################################################################


