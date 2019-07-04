#!/usr/bin/env python
"""cs-tools.py: Script to split and join reads in input file."""
__author__ = "Matthias Hirsch-Hoffmann, Devang Mehta"
__copyright__ = "Copyright 2019, Devang Mehta & Matthias Hirsch-Hoffmann"
__version__ = "2.0.0"
__maintainer__ = "Devang Mehta"
__email__ = "devangmehta@ualberta.ca"
__status__ = "Production"
import sys,os
import click
import logging
import tempfile
import json
import glob
import time
import numpy as np
from Bio import SeqIO

########################################################################################################
@click.group()
def cli():
    pass
########################################################################################################
# SPLIT ################################################################################################
########################################################################################################
@cli.command('split',help='will split input-file into --numjobs jobs.')
@click.argument('configfile', type=click.Path(exists=True,readable=True))
@click.argument('inputfile', type=click.Path(exists=True,readable=True))
@click.option('--format', default='fastq', type=click.Choice(['fasta','fastq','tab','gb']), help='Input-file format (default=FASTQ).',prompt=True)
@click.option('--numjobs', default=1, type=click.INT , help='number of jobs for input-file.',prompt=True)
@click.option('--cluster', default='', type=click.STRING , help='Cluster submit command.',prompt=True)
@click.option('--mode', default='full', type=click.Choice(['full','deconcat']), help='Type of operation (default=full).',prompt=True)
def split(configfile,inputfile,format,numjobs,cluster, mode):
    if not os.path.exists(inputfile+".dir"):
        os.makedirs(inputfile+".dir")
        #outputfiles dict
        outputfile={}
        for i in range (0,numjobs):
            outputfile[i]={}
            outputfile[i]['ofile']=inputfile+".dir/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+str(i)+".fa"
            outputfile[i]['handle']=open(outputfile[i]['ofile'],'wt')
        i=0 #file-counter
        #process input file, sequence by sequence
        for record in SeqIO.parse(inputfile, format):
            #write sequence
            SeqIO.write(record,outputfile[i]['handle'],'fasta')
            #click.echo('write to '+outputfile[i]['ofile'])
            i+=1 #increase job-counter
            if i == numjobs:
                i=0 #reset
        #create execution cmd
        for i in range (0,numjobs):
            outputfile[i]['handle'].close()
            if not cluster=='':
                if mode=='full':
                    outputfile[i]['cluster']=cluster +" \"python3 ./ciderseq.py "+configfile+" "+outputfile[i]['ofile']+" --format fasta\""
                else:
                    outputfile[i]['cluster'] = cluster + " \"python3 ./ciderseq.py " + configfile + " " + outputfile[i]['ofile'] + " --format fasta --no-separation --no-alignment --no-annotation --no-phasing\""
            else:
                if mode == 'full':
                    outputfile[i]['cluster']="python3 ./ciderseq.py "+configfile+" "+outputfile[i]['ofile']+" --format fasta &"
                else:
                    outputfile[i]['cluster'] = "python3 ./ciderseq.py " + configfile + " " + outputfile[i]['ofile'] + " --format fasta --no-separation --no-alignment --no-annotation --no-phasing &"
            click.echo(outputfile[i]['cluster'])
        #query execution
        if click.confirm('Do you want to execute above commands?'):
            #execute jobs
            click.echo('submitting jobs..')
            for i in range (0,numjobs):
                click.echo('.')
                os.system(outputfile[i]['cluster'])
                #delay execution on next command because of directory creation in ciderseq.py
                time.sleep(1)
    else:
        print('Output directory already exists. Use \'rm -rf '+inputfile+'.dir\' to remove split-target directory.')
        sys.exit(1)
########################################################################################################
# JOIN #################################################################################################
########################################################################################################
@cli.command('join',help='will join splitted input-file.')
@click.argument('configfile', type=click.Path(exists=True,readable=True))
@click.argument('inputfile', type=click.Path(exists=True,readable=True))
@click.option('--clean', is_flag=True, help='will cleanup split folder.',prompt=True)
@click.option('--mode', default='full', type=click.Choice(['full','deconcat']), help='Type of operation (default=full).',prompt=True)
def join(configfile,inputfile,clean,mode):
    if mode == 'full':
        #print(inputfile+".dir")
        if os.path.exists(inputfile+".dir"):
            #read config-file
            settings=config(configfile)
            #clean-up
            path=os.path.dirname(os.path.abspath(inputfile))
            #print(path)
            fname=os.path.splitext(os.path.basename(inputfile))[0]
            #print(fname)
            #identify genomes
            genomes=[]
            for genome in settings['phase']['phasegenomes']:
                genomes.append(genome)
    ##CLEANUP AND CREATION OF SUMMARY OUTPUT DIRECTORY #####################################################
            #check outputdirectory
            if not os.path.isdir(path+"/"+settings['outputdir']):
                try:
                    os.mkdir(path+"/"+settings['outputdir'])
                except:
                    print("\nCould not create directory:"+path+"/"+settings['outputdir']+".\n")
                    sys.exit(1)
            else:
                #cleanup
                cleanfile=path+"/"+settings['outputdir']+"/"+fname
                clean_join(cleanfile,['log'])
            #create summary folders
            steps=['separate','align','deconcat','annotate','phase']
            for step in steps:
                tmppath=path+"/"+settings[step]['outputdir']
                #check and create if necessary
                if not os.path.isdir(tmppath):
                    try:
                        os.mkdir(tmppath)
                    except:
                        print("\nCould not create directory:"+tmppath+".\n")
                        sys.exit(1)
                else:
                    #cleanup
                    cleanfile=path+"/"+settings[step]['outputdir']+"/"+fname
                    ftype=[]
                    #clean up steps
                    if step=='separate':
                        ftype=['nohit.fa']
                        for genome in genomes:
                            ftype.append(genome+'.fa')
                        clean_join(cleanfile,ftype)
                    elif step=='align':
                        for genome in genomes:
                            ftype.append(genome+'.fa')
                        clean_join(cleanfile,ftype)
                    elif step=='deconcat':
                        for genome in genomes:
                            ftype.append(genome+'.fa')
                            ftype.append(genome+'.stat')
                        clean_join(cleanfile,ftype)
                    elif step=='annotate':
                        for genome in genomes:
                            ftype.append(genome+'.json')
                        clean_join(cleanfile,ftype)
                    elif step=='phase':
                        for genome in genomes:
                            for fformat in settings[step]['outputformat']:
                                ftype.append(genome+'.'+fformat)
                        clean_join(cleanfile,ftype)

    ##COLLECT AND CREATE SUMMARY FILES #####################################################################
            #collect all fasta files in split-directory, as we don't know how
            #many were created during split
            files=glob.glob(inputfile+".dir/*.fa")
            for ffile in files:
                #print(ffile)
                #set outputfilename
                fname=os.path.splitext(os.path.basename(ffile))[0]
                #print("fname is: "+fname)
                for step in steps:
                    if step=='separate':
                        join_text(inputfile,fname+".fa",settings,"log")
                    elif step=='align':
                        join_SeqIO(inputfile,fname,settings['separate'],'nohit.fa','fasta')
                        for genome in genomes:
                            join_SeqIO(inputfile,fname,settings['separate'],genome+'.fa','fasta')
                    elif step=='deconcat':
                        for genome in genomes:
                            join_SeqIO(inputfile,fname,settings['deconcat'],genome+'.fa','fasta')
                            join_text(inputfile,fname,settings['deconcat'],genome+'.stat')
                    elif step=='annotate':
                        for genome in genomes:
                            join_json(inputfile,fname,settings['annotate'],genome+'.json')
                    elif step=='phase':
                        for genome in genomes:
                            for fformat in settings['phase']['outputformat']:
                                join_SeqIO(inputfile,fname,settings['phase'],genome+'.'+fformat,fformat)
    ### REMOVE SPLIT DIRECOTRY OF FLAG WAS SET #############################################################
            if clean:
                #remove split-directory
                for root, dirs, files in os.walk(inputfile+".dir", topdown=False):
                    for name in files:
                    #print(os.path.join(root, name))
                        os.remove(os.path.join(root, name))
                    for name in dirs:
                        #print(os.path.join(root, name))
                        os.rmdir(os.path.join(root, name))
                os.rmdir(inputfile+".dir")
        else:
            print('Output directory missing.')
        sys.exit(1)
    else:
        print(inputfile + ".dir")
        if os.path.exists(inputfile + ".dir"):
            # read config-file
            settings = config(configfile)
            # clean-up
            path = os.path.dirname(os.path.abspath(inputfile))
            print(path)
            fname = os.path.splitext(os.path.basename(inputfile))[0]
            print(fname)
            # identify genomes
            genomes = []
            for genome in settings['phase']['phasegenomes']:
                genomes.append(genome)
            ##CLEANUP AND CREATION OF SUMMARY OUTPUT DIRECTORY #####################################################
            # check outputdirectory
            if not os.path.isdir(path + "/" + settings['outputdir']):##output directory is logs
                try:
                    os.mkdir(path + "/" + settings['outputdir'])
                    print(path + "/" + settings['outputdir'])
                except:
                    print("\nCould not create directory:" + path + "/" + settings['outputdir'] + ".\n")
                    sys.exit(1)
            else:
                # cleanup
                cleanfile = path + "/" + settings['outputdir'] + "/" + fname
                print('cleanfile1 is'+' '+cleanfile)
                clean_join(cleanfile, ['log'])
            # create summary folders

            step = 'deconcat'
            tmppath = path + "/" + settings[step]['outputdir']
                # check and create if necessary
            print("tmppath is"+tmppath)
            if not os.path.isdir(tmppath):
                try:
                        os.mkdir(tmppath)
                except:
                        print("\nCould not create directory:" + tmppath + ".\n")
                        sys.exit(1)
            else:
                # cleanup
                cleanfile = path + "/" + settings[step]['outputdir'] + "/" + fname
                print('cleanfle2'+" "+cleanfile)
                ftype = []
                # clean up steps
                if step == 'deconcat':

                    ftype.append('.fa')
                    ftype.append('.stat')
                    clean_join(cleanfile, ftype)
                
            ##COLLECT AND CREATE SUMMARY FILES #####################################################################
            # collect all fasta files in split-directory, as we don't know how
            # many were created during split
            files = glob.glob(inputfile + ".dir/*.fa")
            for ffile in files:

                # set outputfilename
                fname = os.path.splitext(os.path.basename(ffile))[0]
                print(fname)

                print(step)
                if step == 'deconcat':
                    join_SeqIO(inputfile, fname, settings['deconcat'],'.fa', 'fasta')
                    join_text(inputfile, fname, settings['deconcat'],'.stat')
                else:
                    print("error")

            ### REMOVE SPLIT DIRECOTRY IF FLAG WAS SET #############################################################
            if clean:
                # remove split-directory
                for root, dirs, files in os.walk(inputfile + ".dir", topdown=False):
                    for name in files:
                        # print(os.path.join(root, name))
                        os.remove(os.path.join(root, name))
                    for name in dirs:
                        # print(os.path.join(root, name))
                        os.rmdir(os.path.join(root, name))
                os.rmdir(inputfile + ".dir")
        else:
            print('Output directory missing.')
        sys.exit(1)
########################################################################################################
# PLOT #################################################################################################
########################################################################################################
@cli.command('plot',help='will create statistical plots on results. needs results from separate and phase step.')
@click.argument('configfile', type=click.Path(exists=True,readable=True))
@click.argument('inputfile', type=click.Path(exists=True,readable=True))


def plot(configfile,inputfile):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    #print('chart')
    #create chart directory
    outputdir=os.path.dirname(inputfile)+"/plots"
    if not os.path.isdir(outputdir):
        try:
            os.mkdir(outputdir)
        except:
            print("\nCould not create directory:"+outputdir+".\n")
            sys.exit(1)
    #read config-file
    settings=config(configfile)
    #identify genomes
    genomes=[]
    for genome in settings['phase']['phasegenomes']:
        genomes.append(genome)
    #for every genome
    for g in genomes:
        ############################ SEQUENCE LENGTH ###############################
        maxlength=0
        #original length
        olength=[]
        fname=os.path.dirname(inputfile)+"/"+settings['separate']['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+g+".fa"
        if os.path.isfile(fname):
            #read result file
            for record in SeqIO.parse(fname, 'fasta'):
                slen = len(record.seq)
                if slen > maxlength:
                    maxlength=slen
                olength.append(slen)
        else:
            print('result file missing:'+fname+"\n")
            sys.exit(1)

        #deconcated length
        dlength=[]
        fname=os.path.dirname(inputfile)+"/"+settings['phase']['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+g+"."+settings['phase']['outputformat'][0]
        if os.path.isfile(fname):
            #read result file
            for record in SeqIO.parse(fname, settings['phase']['outputformat'][0]):
                slen = len(record.seq)
                if slen > maxlength:
                    maxlength=slen
                dlength.append(slen)
        else:
            print('result file missing:'+fname+"\n")
            sys.exit(1)
        #define chart length and qty bins
        maxlength=(int(maxlength/200)+1)*200
        bins=maxlength/200

        fig, ax = plt.subplots()
        oarray=np.asarray(olength)
        darray=np.asarray(dlength)
        plt.xlabel('Sequence length')
        plt.ylabel('Number of sequences')
        plt.title('Sequence length before and after DeConcatenation: '+g)
        t=ax.hist([oarray,darray],int(bins),histtype='bar', range=(0,maxlength),label=['before','after'])
        maxy=0
        for i in t[0]:
            if max(i) > maxy:
                maxy=max(i)+1
        ax.set_yticks(range(int(maxy)))
        ax.legend()
        plt.savefig(outputdir+"/"+os.path.splitext(os.path.basename(inputfile))[0]+'.seqlength.'+g+'.png')   # save the figure to file
        plt.close()

    ############################ Number of Frameshifts ########################################
    #protein list
    proteins=[]
    #define result matrix
    for g in genomes:
        for k in sorted(settings['phase']['phasegenomes'][g]['proteins']):
            proteins.append(k)
    #print(proteins)
    #define result dict
    frameshift={}
    for k in sorted(proteins):
        frameshift[k]=0
    #add value of qty sequences
    qty_seq=0
    #read annotation files
    for g in genomes:
        #generate filename
        fname=os.path.dirname(inputfile)+"/"+settings['annotate']['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+g+".json"
        #check if file exists
        if os.path.isfile(fname):
            #read json result file
            results={}
            with open(fname) as result_file:
                results=json.load(result_file)
                #loop through result
                for seq in results:
                    for id in seq:
                        complete=1
                        reverted=0
                        for k in sorted(settings['phase']['phasegenomes'][g]['proteins']):
                            #print(k)
                            if not k in seq[id]['proteins']: #check if protein from config exists in result
                                complete=0 #incomplete - key missing
                            elif not seq[id]['proteins'][k]['strand'] in [-1,1]:
                                complete=0 #incomplete - contradicting strands
                            else:
                                #check strand and set reverse flag
                                if int(seq[id]['proteins'][k]['strand']) != int(settings['phase']['phasegenomes'][g]['proteins'][k]['strand']):
                                    reverted+=1 #increase reverse, at the end reverse must have the length of proteins = reverse all

                        if complete==1 and (reverted==0 or reverted==len(settings['phase']['phasegenomes'][g]['proteins'])):
                            qty_seq+=1
                            for k in sorted(settings['phase']['phasegenomes'][g]['proteins']):
                                #print k, len(results[seq]['proteins'][k]['hsps'])-1
                                frameshift[k]+=len(seq[id]['proteins'][k]['hsps'])-1
        else:
            print('result file missing:'+fname+"\n")
            sys.exit(1)

    fig, ax = plt.subplots()

    plt.xlabel('Protein')
    plt.ylabel('Number of frameshifts')
    plt.title('Number of frameshifts')

    y=[]
    for p in sorted(proteins):
        y.append(frameshift[p])

    ax.set_yticks(range(max(y)+1))
    ax.set_yticklabels(range(max(y)+1))
    ax.set_xticks(range(len(frameshift)))
    ax.set_xticklabels(sorted(proteins))
    ax.bar(range(len(frameshift)), y, 1/1.5)

    plt.savefig(outputdir+"/"+os.path.splitext(os.path.basename(inputfile))[0]+'.frameshift.png')   # save the figure to file
    plt.close()

    ############################ deconcat-stats ########################################
    #read the stats file
    for g in genomes:
        fname=os.path.dirname(inputfile)+"/"+settings['deconcat']['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+g+".stat"
        if os.path.isfile(fname):
            #read stat-file and process
            qtyseq=0
            deconcatrounds=[]
            for i in range(26):
                deconcatrounds.append(0)

            deconcatscores=[]
            deconcatcases={'1a':0,'1b':0,'1c':0,'1d':0,'2':0,'3':0,'4':0,'5':0}
            with open(fname) as ffile:
                for line in ffile:
                    l=line.replace('\n','').split('\t')
                    qtyseq+=1
                    deconcatrounds[int(l[1])]+=1
                    deconcatscores.append(float(l[2]))
                    deconcatcases['1a']+=int(l[3])
                    deconcatcases['1b']+=int(l[4])
                    deconcatcases['1c']+=int(l[5])
                    deconcatcases['1d']+=int(l[6])
                    deconcatcases['2'] +=int(l[7])
                    deconcatcases['3'] +=int(l[8])
                    deconcatcases['4'] +=int(l[9])
                    deconcatcases['5'] +=int(l[10])
            #rounds ###########################################################################
            fig, ax = plt.subplots()
            plt.xlabel('Number of DeConcat rounds')
            plt.ylabel('Number of sequences')
            plt.title(g)

            ax.set_yticks(range(max(deconcatrounds)+1))
            ax.set_xticks(range(26))
            ax.set_xticklabels(range(26))
            ax.bar(range(26), deconcatrounds, 1/1.5)

            plt.savefig(outputdir+"/"+os.path.splitext(os.path.basename(inputfile))[0]+'.deconcatrounds.'+g+'.png')   # save the figure to file
            plt.close()
            #score ###########################################################################
            fig, ax = plt.subplots()
            plt.xlabel('Final DeConcat Score')
            plt.ylabel('Number of sequences')
            plt.title(g)

            ax.set_xticks(range(26))
            ax.set_xticklabels(range(26))

            t=ax.hist(deconcatscores ,26,histtype='bar',range=[0,25])

            maxy=0
            for i in t[0]:
                #print(i)
                if int(i) >= maxy:
                    #print(i)
                    maxy=(int(i)+1)
                    #print('....'+str(maxy))

            ax.set_yticks(range(int(maxy)))
            ax.set_yticklabels(range(int(maxy)))

            plt.savefig(outputdir+"/"+os.path.splitext(os.path.basename(inputfile))[0]+'.deconcatscore.'+g+'.png')   # save the figure to file
            plt.close()
            #cases ###########################################################################
            fig, ax = plt.subplots()
            plt.xlabel('Alignment Cases')
            plt.ylabel('Number of sequences')
            plt.title(g)

            x=[]
            for p in sorted(deconcatcases):
                x.append(p)
            y=[]
            for p in sorted(deconcatcases):
                y.append(deconcatcases[p])

            ax.set_yticks(range(max(y)+1))
            ax.set_yticklabels(range(max(y)+1))
            ax.set_xticks(range(len(deconcatcases)))
            ax.set_xticklabels(sorted(x))
            ax.bar(range(len(deconcatcases)), y, 1/1.5)

            plt.savefig(outputdir+"/"+os.path.basename(inputfile)+'.deconcatcases.'+g+'.png')   # save the figure to file
            plt.close()
        else:
            print('result file missing:'+fname+"\n")
            sys.exit(1)

    sys.exit(1)
########################################################################################################
# SUB FUNCTIONS
########################################################################################################
def clean_join(inputfile,filetypes):
    #remove all existing summary files
    for ftype in filetypes:
        summaryfile=inputfile+"."+ftype
        #print('remove:'+summaryfile)
        if os.path.isfile(summaryfile):
            #remove file
            os.remove(summaryfile)
########################################################################################################
def join_text(inputfile,fname,settings,fileext):
    resultfile=inputfile+".dir/"+settings['outputdir']+"/"+fname+"."+fileext
    #print(resultfile)
    #print("join_text_inputfile is "+inputfile)
    if os.path.isfile(resultfile):
        #print('read')
        rf = open(resultfile,'r')
        data = rf.read()
        rf.close()
        #print(os.path.splitext(os.path.basename(inputfile))[0])
        #print("fileext is "+fileext)
        summaryfile=os.path.dirname(inputfile)+"/"+settings['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+fileext
        #print(summaryfile)
        #print('write')
        sf = open(summaryfile,"at")
        sf.write(data)
        sf.close()
########################################################################################################
def join_SeqIO(inputfile,fname,settings,fileext,filetype):
    resultfile=inputfile+".dir/"+settings['outputdir']+"/"+fname+"."+fileext
    #print("resultfile is "+resultfile)
    if os.path.isfile(resultfile):
        summaryfile=os.path.dirname(inputfile)+"/"+settings['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+fileext
        #print("summaryfile is"+" "+summaryfile)
        #append to summary file
        with open(summaryfile, "at") as output_handle:
            ##read result file
            for record in SeqIO.parse(resultfile, filetype):
                #write record to summary
                SeqIO.write(record, output_handle, filetype)
    #else:
    #    print("error")
########################################################################################################
def join_json(inputfile,fname,settings,fileext):
    #dictionary for annotations to append and written at the end
    annotation=[]
    #read summary file
    summaryfile=os.path.dirname(inputfile)+"/"+settings['outputdir']+"/"+os.path.splitext(os.path.basename(inputfile))[0]+"."+fileext
    #print(summaryfile)
    if os.path.isfile(summaryfile):
        with open(summaryfile) as ffile:
            tmp=json.load(ffile)
            for r in tmp:
                annotation.append(r)
        ffile.close()
    #print(annotation)
    #append content of resultfile
    resultfile=inputfile+".dir/"+settings['outputdir']+"/"+fname+"."+fileext
    #print(resultfile)
    if os.path.isfile(resultfile):
        with open(resultfile) as ffile:
            tmp=json.load(ffile)
            for r in tmp:
                #print(r)
                annotation.append(r)
        ffile.close()
    #write new summary file
    fout = open(summaryfile,'wt') #write new
    #dump json output
    json.dump(annotation,fout)
    #close outputfile-handler
    fout.close()
########################################################################################################
def config(configfile):
    #check config file
    if configfile is None:
        print("\nConfig file required.\n")
        sys.exit(1)
    #read configfile for output directroy
    settings={}
    with open(configfile) as config_file:
        try:
            settings=json.load(config_file)
        except ValueError:
            print("\nError in config file.\n")
            sys.exit(1)
    return settings
########################################################################################################
########################################################################################################
########################################################################################################
if __name__ == '__main__':
    cli()
########################################################################################################
