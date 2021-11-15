#! /usr/bin/env python

""" a wrapper to run signalP on large protein dataset """

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
import argparse
import shutil

def runSignalP(fasta_filename,output_filename,log_filename,option) :
    #cmd = "/data7/proteinfams/SignalP/signalp-4.1/signalp -f short "+option+" "+fasta_filename+" > "+output_filename+" 2>>"+log_filename
    cmd = 'cd /env/cns/proj/agc/home/rmeheust/programs/signalp-5.0b/bin ; ./signalp -format short '+option+' -fasta '+fasta_filename+' -prefix '+output_filename+' > '+log_filename+' 2>>'+log_filename
    print(cmd)
    status = os.system(cmd)
    return status


if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='a wrapper to run signalP on large protein dataset')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILENAME')
    parser.add_argument('output_filename', help='the path of the OUTPUT_FILENAME')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs (default: 1)')
    parser.add_argument('--option',type=str,default='gram+',help='euk|gram+|gram-|arc (default: gram+)')

    args = parser.parse_args()

    arg2option = {'euk' : '-org euk' , 'gram+' : '-org gram+' , 'gram-' : '-org gram-' , 'arc' : '-org arch' }
    if args.option not in arg2option :
        sys.exit(args.option+' wrong argument, exit')
    else:
        option = arg2option[ args.option ]

    if os.path.exists(args.fasta_filename) :
        initial_fasta_filename = os.path.abspath(args.fasta_filename)
    else:
        sys.exit(args.fasta_filename+' does not exist, exit')

    cpu = args.cpu
    final_output_filename = os.path.abspath( args.output_filename )
    directory = '/'.join( os.path.abspath(args.output_filename).split('/')[:-1] )
    
    if os.path.exists(directory+'/'+'signalP_tmp') :
        sys.exit(directory+'/'+'signalP_tmp already exists, remove it first')
    else :
        os.mkdir(directory+'/'+'signalP_tmp')
        os.mkdir(directory+'/'+'signalP_tmp/fasta')
        os.mkdir(directory+'/'+'signalP_tmp/log')
        os.mkdir(directory+'/'+'signalP_tmp/output')
    
    print('fasta filename: '+initial_fasta_filename)
    print('output filename: '+final_output_filename)
    print('cpu: '+str(cpu))
    print('working directory :'+directory)
    print('option: '+option)

   
    ###################
    # splitting fasta #
    ###################
    
    # Signalp allow no more than 10000 sequences as inpout, so we have to split the dataset into fasta bins of <10,000 sequences !

    cpt2defline = dict()
    binSeq = 10000
    cpt = nb = 0
    binNb = 1
    output_fasta_filename = directory+'/'+"signalP_tmp/fasta"+"/"+str(binNb)+".faa"
    output = open(output_fasta_filename,"w")
    for seq_record in SeqIO.parse(initial_fasta_filename , "fasta") :
        cpt += 1
        nb += 1
        cpt2defline[ str(nb) ] = seq_record.id
        if cpt < binSeq : 
            SeqIO.write( SeqRecord(seq=seq_record.seq,id=str(nb),description="") ,output,"fasta")
        else :
            output.close()
            binNb += 1
            cpt = 1
            output_fasta_filename = directory+'/'+"signalP_tmp/fasta"+"/"+str(binNb)+".faa"
            output = open(output_fasta_filename,"w")
            SeqIO.write( SeqRecord(seq=seq_record.seq,id=str(nb),description="")  ,output,"fasta")
    output.close()


    ###################
    # running signalP #
    ###################


    results = list()
    pool = mp.Pool(processes=cpu,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory
    for root, dirs, files in os.walk(directory+'/'+"signalP_tmp/fasta"):
        for filename in files :
            fasta_filename = root+'/'+filename
            #        print(fasta_filename)
            output_filename = directory+'/'+"signalP_tmp/output"+'/'+filename.replace('.faa','.signalP')
            log_filename = directory+'/'+"signalP_tmp/log"+'/'+filename.replace('.faa','.log')
            results.append( pool.apply_async( runSignalP, args= (fasta_filename,output_filename,log_filename,option,) ))
    pool.close() # Prevents any more tasks from being submitted to the pool
    pool.join() # Wait for the worker processes to exit



    ####################
    # writting signalP #
    ####################

    output = open(final_output_filename , 'w' )
    
    for root, dirs, files in os.walk(directory+'/'+"signalP_tmp/output"):
        for filename in files :
            file = open(root+'/'+filename,"r")
            for line in file :
                line = line.rstrip()
                if re.match("#",line) :
                    continue
                else :
                    liste = line.split()
                    defline = cpt2defline[str(liste[0])]
                    result = liste[1]
                    output.write(defline+"\t"+result+"\n")
            file.close()        
    output.close()

    
    ##########################
    # removing tmp directory #
    ##########################
    
    print('removing tmp directory '+directory+'/'+'signalP_tmp')
    shutil.rmtree(directory+'/'+'signalP_tmp')

    
    sys.exit()
