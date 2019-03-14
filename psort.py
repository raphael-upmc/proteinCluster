#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor,wait
import argparse
import shutil


''' this script parallelizes the psort softare  '''

def runPsort(filename,directory,option) :
    fasta_filename = directory+'/'+'psort_tmp/fasta'+'/'+filename
    psortb_filename = directory+'/'+'psort_tmp/output'+'/'+filename.replace('.faa','.psortb')
    log_filename = directory+'/'+'psort_tmp/log'+'/'+filename.replace('.faa','.log')
    cmd = '/usr/bin/psort --output long '+option+' '+fasta_filename+' >'+psortb_filename+' 2>'+log_filename        
    print(cmd)
    status=os.system(cmd)
    return status




def parsePsort(output_filename,directory) :
    print('writting results in '+output_filename)
    tag_header = 0
    output = open(output_filename,'w')
    for root, dirs, files in os.walk(directory):
        for filename in files :
            file = open(root+'/'+filename,'r')
            header = next(file).rstrip().split('\t')
            if tag_header == 0 :
                output.write('\t'.join(header)+'\n')
                tag_header = 1
                
            for line in file :
                line = line.rstrip()
                liste = line.split('\t')
                output.write('\t'.join(liste)+'\n')
            file.close()
    output.close()





if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='this script parallelize the psort softare')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILENAME')
    parser.add_argument('output_filename', help='the path of the OUTPUT_FILENAME')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs (default: 1)')
    parser.add_argument('--option',type=str,default='gram+',help='archaea|gram+|gram- (default: gram+)')

    args = parser.parse_args()

    arg2option = {'archaea' : '--archaea' , 'gram+' : '--positive' , 'gram-' : '--negative' }
    if args.option not in arg2option :
        sys.exit(args.option+' wrong argument, exit')
    else:
        option = arg2option[ args.option ]
        
    if os.path.exists(args.fasta_filename) :
        fasta_filename = os.path.abspath(args.fasta_filename)
    else:
        sys.exit(args.fasta_filename+' does not exist, exit')

    directory = '/'.join( os.path.abspath(args.output_filename).split('/')[:-1] )
    if not os.path.exists(directory) :
        sys.exit(directory+' does not exist')



    cpu = args.cpu
    output_filename = os.path.abspath( args.output_filename )

    if os.path.exists(output_filename) :
        sys.exit(output_filename+' already exists, remove it first')

        
    if os.path.exists(directory+'/'+'psort_tmp') :
        sys.exit(directory+'/'+'psort_tmp'+' already exists, remove it first')
    else :
        os.mkdir(directory+'/'+'psort_tmp')
        os.mkdir(directory+'/'+'psort_tmp/fasta')
        os.mkdir(directory+'/'+'psort_tmp/log')
        os.mkdir(directory+'/'+'psort_tmp/output')
    

    print('fasta filename: '+fasta_filename)
    print('output filename: '+output_filename)
    print('cpu: '+str(cpu))
    print('working directory :'+directory)
    print('option: '+option)
    
    ###################
    # splitting fasta #
    ###################

    nb = 0
    for seq_record in SeqIO.parse(fasta_filename , "fasta") :
        nb += 1
        
    binSeq = ( nb/cpu ) + 1
    cpt = 0
    binNb = 1
    output_fasta_filename = directory+'/'+'psort_tmp/fasta'+'/'+str(binNb)+".faa"
    output = open(output_fasta_filename,"w")
    for seq_record in SeqIO.parse(fasta_filename , "fasta") :
        cpt += 1
        if cpt < binSeq : 
            SeqIO.write(seq_record,output,"fasta")
        else :
            output.close()
            binNb += 1
            cpt = 1
            output_fasta_filename = directory+'/'+'psort_tmp/fasta'+'/'+str(binNb)+'.faa'
            output = open(output_fasta_filename,"w")
            SeqIO.write(seq_record,output,"fasta")
    output.close()


    
    #################
    # running psort #
    #################

    results = list()
    pool = ProcessPoolExecutor(cpu) # start 20 worker processes and 1 maxtasksperchild in order to release memory
    for root, dirs, files in os.walk(directory+'/'+'psort_tmp/fasta'):
        for filename in files :
            future = pool.submit( runPsort,filename,directory,option )
            results.append( future )
    wait(results) # Prevents any more tasks from being submitted to the pool

    for elt in results :
        print( elt.result() )

    pool.shutdown()
    
    ####################
    # writting results #
    ####################

    parsePsort(output_filename,directory+'/'+'psort_tmp/output')

    ##########################
    # removing tmp directory #
    ##########################
    
#    print('removing tmp directory '+directory+'/'+'psort_tmp')
#    shutil.rmtree(directory+'/'+'psort_tmp')
