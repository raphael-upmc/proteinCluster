#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import shutil
import argparse

def runPRISM(fasta_filename,output_directory,log_filename) :
    cmd = 'java -jar /home/meheurap/programs/prism-releases-master/prism.jar -a -p -f '+fasta_filename+' -tt -sug -res -rib -w 10000 -r /home/meheurap/programs/prism-releases-master/prism/WebContent/ -o '+output_directory+' -web >'+log_filename
    print(cmd)
    status = os.system(cmd)

    return status

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='this script parallelize the psort softare')
    parser.add_argument('fastaList_filename', help='the path of the FASTALIST_FILENAME which contains the list of genomes paths to analyzed')
    parser.add_argument('output', help='the path of the OUTPUT directory')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs (default: 6)')


    args = parser.parse_args()

    if os.path.exists(args.fastaList_filename) :
        fastaList_filename = os.path.abspath(args.fastaList_filename)
    else:
        sys.exit(args.fastaList_filename+' does not exist, exit')


    directory = os.path.abspath(args.output)
    if os.path.exists(directory) :
        sys.exit(directory+' already exists, remove it first')
    else :
        os.mkdir(directory)
        os.mkdir(directory+'/'+'logs')
        os.mkdir(directory+'/'+'results')

    try :
        cpu = int(args.cpu)
    except :
        cpu = 6


    print('fasta filenames list: '+fastaList_filename)
    print('output directory: '+directory)
    print('cpu: '+str(cpu))


    filenameList = set()
    file = open(fastaList_filename,'r')
    for line in file :
        fasta_filename = line.rstrip()
        if not os.path.exists( os.path.abspath(fasta_filename) ) :
            print(fasta_filename+' does not exist, skipped')
        else :
            filenameList.add( os.path.abspath(fasta_filename) )
    file.close()
    
    sys.exit()

    #################
    # running PRISM #
    #################
    
    results = list()
    pool = mp.Pool(processes=cpu,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory
    for fasta_filename in filenameList :
        output_directory = directory+'/'+'results'+'/'+fasta_filename
        log_filename = directory+'/'+'logs'+'/'+fasta_filename+'.log'
        
        if os.path.exists(output_directory+'/'+'geneclusters.txt') :
            print(output_directory+' already exists, skipped this genome')
            continue
        else :
            if os.path.exists(output_directory) :
                shutil.rmtree(output_directory)
            os.mkdir(output_directory) # need to create the output directory before running the cmd line ??
            results.append( pool.apply_async( runPRISM, args= (fasta_filename,output_directory,log_filename,) ))
            break
    pool.close() # Prevents any more tasks from being submitted to the pool
    pool.join() # Wait for the worker processes to exit

    print(results)


        
