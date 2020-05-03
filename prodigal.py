#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import argparse
from concurrent.futures import ProcessPoolExecutor,wait

def runningProdigal(input_filename,output_filename,protein_filename,gene_filename) :
#    cmd = '/shared/software/bin/prodigal -i '+input_filename+' -o '+output_filename+' -a '+protein_filename+' -d '+gene_filename+' -m -p single >/dev/null 2>/dev/null'
    cmd = '/shared/software/bin/prodigal -i '+input_filename+' -a '+protein_filename+' -m -p single >/dev/null 2>/dev/null'
    status = os.system(cmd)
    return cmd,status

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='extracting the proteins for a contigs file using prodigal')
    parser.add_argument('directory', help='the path of the DIRECTORY that contains the contigs file')
    parser.add_argument('--cwd', help='the path of the CWD where the store will be stored')
    parser.add_argument('--cpu',type=int,default=6,help='number of CPUs (default: 6)')
    parser.add_argument('--force',action='store_true',default=False,help='force prodigal even if the output folders already exist  (default: False)')
    args = parser.parse_args()

    if os.path.exists(args.directory) :
        directory = os.path.abspath(args.directory)
    else:
        sys.exit(args.directory_filename+' does not exist, exit')


    if args.cwd :
        if os.path.exists(args.cwd) :
            cwd = os.path.abspath(args.cwd)
        else:
            sys.exit(args.cwd+' does not exist, exit')
    else:
        cwd = os.getcwd()

    cpu = args.cpu

    print('genome directory: '+directory)
    print('working directory: '+cwd)
    print('number of CPUs: '+str(cpu))
    print()
    
    if not args.force :
        folder = cwd+"/proteins"
        if os.path.exists(folder) :
            sys.exit(folder+" already exists, remove it first")
        os.mkdir(folder)

        folder = cwd+"/genes"
        if os.path.exists(folder) :
            sys.exit(folder+" already exists, remove it first")
        os.mkdir(folder)

        folder = cwd+"/prodigal"
        if os.path.exists(folder) :
            sys.exit(folder+" already exists, remove it first")
        os.mkdir(folder)



    # parallelizing
    results = list()
    pool = ProcessPoolExecutor(args.cpu) # start 20 worker processes and 1 maxtasksperchild in order to release memory    
    cpt = 0
    for root, dirs, files in os.walk(directory):
        for filename in files :
            cpt += 1
            input_filename = root+'/'+filename
            
            output_filename = cwd+'/'+'prodigal'+'/'+filename+'.prodigal.genes'
            gene_filename = cwd+'/'+'genes'+'/'+filename+'.genes.fna'
            protein_filename = cwd+'/'+'proteins'+'/'+filename+'.proteins.faa'
            if not os.path.exists(protein_filename) :
                future = pool.submit( runningProdigal,input_filename,output_filename,protein_filename,gene_filename)
                results.append(future)
            else:
                continue

    wait(results)

    error = 0
    for elt in results :
        cmd,status = elt.result()
        print(str(status)+'\t'+cmd)
        if status :
            print('\t'+'==>'+' '+'Error')
            error += 1

    pool.shutdown()

