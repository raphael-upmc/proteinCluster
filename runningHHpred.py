#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,wait

''' the script runs the HHpred software on the pfam db and parses the results '''



def runningHhblits(hhm_filename,hhblits_database,hhr_filename) :
    basename = os.path.basename(hhm_filename).split('.')[0]
    cmd = '/home/meheurap/programs/hhsuite-3.0-beta.3-Linux/bin/hhblits -i '+hhm_filename+' -o '+hhr_filename+' -d '+hhblits_database+'  -v 0 -p 50 -E 0.001 -z 1 -Z 32000 -B 0 -b 0 -n 2 -cpu 1'
    status = os.system(cmd)
    if status == 0 :
        return basename,True
    else:
        return basename,False


hhm_directory = '/data7/proteinfams/DPANN/DPANN_protein_clustering_20190405/hhblits/hhm'
hhblits_database = '/data7/proteinfams/HHpred/pfam'


error = 0
results = list()
pool = ProcessPoolExecutor(48) # start 20 worker processes and 1 maxtasksperchild in order to release memory

cpt = 0
for root, dirs, files in os.walk(hhm_directory):
    for filename in files :
        cpt += 1
        hhm_filename = root+'/'+filename
        hhm_basename = os.path.basename(hhm_filename).split('.')[0]
        hhblits_database_basename = os.path.basename(hhblits_database)        
        hhr_filename = '/data7/proteinfams/HHpred/results/'+hhm_basename+'__'+hhblits_database_basename+'.hhr'

        future = pool.submit( runningHhblits,hhm_filename,hhblits_database,hhr_filename )
        results.append(future)
        # if cpt == 100 :
        #     break

    wait(results)

    for elt in results :
        subfamily,result = elt.result()
        if not  result :
            print('\t'+subfamily+' '+'==>'+' '+'Error' )
            error += 1
    pool.shutdown()
    
    if error != 0 :
        print('\n'+str(error)+' hhm failed to run hhblits\n')
