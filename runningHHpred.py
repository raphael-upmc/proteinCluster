#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,wait
import json
import argparse

''' the script runs the HHpred software on the pfam db and parses the results '''

def subfamily2familyFunction(tsv_filename) :
    subfamily2family = dict()
    file = open(tsv_filename,'r')
    next(file)
    for line in file :
        line = line.rstrip()
        orf,subfamily,family = line.split('\t')
        subfamily2family[ subfamily ] = family
    file.close()
    return subfamily2family


def readingHhrFile(hhr_filename) :
    #print(hhr_filename)
    bestProb = 0
    bestHit = ''

    tag = 0
    file = open(hhr_filename,'r')
    for line in file :
        line = line.strip()
        if re.match(r'Query',line) :
            query = line.split()[1]
            continue
        
        if re.match(r'Match_columns',line) :
            qlen = float( line.split()[1] )
            continue

        if re.match(r'No Hit',line) :
            tag = 1
            continue
        
        if tag == 1 and line != '':
            liste = line.split()
            # print(liste)
            if not re.search(r'\(',liste[-2]) :
                target = liste[1]
                probs = float(liste[-9])
                qcoord = liste[-3]
                qstart = float(qcoord.split('-')[0])
                qend = float(qcoord.split('-')[1])            
                qcover = ( qend - qstart + 1.0 ) / qlen            
                scoord = liste[-2]
                sstart = float(scoord.split('-')[0])
                send = float(scoord.split('-')[1])            
                slen = float(liste[-1].replace('(','').replace(')',''))
                scover = ( send - sstart + 1 ) / slen
            else: # special case when liste[9] = 120-198(200)
                target = liste[1]
                probs = float(liste[-9])
                qcoord = liste[-3]
                qstart = float(qcoord.split('-')[0])
                qend = float(qcoord.split('-')[1])            
                qcover = ( qend - qstart + 1.0 ) / qlen            
                scoord,slen = liste[-2].split('(')
                sstart = float(scoord.split('-')[0])
                send = float(scoord.split('-')[1])            
                slen = float(slen.replace('(','').replace(')',''))
                scover = ( send - sstart + 1 ) / slen

            if probs >= bestProb :
                bestProb = probs
                bestHit = [query,target,str(probs),str(qcover),str(scover)]
    file.close()
    return bestHit


def runningHhblits(hhm_filename,hhblits_database,hhr_filename) :
    basename = os.path.basename(hhm_filename).split('.')[0]
    cmd = '/home/meheurap/programs/hhsuite-3.0-beta.3-Linux/bin/hhblits -i '+hhm_filename+' -o '+hhr_filename+' -d '+hhblits_database+'  -v 0 -p 50 -E 0.001 -z 1 -Z 32000 -B 0 -b 0 -n 2 -cpu 1'
    status = os.system(cmd)
    if status == 0 :
        return basename,True
    else:
        return basename,False




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run Pfam annotation using hhblits on the subfamilies')
    parser.add_argument('config_filename', help='the path of the CONFIG_FILENAME created by the subfamilies.py script')
    parser.add_argument('output_filename', help='the path of the OUTPUT_FILENAME')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs used by hhblits (default: 1)')

    args = parser.parse_args()
    
    if not os.path.exists(args.config_filename) :
        sys.exit(args.config_filename+' does not exist, exit')
    else:
        config_filename = os.path.abspath(args.config_filename)
        
    with open(config_filename) as f:
        data = json.load(f)

    cwd = data['directory']
    tsv_filename = os.path.abspath(cwd+'/'+'orf2subfamily2family.tsv')
    subfamily2family = subfamily2familyFunction(tsv_filename)
    output_filename = os.path.abspath(args.output_filename)
    
    hhm_directory = cwd+'/hhblits/hhm'
    hhblits_database = '/data7/proteinfams/HHpred/pfam'


    error = 0
    results = list()
    pool = ProcessPoolExecutor(args.cpu) # start 20 worker processes and 1 maxtasksperchild in order to release memory

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
            # if cpt == 20 :
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


    output = open(output_filename,'w')
    for root, dirs, files in os.walk('/data7/proteinfams/HHpred/results/'):
        for filename in files :
            hhr_filename = root+'/'+filename
            besthit = readingHhrFile(hhr_filename)
            subfamily = besthit[0]
            family = subfamily2family[subfamily]
            output.write( family+'\t'+'\t'.join(besthit)+'\n' )
    output.close()
