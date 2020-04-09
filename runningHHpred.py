#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,wait
import json
import argparse

''' the script runs the HHpred software on the pfam db and parses the results '''


class Pfam :
    """ a pfam class """
    def __init__(self,accession,name,clanAccession,clanName,description):
        self.accession = accession
        self.name = name
        self.clanAccession = clanAccession
        self.clanName = clanName
        self.description = description

        
def readingPfamDescription():
    pfamAccession2pfamObject = dict()
    info_filename = "/groups/banfield/projects/multienv/proteinfams/NCBI_balanced_dataset/ncbiGenomeDbComprehensive20171212/PFAM/Pfam-A.clans.tsv"
    file = open(info_filename,"r")
    for line in file :
        line = line.rstrip()
        accession,clanAccession,clanName,name,description = line.split("\t")
        pfamAccession2pfamObject[accession] = Pfam(accession,name,clanAccession,clanName,description)
    file.close()
    return pfamAccession2pfamObject


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


def readingHhrFile(hhr_filename,hhblits_database_basename,accession2desc) :
    bestProb = 0
    bestHit = ''

    tag = 0
    file = open(hhr_filename,'r')

    try :
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
                if hhblits_database_basename == 'pfam' :
                    target = accession2desc[liste[1]]+' ('+liste[1]+')'
                elif hhblits_database_basename == 'scop40' :
                    target = accession2desc[liste[1]]+' ('+liste[1]+')'
                elif hhblits_database_basename == 'arCOG' :
                    target = accession2desc[liste[1]]
                else:
                    target = liste[1]
                    
                if not re.search(r'-',liste[-1]) :
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
                    probs = float(liste[-8])
                    qcoord = liste[-2]
                    qstart = float(qcoord.split('-')[0])
                    qend = float(qcoord.split('-')[1])            
                    qcover = ( qend - qstart + 1.0 ) / qlen            
                    scoord,slen = liste[-1].split('(')
                    sstart = float(scoord.split('-')[0])
                    send = float(scoord.split('-')[1])            
                    slen = float(slen.replace('(','').replace(')',''))
                    scover = ( send - sstart + 1 ) / slen

                if probs >= bestProb :
                    bestProb = probs
                    bestHit = [query,target,str(probs),str(qcover),str(scover)]
        file.close()
    except :
        #print(liste)
        return bestHit,False



    return bestHit,True


def runningHhblits(hhm_filename,hhblits_database,hhr_filename) :
    basename = os.path.basename(hhm_filename).split('.')[0]
    cmd = '/groups/banfield/users/meheurap/programs/hhsuite-3.0-beta.3-Linux/bin/hhblits -i '+hhm_filename+' -o '+hhr_filename+' -d '+hhblits_database+'  -v 0 -p 50 -E 0.001 -z 1 -Z 32000 -B 0 -b 0 -n 2 -cpu 1'
    status = os.system(cmd)
    if status == 0 :
        return basename,True
    else:
        return basename,False




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run Pfam hhm annotations using hhblits on the subfamilies')
    parser.add_argument('config_filename', help='the path of the CONFIG_FILENAME created by the subfamilies.py script')
    parser.add_argument('hhm_database', help='the path of the HHM_DATABASE')
    parser.add_argument('output_filename', help='the path of the OUTPUT_FILENAME')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs used by hhblits (default: 1)')
    parser.add_argument('--hhr_directory',help='the output directory where the results will be store (default: ./hhr_results)')
    parser.add_argument('--probs',type=float,default=95,help='hhblits probability threshold (default: 95)')
    
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

    if not os.path.exists(args.hhm_database+'_a3m.ffdata') :
        sys.exit(args.hhm_database+' does not exist, exit')
    else:
        hhblits_database = os.path.abspath(args.hhm_database)
        hhblits_database_basename = os.path.basename(hhblits_database)
        
        
#    hhblits_database = '/groups/banfield/projects/multienv/proteinfams/HHpred_db/pfam'

    if args.hhr_directory == None :
        hhr_directory = cwd+'/'+'hhr_results'
    else:    
        hhr_directory = args.hhr_directory



    if args.probs < 0 or args.probs > 100 :
        sys.exit('probability threshold should be between 0 and 100! ('+str(args.probs)+')')
        
    print('running HHblits on '+hhblits_database+' '+hhblits_database_basename)
    error = 0
    results = list()
    pool = ProcessPoolExecutor(args.cpu) # start 20 worker processes and 1 maxtasksperchild in order to release memory

    cpt = 0
    for root, dirs, files in os.walk(hhm_directory):
        for filename in files :
            cpt += 1

            subfamily = filename.replace('.hhm','')

            if subfamily not in data['clusters'] :
                continue
            print(subfamily)            
            hhm_filename = root+'/'+filename
            hhm_basename = os.path.basename(hhm_filename).split('.')[0]

            
            hhr_filename = hhr_directory+'/'+hhm_basename+'__'+hhblits_database_basename+'.hhr'
            if os.path.exists(hhr_filename) :
                continue
            else:
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
    print('done')

    

    print('writting the final output...')

    # get accession description
    accession2desc = dict()
    if hhblits_database_basename == 'pfam' :
        desc_filename = hhblits_database + '_hhm.ffdata'        
        file = open(desc_filename,'r')
        for line in file :
            line = line.rstrip()
            if re.match('NAME',line) :
                liste = line.split(' ; ')
                accession = liste[0].split()[-1]
                desc = liste[-1]
                accession2desc[ accession ] = desc
            else:
                continue
        file.close()
    elif hhblits_database_basename == 'scop40' :
        desc_filename = hhblits_database + '_hhm.ffdata'        
        file = open(desc_filename,'r')
        for line in file :
            line = line.rstrip()
            if re.match('NAME ',line) :
                liste = line.split()
                #print(liste)
                accession = liste[1]
                desc = liste[2]
                accession2desc[ accession ] = desc
                #print(accession+'\t'+desc)
            else:
                continue
        file.close()
    elif hhblits_database_basename == 'arCOG' :
        desc_filename = '/groups/banfield/projects/multienv/proteinfams/HHpred_db/arCOG_2157_annotations.tsv'
        file = open(desc_filename,'r')
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
            #print(liste)
            accession = liste[1]
            letter = liste[2]
            if len(liste) == 3:
                desc = ''
                #print(liste)
            else:
                desc = liste[3]
            if desc == '' :
                accession2desc[ accession ] = 'no_description'+' ('+accession+')'+' ('+letter+')'
            else:
                accession2desc[ accession ] = desc+' ('+accession+')'+' ('+letter+')'
        file.close()
        
    output = open(output_filename,'w')
    for root, dirs, files in os.walk(hhr_directory):
        for filename in files :
            hhr_filename = root+'/'+filename
            besthit,result = readingHhrFile(hhr_filename,hhblits_database_basename,accession2desc)
            if len(besthit) == 0 :
                print(hhr_filename+' ==> error')
                continue
            subfamily = besthit[0]
            probs = besthit[2]
            if float(probs) < args.probs :
                continue
            family = subfamily2family[subfamily]
            output.write( family+'\t'+'\t'.join(besthit)+'\n' )
    output.close()
    print('done')
