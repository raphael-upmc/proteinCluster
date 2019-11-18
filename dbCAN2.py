#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
#from annotation import *
from Bio import SeqIO
import argparse

def mergingResults(hmm,hotpep,diamond) :
    result = defaultdict(int)
    if len(hmm) != 0 :
        hmm_res = list()
        for elt in sorted(hmm,key=lambda x:x[1]) :
            hmm_res.append(elt[0])
        result[ '+'.join(hmm_res) ] += 1

    if len(hotpep) != 0 :
        hotpep_res = list()
        for elt in hotpep :
            hotpep_res.append(elt)
        result[ '+'.join(hotpep_res) ] += 1


    if len(diamond) != 0 :
        diamond_res = list()
        for elt in sorted(diamond,key=lambda x:x[1]) :
            diamond_res.append(elt[0])
        result[ '+'.join(diamond_res) ] += 1



    if len(result) == 3 : # return hmm result
        return '+'.join(hmm_res)
    elif len(result) == 2 : # return choose the best or if equal hmm > diamond > hotpep
        best,nb = sorted(result.items(),key=lambda x:x[1],reverse=True)[0]
        if nb == 2 :
            return best
        else:
            if len(hmm) != 0 :
                return '+'.join(hmm_res)
            else:
                return '+'.join(diamond_res)
    elif len(result) == 1 : # return the only result
        return list(result.keys())[0]
    else :
        print(result)
        return 'error'

    
def parsingDiamond(filename) :

    orf2cazy = defaultdict(list)
    file = open(filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        orf = liste[0]
        subject = liste[1]
        cazy = subject.split('|')[1]
        start = int(liste[6])
        end = int(liste[7])
        orf2cazy[orf].append( (cazy,start,end) )
    file.close()
    return orf2cazy
        
def parsingHmm(filename) :
    orf2carbohydrate = defaultdict(list)
    file = open(filename,'r')
    header =next(file).rstrip()
#    print(header.split('\t'))
    
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
#        print(liste)
        hmm = liste[0].replace('.hmm','')
        orf = liste[2]
        start = int(liste[7])
        end = int(liste[8])
        orf2carbohydrate[ orf ].append( (hmm,start,end) )
    file.close()

    for orf,liste in orf2carbohydrate.items() :
        newList = sorted(liste,key= lambda x:x[1] )
        orf2carbohydrate[ orf ] = newList
#        print( orf+'\t'+str( newList ) )
    return orf2carbohydrate

def parsingHotpep(filename) :
    orf2carbohydrate = defaultdict(list)
    file = open(filename,'r')
    header =next(file).rstrip()
    
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        cazy = liste[0]
        orf = liste[2]
        orf2carbohydrate[ orf ].append( cazy )
    file.close()
    return orf2carbohydrate


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='annotating protein sequences with CAZY using dbCAN2')
    parser.add_argument('protein_filename', help='the path of the FASTA_PROTEIN_FILE')
    parser.add_argument('orf2bin_filename',help='the path of the ORF2BIN_FILE, this file is a tab-separated file, first col is the orf, second col is the genome. First line is skipped')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    parser.add_argument('--cpu',type=int,default=6,help='number of CPUs (default: 6)')

    args = parser.parse_args()
    
    if os.path.exists(args.protein_filename) :
        protein_filename = os.path.abspath(args.protein_filename)
    else:
        sys.exit(args.protein_filename+' does not exist, exit')

    if os.path.exists(args.orf2bin_filename) :
        orf2bin_filename = os.path.abspath(args.orf2bin_filename)
    else:
        sys.exit(args.orf2bin_filename+' does not exist, exit')


    cpu = args.cpu

    cazy_output_filename = os.path.abspath(args.output_filename)
    cwd = '/'.join(cazy_output_filename.split('/')[:-1])
    
    print('protein_filename: '+protein_filename)
    print('orf2bin_filename: '+orf2bin_filename)
    print('output_filename: '+cazy_output_filename)
    print('number of CPUs: '+str(cpu))
    print('current working directory: '+cwd)

    if os.path.exists(cwd+"/Results") :
        sys.exit(cwd+"/Results already exists, remove it first")
    
    orf2bin = dict()
    file = open(orf2bin_filename,'r')
    for line in file :
        line = line.rstrip()
        orf,genome = line.split('\t')
        orf2bin[ orf ] = genome
    file.close()

    ## ! when no results in hotpep, no file is created....
    genome2seqList = defaultdict(list)
    for record in SeqIO.parse(protein_filename,'fasta') :
        genome = orf2bin [ record.id ]
        genome2seqList[ genome ].append( record )
    print(len(genome2seqList))

    ##################
    # running dbCAN2 #
    ##################
    
    index2genome = dict()
    cpt = 0
    for genome, seqList in genome2seqList.items() :
        cpt += 1
        output_filename = cwd+'/'+'tmp'+'.faa'
        SeqIO.write(seqList,output_filename,'fasta')
        print(str(cpt)+'\t'+genome)
        index2genome[ str(cpt) ] = genome
        cmd = "/data7/proteinfams/dbCAN2/Tools/run_dbcan.py "+output_filename+" protein --dia_eval 1e-20 --dia_cpu "+str(cpu)+" --hmm_eval 1e-20 --hmm_cpu "+str(cpu)+" --hotpep_cpu "+str(cpu)+" --out_dir "+cwd+"/Results --db_dir /data7/proteinfams/dbCAN2/Tools/example/db --out_pre "+str(cpt)+'.'+' >/dev/null'
        print(cmd)
        status = os.system(cmd)
        print('done with status: '+str(status))

        
    ##########################
    # parsing dbCan2 results #
    ##########################
    
    genome2filename = defaultdict(list)
    for root, dirs, files in os.walk(cwd+'/Results/'):
        for filename in files :
            index = filename.replace('.uniInput','')
            index = index.replace('.diamond.out','')
            index = index.replace('.hmmer.out','')
            index = index.replace('.signalp.out','')
            index = index.replace('.Hotpep.out','')
            index = index.replace('.signalp.neg','')
            index = index.replace('.signalp.pos','')

            genome = index2genome[index]
            
            if filename == index+'.uniInput' :
                continue
            elif filename == index+'.diamond.out' :
                genome2filename[ genome ].append(root+'/'+filename)
            elif filename == index+'.hmmer.out' :
                genome2filename[ genome ].append(root+'/'+filename)
                continue
            elif filename == index+'.signalp.out' :
                continue
            elif filename == index+'.signalp.neg' :
                continue
            elif filename == index+'.signalp.pos' :
                continue
            elif filename == index+'.Hotpep.out' :
                genome2filename[ genome ].append(root+'/'+filename)
                continue
            else:
                print(genome+'\t'+filename)
                sys.exit('ERROR '+filename)

    output = open(cazy_output_filename,'w')
    for genome,files in genome2filename.items() :
        print(genome+'\t'+str(len(files)))
        if len(files) < 3 :
            print(files)
            sys.exit('error')
        for filename in files :
            if re.search(r'.diamond.out$',filename) :
                orf2diamond = parsingDiamond(filename)
            elif re.search(r'.hmmer.out$',filename) :
                orf2hmm = parsingHmm(filename)
            elif re.search(r'.Hotpep.out$',filename) :
                orf2hotpep = parsingHotpep(filename)
            else:
                sys.exit('ERROR '+filename)


        orfSet = set()
        for orf in orf2hotpep :
            orfSet.add(orf)    
        for orf in orf2hmm :
            orfSet.add(orf)
        for orf in orf2diamond :
            orfSet.add(orf)

        for orf in orfSet :
            #    print(orf+'\t'+str(orf2hotpep[orf])+'\t'+str(sorted(orf2hmm[ orf ],key=lambda x:x[1]))+'\t'+str(sorted(orf2diamond[ orf ],key=lambda x:x[1])))
            result = mergingResults(orf2hmm[orf],orf2hotpep[orf],orf2diamond[orf])
            if result == 'error':
                sys.exit('error')
            output.write(orf+'\t'+result+'\n')
    output.close()

    
