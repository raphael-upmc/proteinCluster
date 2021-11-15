#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor,wait


def rpl2pfam() :
    rp2pfam = {'RPL14' : 'PF00238' , 'RPL15' : 'PF00828' , 'RPL16' : 'PF00252' , 'RPL18' : 'PF00861' , 'RPL22' : 'PF00237' , 'RPL24' : 'PF17136' , 'RPL2' : 'PF03947' , 'RPL3' : 'PF00297' , 'RPL4' : 'PF00573' , 'RPL5' : 'PF00673' , 'RPL6': 'PF00347' , 'RPS10' : 'PF00338' , 'RPS17' : 'PF00366' , 'RPS19' : 'PF00203' , 'RPS3' : 'PF00189' , 'RPS8' : 'PF00410' }

    pfam2rp = dict()
    for rp,pfam in rp2pfam.items() :
        pfam2rp[pfam] = rp

    return rp2pfam, pfam2rp

def buildingHmmDb(pfam2rp) :
    output_filename = 'rp.hmm'
    output = open(output_filename,'w')

    pfamList = list()
    name2accession = dict()
    pfam2desc = dict()
    pfamList_filename = '/env/cns/db/Pfam/Pfam_latest/Pfam-A.hmm'
    file = open(pfamList_filename,'r')
    for line in file :
        line = line.rstrip()
        pfamList.append(line)

        liste = line.split()
        if re.match(r'NAME',line) :
            name = ' '.join(liste[1:])
            accession = ''
            desc = ''

        if re.match(r'ACC',line) :
            accession = ' '.join(liste[1:])
            name2accession[ name ] = accession            
            
        if re.match(r'DESC',line) :
            desc = ' '.join(liste[1:])
            pfam2desc[accession] = desc+' ('+accession+')'

        if line == '//' :
            if accession.split('.')[0] in pfam2rp :
                print(pfam2desc[accession])
                output.write('\n'.join(pfamList)+'\n')
                pfamList = []
            else:
                pfamList = []
                continue
                print(pfam2desc[accession])
    file.close()
    output.close()

def HMMsearch() :
    hmm_database = ''
    directory = ''
    for directory in directoryList :
        for root, dirs, files in os.walk(directory):
            for filename in files :
                basename = filename.replace('.hmm','')
                hmm_filename = root+'/'+filename
                domtblout_filename = '/data7/proteinfams/genomicContext/domtblout/'+basename+'.domtblout'
                if not os.path.exists(domtblout_filename) :
                    cmd = 'hmmsearch -E 1e-10 --cpu 1 --domtblout '+domtblout_filename+' '+hmm_filename+' '+fasta_filename+' >/dev/null 2>/dev/null'
                    status = os.system(cmd)
                    print(str(status)+'\t'+cmd)

def readingHMM(domtblout_filename) :
    orf2hmm = defaultdict(list)
    file = open(domtblout_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search(r'^#',line) :
            continue

        liste = line.split()
        orf = liste[0]
        length = liste[2]
        hmm = liste[3]
        ko = liste[3].split(".")[0]
        hmmLength = liste[5]
        evalue = liste[6]
        bitscore = liste[7]
        cEvalue = liste[11] # conditional Evalue
        iEvalue = liste[12] # independant Evalue

        hmmS = liste[15]
        hmmE = liste[16]

        aliS = liste[17]
        aliE = liste[18]

        envS = liste[19]
        envE = liste[20]

        orfCover = float( int(envE) - int(envS) + 1 ) / float(length)
        hmmCover = float( int(hmmE) - int(hmmS) + 1 ) / float(hmmLength)
        orf2hmm[ orf ].append( hmm , float(evalue) , float(bitscore) , orfCover , hmmCover , length, hmmLength )
    file.close()
    return orf2hmm


if __name__ == "__main__":

    rp2pfam, pfam2rp = rpl2pfam()
    buildingHmmDb(pfam2rp)

    # parser = argparse.ArgumentParser(description='extracting the 16 ribosomal proteins')
    # parser.add_argument('protein_filename', help='the path of the FASTA_PROTEIN_FILE')
    # parser.add_argument('orf2bin_filename',help='the path of the ORF2BIN_FILE, this file is a tab-separated file, first col is the orf, second col is the genome. First line is skipped')
    # parser.add_argument('feature_filename',help='the path of the FEATURE_FILE, this file is a tab-separated file.')
    # parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    # parser.add_argument('output_summary_filename',help='the path of the OUTPUT_SUMMARY_FILE, results will be stored in this file')
    # parser.add_argument('--cpu',type=int,default=6,help='number of CPUs (default: 6)')
    # parser.add_argument('--rp14',action='store_true',default=False,help='only consider 14RPs for Archaea (default: False)')
    
    # args = parser.parse_args()
    
    # if os.path.exists(args.protein_filename) :
    #     protein_filename = os.path.abspath(args.protein_filename)
    # else:
    #     sys.exit(args.protein_filename+' does not exist, exit')

    # if os.path.exists(args.orf2bin_filename) :
    #     orf2bin_filename = os.path.abspath(args.orf2bin_filename)
    # else:
    #     sys.exit(args.orf2bin_filename+' does not exist, exit')

    # if os.path.exists(args.feature_filename) :
    #     feature_filename = os.path.abspath(args.feature_filename)
    # else:
    #     sys.exit(args.feature_filename+' does not exist, exit')


    # cpu = args.cpu
    
    # output_summary_filename = os.path.abspath(args.output_summary_filename)
    
    # output_aln_filename = os.path.abspath(args.output_filename)
    # cwd = '/'.join(output_aln_filename.split('/')[:-1])
    
    # print('protein_filename: '+protein_filename)
    # print('orf2bin_filename: '+orf2bin_filename)
    # print('feature_filename: '+feature_filename)
    # print('output_aln_filename: '+output_aln_filename)
    # print('output_summary_filename: '+output_summary_filename)
    # print('number of CPUs: '+str(cpu))
    # print('current working directory: '+cwd)
    # print('14RP: '+str(args.rp14))
    
    # folder = cwd+"/16RP_results"
    # if os.path.exists(folder) :
    #     sys.exit(folder+" already exists, remove it first")

    # os.mkdir(folder)



