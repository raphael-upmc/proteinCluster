#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import argparse


def parsingBlastpOutput(blastp_filename) :
    orf2score = defaultdict(float)
    orf2annot = dict()
    file = open(blastp_filename,'r')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        orf = liste[0]
        score = float( liste[4] )
        accession = liste[1].split('|')[2]
        if accession in accession2annot :
            if score > orf2score[orf] :
                orf2score[orf] = score
                orf2annot[orf] = accession2annot[ accession ] 
        else:
            print(accession+'\t'+liste[1])
    file.close()
    return orf2annot,orf2score


def tcdbAnnotation(tcdb_annot_filename,tcdb_fasta_filename) :
    accession2annot = dict()
    file = open(tcdb_annot_filename,'r')
    for line in file :
        line = line.rstrip()
        try :
            accession = line.split()[0]
            enzyme = line.split(';')[1].strip()
            enzyme = '.'.join( enzyme.strip().split('.')[:2])
            annot = ';'.join(line.split(';')[2:]).strip()
            accession2annot[ accession ] = enzyme+'|'+annot
        except :
            print(line)
    file.close()

    file = open(tcdb_fasta_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search('>',line) :
            liste = line.split('|')
            accession = liste[2].strip()
            annot = liste[-1]
            enzyme = liste[-2]        
            enzyme = annot.split()[0]
            enzyme = '.'.join( enzyme.strip().split('.')[:2])
            
            annot = annot.split('[')[0]
            annot = annot.replace(enzyme,'')
            annot = annot.strip()        
            if accession not in accession2annot :
                accession2annot[ accession ] = enzyme+'|'+annot
        else:
            continue
    file.close()    
    return accession2annot



if __name__ == "__main__":

    tcdb_db_filename = '/groups/banfield/users/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb'
    tcdb_fasta_filename = '/groups/banfield/users/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb.faa'
    tcdb_annot_filename = '/groups/banfield/users/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb.dr'


    parser = argparse.ArgumentParser(description='running  and parsing the TCDB')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILE')
    parser.add_argument('blastp_filename',help='the path of the BLASTP_FILE')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')

    args = parser.parse_args()

    if os.path.exists(args.fasta_filename) :
        query_fasta_filename = os.path.abspath(args.fasta_filename)
    else:
        sys.exit(args.fasta_filename+' does not exist, exit')

    # output_filename = 'tcdb.annot'
    # blastp_filename = '/home/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb.blastout'    
    # query_fasta_filename = '/data7/proteinfams/3.6k.PF/mmseqsProteinClustering/3600genomes.4pub.cleaned.faa'

    blastp_filename = args.blastp_filename
    output_filename = args.output_filename
    
    if not os.path.exists(tcdb_db_filename+'.pin') :
        cmd = 'makeblastdb -in '+tcdb_fasta_filename+' -dbtype prot -out '+tcdb_db_filename
        print(cmd)
        status = os.system(cmd)
    else:
        print(tcdb_db_filename+' already exists, no need to run makeblastdb')

        
    if not os.path.exists(blastp_filename) :
        cmd = 'blastp -query '+query_fasta_filename+' -db '+tcdb_db_filename+' -out '+blastp_filename+' -max_target_seqs 5000 -evalue 1e-20 -num_threads 6  -outfmt \"6 qseqid sseqid pident evalue bitscore qstart qend qlen sstart send slen\" '
        print(cmd)
        status = os.system(cmd)
    else:
        print(blastp_filename+' already exists, no need to run blastp')

    print('reading and storing tcdb annotation...')
    accession2annot = tcdbAnnotation(tcdb_annot_filename,tcdb_fasta_filename)

    print('reading and parsing '+blastp_filename+'...')
    orf2annot,orf2score = parsingBlastpOutput(blastp_filename)

    print('writting output '+output_filename+'...')
    output = open(output_filename,'w')
    output.write('orf'+'\t'+'annot'+'\n')
    for orf,annot in orf2annot.items() :
        output.write(orf+'\t'+annot+'\n')
    output.close()

    sys.exit()




    # makeblastdb -in /data9/genasci/genasci_metabolism/analysis_after_filtering_drep/hmm/tcdb/tcdb.faa -dbtype prot -out /home/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb.faa
