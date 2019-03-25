#! /usr/bin/env python

import os,sys,re
from collections import defaultdict


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
            annot = ';'.join(line.split(';')[2:]).strip()
            accession2annot[ accession ] = annot
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
            enzyme = enzyme.strip()
            annot = annot.split('[')[0]
            annot = annot.replace(enzyme,'')
            annot = annot.strip()        
            if accession not in accession2annot :
                accession2annot[ accession ] = annot
        else:
            continue
    file.close()    
    return accession2annot



if __name__ == "__main__":
    
    tcdb_db_filename = '/home/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb'
    tcdb_fasta_filename = '/data9/genasci/genasci_metabolism/analysis_after_filtering_drep/hmm/tcdb/tcdb.faa'
    query_fasta_filename = '/data7/proteinfams/3.6k.PF/mmseqsProteinClustering/3600genomes.4pub.cleaned.faa'
    blastp_filename = 'tcdb.blastout'
    tcdb_annot_filename = '/home/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb.dr'

    if not os.path.exists(blastp_filename) :
        cmd = 'blastp -query '+query_fasta_filename+' -db '+tcdb_db_filename+' -out '+blastp_filename+' -max_target_seqs 5000 -evalue 1e-20 -num_threads 6  -outfmt \"6 qseqid sseqid pident evalue bitscore qstart qend qlen sstart send slen\" '
        print(cmd)
        status = os.system(cmd)

    accession2annot = tcdbAnnotation(tcdb_annot_filename,tcdb_fasta_filename)
    orf2annot,orf2score = parsingBlastpOutput(blastp_filename)

    output = open(output_filename,'w')
    output.write('orf'+'\t'+'annot'+'\n')
    for orf,annot in orf2annot.items() :
        output.write(orf+'\t'+annot+'\n')
    output.close()

    sys.exit()



# makeblastdb -in /data9/genasci/genasci_metabolism/analysis_after_filtering_drep/hmm/tcdb/tcdb.faa -dbtype prot -out /home/meheurap/proteinCluster/natureCommRevision/mmseqsClustering/new_cov0.5_prob0.95/annotation/tcdb.faa
