#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO




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

    output_filename = os.path.abspath(args.output_filename)
    cwd = '/'.join(cazy_output_filename.split('/')[:-1])
    
    print('protein_filename: '+protein_filename)
    print('orf2bin_filename: '+orf2bin_filename)
    print('output_filename: '+cazy_output_filename)
    print('number of CPUs: '+str(cpu))
    print('current working directory: '+cwd)

    if os.path.exists(cwd+"/16RP_results") :
        sys.exit(cwd+"/Results already exists, remove it first")
    
    orf2bin = dict()
    file = open(orf2bin_filename,'r')
    header = next(file)
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

output = open("rp16_table.tsv",'w')
output.write('genome\tscaffold\tL15\tL18\tL6\tS8\tL5\tL24\tL14\tS17\tL16\tS3\tL22\tS19\tL2\tL4\tL3\tS10'+'\n')

cpt = 0
for genome,seqList in genome2seqList.items() :
    cpt += 1
    fasta_filename = 'tmp'+'.faa'
    SeqIO.write(seqList,fasta_filename,'fasta')
    cmd = "/home/meheurap/.pyenv/shims/rp16.py -f "+fasta_filename+" -d /home/cbrown/databases/rp16/Laura/ 1>tmp.tsv 2>>error.log"
    status = os.system(cmd)
    print(str(cpt)+'\t'+genome+'\t'+str(status))


    if status == 0 :
        file = open('tmp.tsv','r')
        next(file)
        for line in file :
            output.write(genome+'\t'+line)
        file.close()


    for usearch_filename in ['tmp-usearch_prot-rpL14_JGI_MDM.filtered.b6','tmp-usearch_prot-rpS8_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL16_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL4_JGI_MDM.filtered.b6','tmp-usearch_prot-rpS3_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL5_JGI_MDM.filtered.b6','tmp-usearch_prot-rpS19_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL24_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL22_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL18_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL15_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL3_JGI_MDM.filtered.b6','tmp-usearch_prot-rpS10_JGI_MDM.filtered.b6','tmp-usearch_prot-rpS17_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL6_JGI_MDM.filtered.b6','tmp-usearch_prot-rpL2_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL14_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpS8_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL16_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL4_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpS3_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL5_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpS19_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL24_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL22_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL18_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL15_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL3_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpS10_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpS17_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL6_JGI_MDM.filtered.b6','tmp-usearch_nucl-rpL2_JGI_MDM.filtered.b6'] :

        if os.path.exists(usearch_filename) :
            os.remove(usearch_filename)

    if status == 0 :
        os.remove(fasta_filename)

output.close()
