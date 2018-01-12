#! /usr/bin/python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO



if os.path.exists('psort_tmp') :
    sys.exit('psort_tmp already exists, remove it first')
else :
    os.mkdir('psort_tmp')
    os.mkdir('psort_tmp/fasta')
    os.mkdir('psort_tmp/log')
    os.mkdir('psort_tmp/output')
    

nb = 150269
cpu = 8

binSeq = ( nb/cpu ) + 1
cpt = 0
binNb = 1
output_filename = "psort_tmp/fasta"+"/"+str(binNb)+".fasta"
output = open(output_filename,"w")
fasta_filename = "cprOnly.fa"
for seq_record in SeqIO.parse(fasta_filename , "fasta") :
    cpt += 1
    if cpt < binSeq : 
        SeqIO.write(seq_record,output,"fasta")
    else :
        output.close()
        binNb += 1
        cpt = 1
        output_filename = "psort_tmp/fasta"+"/"+str(binNb)+".fasta"
        output = open(output_filename,"w")
        SeqIO.write(seq_record,output,"fasta")

output.close()

for root, dirs, files in os.walk("psort_tmp/fasta"):
    for filename in files :
        fasta_filename = root+"/"+filename
        psortb_filename = "psort_tmp/output"+"/"+filename.replace(".fasta",".psortb")
        log_filename = "psort_tmp/log"+"/"+filename.replace(".fasta",".log")
        cmd = "nohup psort --output long --positive "+fasta_filename+" >"+psortb_filename+" 2>"+log_filename+" &"
        print(cmd)
#        os.system(cmd)
