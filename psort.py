#! /usr/bin/python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import multiprocessing as mp

''' this script parallelize the psort softare  '''

def runPsort(filename) :
    fasta_filename = "psort_tmp/fasta"+"/"+filename
    psortb_filename = "psort_tmp/output"+"/"+filename.replace(".fasta",".psortb")
    log_filename = "psort_tmp/log"+"/"+filename.replace(".fasta",".log")
    cmd = "psort --output long --positive "+fasta_filename+" >"+psortb_filename+" 2>"+log_filename        
    print(cmd)
    os.system(cmd)
    return "done"




def parsePsort(output_filename) :
    print(output_filename)
    tag_header = 0
    output = open(output_filename,'w')
    for root, dirs, files in os.walk("psort_tmp/output"):
        for filename in files :
            file = open(root+'/'+filename,'r')
            header = next(file).rstrip().split('\t')
            if tag_header == 0 :
                output.write('\t'.join(header)+'\n')
                tag_header = 1
                
            for line in file :
                line = line.rstrip()
                liste = line.split('\t')
                output.write('\t'.join(liste)+'\n')
            file.close()
    output.close()


fasta_filename = sys.argv[1]
cpu = int(sys.argv[2])
output_filename = sys.argv[3]




if os.path.exists('psort_tmp') :
    sys.exit('psort_tmp already exists, remove it first')
else :
    os.mkdir('psort_tmp')
    os.mkdir('psort_tmp/fasta')
    os.mkdir('psort_tmp/log')
    os.mkdir('psort_tmp/output')
    

nb = 898

###################
# splitting fasta #
###################

binSeq = ( nb/cpu ) + 1
cpt = 0
binNb = 1
output_fasta_filename = "psort_tmp/fasta"+"/"+str(binNb)+".fasta"
output = open(output_fasta_filename,"w")
for seq_record in SeqIO.parse(fasta_filename , "fasta") :
    cpt += 1
    if cpt < binSeq : 
        SeqIO.write(seq_record,output,"fasta")
    else :
        output.close()
        binNb += 1
        cpt = 1
        output_fasta_filename = "psort_tmp/fasta"+"/"+str(binNb)+".fasta"
        output = open(output_fasta_filename,"w")
        SeqIO.write(seq_record,output,"fasta")

output.close()


#################
# running psort #
#################


results = list()
pool = mp.Pool(processes=cpu,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory
for root, dirs, files in os.walk("psort_tmp/fasta"):
    for filename in files :
        results.append( pool.apply_async( runPsort, args= (filename,) ))
pool.close() # Prevents any more tasks from being submitted to the pool
pool.join() # Wait for the worker processes to exit


####################
# writting results #
####################


parsePsort(output_filename)
