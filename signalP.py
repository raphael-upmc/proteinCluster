#! /usr/bin/env python

""" a wrapper to run signalP on large protein dataset """

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp

def runSignalP(fasta_filename,output_filename,log_filename) :
    cmd = "/data7/proteinfams/SignalP/signalp-4.1/signalp -f short -t gram- "+fasta_filename+" > "+output_filename+" 2>>"+log_filename
    print(cmd)
    status = os.system(cmd)
    return status




initial_fasta_filename = sys.argv[1]
final_output_filename = sys.argv[2]
cpu = int(sys.argv[3])


if os.path.exists('signalP_tmp') :
    sys.exit('signalP_tmp already exists, remove it first')
else :
    os.mkdir('signalP_tmp')
    os.mkdir('signalP_tmp/fasta')
    os.mkdir('signalP_tmp/log')
    os.mkdir('signalP_tmp/output')
    

###################
# splitting fasta #
###################
    
# Signalp allow no more than 10000 sequences as inpout, so we have to split the dataset into fasta bins of <10,000 sequences !

cpt2defline = dict()
binSeq = 10000
cpt = nb = 0
binNb = 1
output_fasta_filename = "signalP_tmp/fasta"+"/"+str(binNb)+".faa"
output = open(output_fasta_filename,"w")
for seq_record in SeqIO.parse(initial_fasta_filename , "fasta") :
    cpt += 1
    nb += 1
    cpt2defline[ str(nb) ] = seq_record.id
    if cpt < binSeq : 
        SeqIO.write( SeqRecord(seq=seq_record.seq,id=str(nb),description="") ,output,"fasta")
    else :
        output.close()
        binNb += 1
        cpt = 1
        output_fasta_filename = "signalP_tmp/fasta"+"/"+str(binNb)+".faa"
        output = open(output_fasta_filename,"w")
        SeqIO.write( SeqRecord(seq=seq_record.seq,id=str(nb),description="")  ,output,"fasta")
output.close()


###################
# running signalP #
###################


results = list()
pool = mp.Pool(processes=cpu,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory
for root, dirs, files in os.walk("signalP_tmp/fasta"):
    for filename in files :
        fasta_filename = root+'/'+filename
#        print(fasta_filename)
        output_filename = "signalP_tmp/output"+'/'+filename.replace('.faa','.signalP')
        log_filename = "signalP_tmp/log"+'/'+filename.replace('.faa','.log')
        results.append( pool.apply_async( runSignalP, args= (fasta_filename,output_filename,log_filename) ))
pool.close() # Prevents any more tasks from being submitted to the pool
pool.join() # Wait for the worker processes to exit



####################
# writting signalP #
####################

output = open(final_output_filename , 'w' )

for root, dirs, files in os.walk("signalP_tmp/output"):
    for filename in files :
        file = open(root+'/'+filename,"r")
        for line in file :
            line = line.rstrip()
            if re.match("#",line) :
                output.write(line+"\n")
            else :

                liste = line.split()
                defline = cpt2defline[int(liste[0])]
                result = liste[9]
                output.write(defline+"\t"+result+"\n")
        file.close()        
output.close()
