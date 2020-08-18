#! /usr/bin/env python

import os,re,sys
from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter
from ete3 import Tree

print(re.escape("ApbE (fmnB)"))
print(re.escape("Fmn_Bind (pplA)"))

feature_filename = '/groups/banfield/projects/multienv/proteinfams/GTDB/gtdb_31k_genomes.feature'

pfamSet = {'PF04205' : 'Fmn_Bind (pplA)','PF02424':'ApbE (fmnB)', 'PF13486' : 'PceA (PF13486)' , 'PF01794' : 'Ferric_reduct (PF01794)' , 'PF10029' : 'DUF2271 (PF10029)' , 'PF12094' : 'DUF3570 (PF12094)', 'PF16357' : 'PepSY_TM_like_2 (PF16357)' , 'PF03929' : 'PepSY_TM (PF03929)' , 'PF12801' : 'Fer4_5 (PF12801)' , 'NQR2_RnfD_RnfE (PF03116) + NAD_binding_1 (PF00175)' : 'NQR2_RnfD_RnfE (PF03116) + NAD_binding_1 (PF00175)'}
keggSet = {'K00351':'K00351 (Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1])', 'K03616':'K03616 (Na+-translocating ferredoxin:NAD+ oxidoreductase subunit B [EC:7.2.1.2])','K00376':'K00376 (nitrous-oxide reductase [EC:1.7.2.4])','K19339':'K19339 (NosR/NirI family transcriptional regulator, nitrous oxide reductase regulator)','K03885' : 'K03885 (NADH dehydrogenase [EC:1.6.99.3])', 'K04084' : 'K04084 (thiol:disulfide interchange protein DsbD [EC:1.8.1.8])'}

system2annot = {

    'Nos' : [
        set(['K00376 (nitrous-oxide reductase [EC:1.7.2.4])','ApbE (fmnB)']) ,
        set(['K00376 (nitrous-oxide reductase [EC:1.7.2.4])','Fmn_Bind (pplA)']) ] ,

    'NQR2_RnfD_RnfE (PF03116) + NAD_binding_1 (PF00175)' : [
        set(['NQR2_RnfD_RnfE (PF03116) + NAD_binding_1 (PF00175)']) ] ,

    
    'Fer4_5' : [
        set(['Fer4_5 (PF12801)','ApbE (fmnB)']) ,
        set(['Fer4_5 (PF12801)','Fmn_Bind (pplA)']) ,
        set(['K19339 (NosR/NirI family transcriptional regulator, nitrous oxide reductase regulator)','ApbE (fmnB)']) ,
        set(['K19339 (NosR/NirI family transcriptional regulator, nitrous oxide reductase regulator)','Fmn_Bind (pplA)']) ],
        
    'Rnf' : [
        set(['Fmn_Bind (pplA)','K03616 (Na+-translocating ferredoxin:NAD+ oxidoreductase subunit B [EC:7.2.1.2])']) ,
        set(['ApbE (fmnB)','K03616 (Na+-translocating ferredoxin:NAD+ oxidoreductase subunit B [EC:7.2.1.2])']) ] ,
    
    'NQR' : [
        set(['K00351 (Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1])','Fmn_Bind (pplA)']) ,
        set(['K00351 (Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1])','ApbE (fmnB)']) ] ,
    
    'OrganohalideReductase' :  [
        set(['PceA (PF13486)','ApbE (fmnB)']) ,
        set(['PceA (PF13486)','Fmn_Bind (pplA)']) ] ,
    'EET' : [
        set(['K03885 (NADH dehydrogenase [EC:1.6.99.3])','Fmn_Bind (pplA)']) ,
        set(['K03885 (NADH dehydrogenase [EC:1.6.99.3])','ApbE (fmnB)']) ] ,

    'Ferric_reduct (PF01794)' : [
        set(['Ferric_reduct (PF01794)','Fmn_Bind (pplA)']) ,
        set(['Ferric_reduct (PF01794)','ApbE (fmnB)']) ,
        set(['Ferric_reduct (PF01794)','DUF2271 (PF10029)']) ,
        set(['Ferric_reduct (PF01794)','DUF3570 (PF12094)']) ] ,

    'dsbD' : [
        set(['K04084 (thiol:disulfide interchange protein DsbD [EC:1.8.1.8])','Fmn_Bind (pplA)']) ,
        set(['K04084 (thiol:disulfide interchange protein DsbD [EC:1.8.1.8])','ApbE (fmnB)']) ,
        set(['K04084 (thiol:disulfide interchange protein DsbD [EC:1.8.1.8])','DUF2271 (PF10029)']) ,
        set(['K04084 (thiol:disulfide interchange protein DsbD [EC:1.8.1.8])','DUF3570 (PF12094)']) ] ,
    
    'PepSY' : [
        set(['PepSY_TM (PF03929)','Fmn_Bind (pplA)']) ,
        set(['PepSY_TM (PF03929)','ApbE (fmnB)']) ,
        set(['PepSY_TM (PF03929)','DUF2271 (PF10029)']) ,
        set(['PepSY_TM (PF03929)','DUF3570 (PF12094)']) ,
        set(['PepSY_TM_like_2 (PF16357)','Fmn_Bind (pplA)']) ,
        set(['PepSY_TM_like_2 (PF16357)','ApbE (fmnB)']) ,
        set(['PepSY_TM_like_2 (PF16357)','DUF2271 (PF10029)']) ,
        set(['PepSY_TM_like_2 (PF16357)','DUF3570 (PF12094)']) ]

}




orf2desc = defaultdict(set) # cases of fusion or multidomain

###########
# SIGNALP #
###########
print('reading signalp....')
directory = '/groups/banfield/projects/multienv/proteinfams/GTDB/annotation/SignalP'
orf2signalp = dict()

for root, dirs, files in os.walk(directory):
    for filename in files :
        signalp_filename = root+'/'+filename
        file = open(signalp_filename,'r')
        for line in file :
            line = line.rstrip()
            if re.match('#',line) :
                continue
            liste = line.split('\t')
            orf = liste[0]
            signalp = liste[1]
            orf2signalp[orf] = signalp
        file.close()

#########
# TMHMM #
#########
print('reading tmhmm....')
directory = '/groups/banfield/projects/multienv/proteinfams/GTDB/annotation/TMHMM'
orf2tmhmm = dict()

for root, dirs, files in os.walk(directory):
    for filename in files :
        tmhmm_filename = root+'/'+filename
        file = open(tmhmm_filename,'r')
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
            orf = liste[0]
            tmhmm = int(liste[4].replace('PredHel=',''))
            orf2tmhmm[orf] = tmhmm
        file.close()

########
# PFAM #
########
print('reading pfam....')
gtdb_directory = '/groups/banfield/projects/multienv/proteinfams/GTDB/annotation/'

name2pfam = dict()
pfam2name = dict()
info_filename = '/groups/banfield/projects/multienv/proteinfams/GTDB/annotation/Pfam-A.clans.version33.0.tsv'
file = open(info_filename,"r")
for line in file :
    line = line.rstrip()
    accession,clanAccession,clanName,name,description = line.split("\t")
    pfam2name[accession] = name
    name2pfam[name] = accession
file.close()

orfSet_PF02424 = set()
orfSet_PF04205 = set()
orfSet_PF02683 = set()
cpt = 0
orf2pfam = defaultdict(list)
pfam_filename = gtdb_directory+'/'+'gtdb2Pfam-A-version33.0.result'
file = open(pfam_filename,"r")
for line in file :
    line = line.rstrip()
    if re.match(r'#',line) :
        continue
    liste = line.split()
    orfName = liste[0]
    pfamName = liste[1]
    pfamAccession = name2pfam[pfamName]
    if pfamAccession == 'PF02683' :
        orfSet_PF02683.add(orfName)

    if pfamAccession in pfamSet :
        orf2desc[orfName].add(pfamSet[pfamAccession])

    if pfamAccession == 'PF02424' :
        orfSet_PF02424.add(orfName)

    if pfamAccession == 'PF04205' :
        orfSet_PF04205.add(orfName)

    start,end = liste[3].split('-')
    cEvalue = liste[5] # conditional Evalue
    orf2pfam[orfName].append( ( int(start) , int(end) , pfamAccession , pfam2name[pfamAccession], cEvalue ) )
file.close()

orf2pfamArchitecture = dict()
for orf,liste in orf2pfam.items() :
    architecture = list()
    for domainInfoList in sorted( liste ,key=itemgetter(0,1) ) :
        architecture.append(domainInfoList[3]+" ("+domainInfoList[2]+")")

        #print(architecture)
    orf2pfamArchitecture[orf] = " + ".join(architecture)
    if orf2pfamArchitecture[orf] in pfamSet : # fusion
        orf2desc[orf].add(pfamSet[ orf2pfamArchitecture[orf] ])



########
# KEGG #
########
print('reading kegg....')
keggAccession2keggObject = dict()
filename = '/shared/db/kegg/kofam/latest/metadata/ko_list'
file = open(filename,"r")
for line in file :
    line = line.rstrip()
    liste = line.split("\t")
    accession = liste[0]
    description = liste[-1]
    keggAccession2keggObject[accession] = accession+' ('+description+') '
file.close()    

orf2kegg = dict()
kegg_filename = '/groups/banfield/projects/multienv/proteinfams/GTDB/annotation/'+'gtdb2KEGG-A-version33.0.result'
file = open(kegg_filename,'r')
next(file)
for line in file :
    line = line.rstrip()
    if re.match('#',line) :
        continue
    liste = line.split()
    orf = liste[0]
    accession = liste[1]
    if accession not in keggAccession2keggObject :
        orf2kegg[orf] = accession
    else:
        orf2kegg[orf] = keggAccession2keggObject[accession]

    if accession in keggSet :
        if accession == 'K03885' : # nad from EET ==> check if 1 transmembrane helix
            if orf2tmhmm[ orf ] > 0 :
                orf2desc[orf].add(keggSet[accession])
                #print(orf+'\t'+str(orf2tmhmm[orf]))
        else:
            orf2desc[orf].add(keggSet[accession])
            #print(orf+'\t'+keggSet[accession])
        
file.close()
print('kegg: '+str(len(orf2kegg)))

print('reading gtdb taxonomy....')
bin2taxonomy = dict()
filename = '/shared/db/gtdb/latest/taxonomy/gtdb_taxonomy.tsv'
file = open(filename,"r")
for line in file :
    line = line.rstrip()
    genomeAccession,taxonomy = line.split("\t")
    bin2taxonomy[ genomeAccession ] = taxonomy
file.close()

print('reading fasta file....')
annot2seqList =  defaultdict(list)
seq2len = dict()
fasta_filename = '/groups/banfield/projects/multienv/proteinfams/GTDB/gtdb_31k_genomes.faa'
for record in SeqIO.parse(fasta_filename,'fasta') :
    seq2len[record.id] = len(record)
    if record.id in orfSet_PF02424 :
        annot2seqList['PF02424'].append(record)
    if record.id in orfSet_PF04205 :
        annot2seqList['PF04205'].append(record)
    if record.id in orfSet_PF02683 :
        annot2seqList['PF02683'].append(record)

print(len(orfSet_PF02683))
SeqIO.write(annot2seqList['PF02424'],'PF02424.faa','fasta')
SeqIO.write(annot2seqList['PF04205'],'PF04205.faa','fasta')
SeqIO.write(annot2seqList['PF02683'],'PF02683.faa','fasta')

orf_filename = 'orf.txt'
output = open(orf_filename,'w')
for orf in orf2desc :
    output.write(orf+'\n')
output.close()





######################################
# running extracting genomic context #
######################################

print('performing genomic context extraction...')
genomicContext_filename = 'extracytoplasmic_flavinylation_system.genomicContext'
cmd = '/home/meheurap/scripts/proteinCluster/extractingGenomicContext.py '+orf_filename+' '+feature_filename+' '+genomicContext_filename+' -k 5'
status = os.system(cmd)
print('status: '+str(status))
print('done')

genome2scaffold2nb2annot = dict()
file = open(genomicContext_filename,'r')
next(file)
for line in file :
    line = line.rstrip()
    liste = line.split('\t')
    centroid = liste[0]
    genome = liste[1]
    orf = liste[-1]
    scaffold = liste[2]
    
    if genome not in genome2scaffold2nb2annot :
        genome2scaffold2nb2annot[genome] = dict()
        
    if scaffold not in genome2scaffold2nb2annot[genome] :
        genome2scaffold2nb2annot[genome][scaffold] = dict()

    coord = liste[3]+'..'+liste[4]+' ('+liste[5]+')'
    genome2scaffold2nb2annot[genome][scaffold][int(liste[3])] = orf+'\t'+scaffold+'\t'+coord    
file.close()


#############################
# detecting the RNF and NqR #
#############################

genome2systems = defaultdict(set)
orfFinalSet = set()
orf2system = defaultdict(set)
for genome,scaffold2nb2annot in genome2scaffold2nb2annot.items() :
    for scaffold, nb2annot in scaffold2nb2annot.items() :    
        nb2annotList = sorted(nb2annot)
        for i in range(len(nb2annotList)) :
            nb = nb2annotList[i]
            annot = nb2annot[nb]
            orf = annot.split('\t')[0]

            
            if orf in orf2pfamArchitecture :
                pfam = orf2pfamArchitecture[orf]
            else:
                pfam = 'Na'

            if orf in orf2kegg :
                kegg = orf2kegg[orf]
            else:
                kegg = 'Na'

            if orf in orf2signalp :
                signalp = orf2signalp[orf]
            else:
                signalp = 'Na'
                
            if orf in seq2len :
                l = str(seq2len[orf])
            else:
                l = 'Na'

            if genome not in bin2taxonomy :
                taxonomy = 'Na'
            else:
                taxonomy = bin2taxonomy[genome]

            if orf in orf2desc :
                desc = ','.join(list(orf2desc[orf]))
            else:
                desc = '-'

            #print(genome+'\t'+'taxonomy'+'\t'+annot+'\t'+str(l)+'\t'+signalp+'\t'+pfam+'\t'+kegg+'\t'+desc)

            if orf in orf2desc : # check 5 genes before and after
                if i - 5 < 0 :
                    start = 0
                else:
                    start = i-5

                if i + 5 < len(nb2annotList) :
                    end = i+5
                else:
                    end = len(nb2annotList) - 1

                j = start
                annotSet = set()
                orf2Set = set()
                while j <= end :
                    orf2 = nb2annot[ nb2annotList[j] ].split('\t')[0]
                    orf2Set.add(orf2)

                    if orf2 in orf2desc :
                        annotSet |= orf2desc[orf2]

                    j += 1

                #print(orf+'\t'+str(annotSet))
                for system,annotList in system2annot.items() :
                    for liste in annotList :
                        if liste.issubset(annotSet) and len(orf2desc[orf].intersection(liste)) > 0 : # check if current orf has the right annotation
                            orf2system[orf].add(system)
                            genome2systems[ genome ].add(system)
                            orfFinalSet |= orf2Set
                            #print(system+'\t'+str(liste))
                        else:
                            continue

                



######################
# writing the output #
######################


system2genomes = defaultdict(set)
output = open('extra-cytoplasmic_flavinylation_system.annot','w')
output2 = open('orf2annot.txt','w')
output2.write('orf\tannot\n')
scaffold2family = defaultdict(list)
for genome,scaffold2nb2annot in genome2scaffold2nb2annot.items() :
    genomeTag = 0
    for scaffold, nb2annot in scaffold2nb2annot.items() :
        scaffoldTag = 0
        orfTag = 0
        for nb,annot in sorted(nb2annot.items()) :

            orf = annot.split('\t')[0]
            if orf not in orfFinalSet :
                orfTag = 1
                continue

            scaffoldTag = 1
            genomeTag = 1
            
            if orf in orf2system :
                system =  ','.join(list(orf2system[orf]))
            else:
                system = '-'
            
            if orf in orf2desc :
                desc = ','.join( list(orf2desc[orf]) )
#                print(orf+'\t'+desc)
                if re.search(r'ApbE\ \(fmnB\)',desc) :
                    output2.write(orf+'\t'+system+'\n')
                if re.search(r'Fmn_Bind\ \(pplA\)',desc) :
                    output2.write(orf+'\t'+system+'\n')
            else:
                desc = '-'
            
            if orf in orf2pfamArchitecture :
                pfam = orf2pfamArchitecture[orf]
            else:
                pfam = 'Na'

            if orf in orf2kegg :
                kegg = orf2kegg[orf]
            else:
                kegg = 'Na'

            if orf in orf2signalp :
                signalp = orf2signalp[orf]
            else:
                signalp = 'Na'
                
            if orf in seq2len :
                l = str(seq2len[orf])
            else:
                l = 'Na'

            if genome not in bin2taxonomy :
                taxonomy = 'Na'
            else:
                taxonomy = bin2taxonomy[genome]

            if orf in orf2tmhmm :
                tmhmm = str(orf2tmhmm[orf])
            else:
                tmhmm = 'Na'
                                        
                
            if orfTag == 1 :
                output.write('\n')
                orfTag = 0
            output.write(genome+'\t'+taxonomy+'\t'+annot+'\t'+str(l)+'\t'+signalp+'\t'+tmhmm+'\t'+pfam+'\t'+kegg+'\t'+desc+'\t'+system+'\n')
            system2genomes[system].add(genome)
        if scaffoldTag == 1 :
            output.write('\n')
    if genomeTag == 1 :
        output.write('\n')
output.close()
output2.close()


# for system,liste in system2genomes.items() :
#     print(system+'\t'+str(len(liste)))





systemList = ['Rnf','NQR','Nos','EET','OrganohalideReductase', 'Ferric_reduct (PF01794)' , 'PepSY' , 'dsbD' , 'NQR2_RnfD_RnfE (PF03116) + NAD_binding_1 (PF00175)' , 'Fer4_5' ]

output = open('genome2systems.matrix','w')
output.write('\t'+'\t'.join(systemList)+'\n')
for asm,liste in genome2systems.items() :
    line = asm
    for system in systemList :
        if system in liste :
            line += '\t'+'1'
        else:
            line += '\t'+'0'
    output.write(line+'\n')
file.close()


otuSet = set()
tree_filename = '/groups/banfield/projects/multienv/proteinfams/GTDB/gtdb_genus_genomes.16RP.aln.treefile'
t = Tree(tree_filename)
for leaf in t:
    #print(leaf.name)
    otuSet.add(leaf.name)


output = open('genome2systems_genera.matrix','w')
output.write('\t'+'\t'.join(systemList)+'\n')
for asm,liste in genome2systems.items() :
    if asm not in otuSet :
        continue
    line = asm
    for system in systemList :
        if system in liste :
            line += '\t'+'1'
        else:
            line += '\t'+'0'
    output.write(line+'\n')
file.close()



