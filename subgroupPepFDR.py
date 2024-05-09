import sys
import re
from itertools import product
from Bio import SeqIO
import numpy as np
import os
import getopt

'''
sort the file according to SpecEvalue first
sort -s -g -t$'\t' -k13,13 file

the script separate novel target and decoy hits, then uses T-TDC method in following paper
to calculate subgroup FDR.
Improved False Discovery Rate Estimation Procedure for Shotgun Proteomics.
Uri Keich, Attila Kertesz-Farkas,and William Stafford Noble. 2015 JPR

'''

def combination_replace(pep, from_aa, to_aa):
    options = [(c,) if c != from_aa else (from_aa, to_aa) for c in pep]
    return list(''.join(o) for o in product(*options))

def check_deamination(pep,dic): #deamination: N (Asn)>- D (Asp)
    pep_list=combination_replace(pep,"D","N")
    k=0
    for pep in pep_list:
        if pep in dic:
            k=1
            break;
    if k==1:
        return True # it is a known peptide when deamination considered
    else:
        return False #it is a novel peptide

def check_match_knownProt(s,protdic):
    check=0
    for seq in protdic: ## double check if it is known peptide
        if s in seq:
            check=1
            break;

    return check

decoy_prefix = "XXX_"
psm_qval = 0.05

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command!")
    print("Example: python subgroupPepFDR.py --input PSM_filename --output output_filename --knownproteins uniprot.fa --decoy_prefix XXX_  --psm_qval 0.01")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input_psm=','output=','splitchar='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--knownproteins': fa_file=arg
        elif opt == '--decoy_prefix': decoy_prefix=arg
        elif opt == '--psm_qval': psm_qval = float(arg)
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()


input=open(input_file,'r')# the input tsv file need to be sorted by SpectEval in assending order first
handle=SeqIO.parse(fa_file,'fasta')

knownPepDic={}

for record in handle:
    protseq = record.seq
    
    pep=row[1].replace("L","I")
    knownPepDic[pep]=row[0]

handle.close()

knownPepDic ={}
knownProtDic = {}
for record in handle2:
    seq=str(record.seq).replace("L","I")
    if seq not in known_seq:
        knownPepDic[seq]=record.id

output=open(output_file,'w')

header= ["NovelPSM-FDR","NovelPep-FDR"]

output.write("\t".join(header)+"\n")

novel_targetcount=0
novel_decoycount=0

novel_prefix = "SEP"
decoy_prefix = "XXX_"

novel_pep={}
novel_decoy_score=[]
for line in input:
    row=line.strip().split('\t')
    pep=re.sub("[\W\d]","",row[8].strip())#strip modifications
    pep2=pep.replace("L","I")
    pro=row[9]

    if novel_prefix not in pro:
        continue;
    elif pep2 in knownPepDic:
        continue;
    elif check_match_knownProt(pep2,):
        continue;
    elif check_deamination(pep2,knownPepDic):
        continue;

    else: #do novel only target-decoy FDR calculation
        if decoy_prefix in acc:
            decoycount+=1
        else:
            targetcount+=1

        FDR=float(decoycount)/float(targetcount)
 
        if pep not in pep_dic:
            novdl_pep[pep] = FDR

        if FDR<float(sys.arg[2]):
            if "REV_" not in pro: #write only target PSMs
                output.write("\t".join(row)+"\n")
        else:
            break;

print "Hits in novel search space: targe,decoy",novel_targetcount,novel_decoycount
print "max novel decoy score is",max(novel_decoy_score)
input.close()
output.close()
