import sys
import re

input=open(sys.argv[1],"r")
output=open(sys.argv[2],"w")

pep_dic={}

for line in input:
    row=line.split("\t")
    seq=re.sub("[^A-Z]","",row[9][2:-2])
    if seq=="Peptide":
        continue;
    if seq=="":
        continue;
    if seq not in pep_dic:
        output.write(">%s\n%s\n" % (seq,seq))
        pep_dic[seq]=1

print(len(pep_dic))

input.close()
output.close()
