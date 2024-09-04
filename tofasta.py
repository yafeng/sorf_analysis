import sys

input=open(sys.argv[1],"r")
output=open(sys.argv[2],"w")

print input.readline().split("\t")
pep_dic={}

for line in input:
    row=line.split("\t")
    header = "smORF_sORFs.org_" + row[0]
    seq=row[13].replace("*","")

    if seq not in pep_dic:
        output.write(">%s\n%s\n" % (header,seq))
        pep_dic[seq]=1

input.close()
output.close()
