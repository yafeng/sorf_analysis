import sys

'''
sort the input file according to SpectEvalue in assending order first
sort -s -g -t$'\t' -k13,13 file
the script uses  T-TDC method in following paper to calculate FDR
Improved False Discovery Rate Estimation Procedure for Shotgun Proteomics.
Uri Keich, Attila Kertesz-Farkas,and William Stafford Noble. 2015 JPR

'''

input=open(sys.argv[1],'r')
decoy_prefix = sys.argv[2]
output=open(sys.argv[3],'w')

header=input.readline().strip().split("\t")
header=header+["PSM-FDR","Pep-FDR"]
output.write("\t".join(header)+"\n")

targetcount=0
decoycount=0

pep_dic={}

for line in input:
    row=line.strip().split('\t')
    pep=row[8]
    acc=row[9]

    if decoy_prefix in acc:
        decoycount+=1
    else:
        targetcount+=1

    FDR=float(decoycount)/float(targetcount)
 
    if pep not in pep_dic:
       pep_dic[pep] = FDR
    
    row.append(str(FDR)) # PSM-FDR
    row.append(str(pep_dic[pep])) #Pep-FDR

    if FDR<0.01:
        if decoy_prefix not in acc: #write only target PSMs
            output.write("\t".join(row)+"\n")
    else:
        break;

print ("target hits",targetcount)
print ("decoy hits", decoycount)

input.close()
output.close()
