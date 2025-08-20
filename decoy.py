#!/usr/bin/env python3

import sys
from Bio import SeqIO
with open(sys.argv[1]) as fp, open(sys.argv[2], 'w') as wfp:
    for target in SeqIO.parse(fp, 'fasta'):
        wfp.write(">%s\n%s\n" % (target.description,target.seq))

        decoyseq = target.seq[::-1]  #reversing target protein sequences
        decoydescription = "decoy_"+target.description
        wfp.write(">%s\n%s\n" % (decoydescription,decoyseq))        
