This is a tutorial for identificaiton of small open reading frame from MS-based proteomics data. 

Required softwares:

1.[MSGFPlus](https://github.com/MSGFPlus/msgfplus/release)

2.[ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)

3.[Mono](https://www.mono-project.com/download/stable/#download-lin)

4.[Nextflow](https://nextflow.io)


# Step 1 - Prepare MS data
**1.1 download MS data from PRIDE proteomics data repository**,here we use the mouse brain proteommics dataset PXD001250 as an example.

```wget  ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/10/PXD001250```

**1.2 convert MS raw data to mzml format**

```cat rawfilelist.txt |parallel -j 24 mono ThermoRawFileParser.exe -i={} -o=/data/home/yz332/MS_data/PXD001250_mzML/ -f=1 -m=1' ```
`-o` :the converted  output file path
`-j` : the numbers of cores to run in parallel

# Step 2 - Prepare small ORFs protein database
**2.1 download sORFs protein sequences** from [sORF.org](http://www.sorfs.org), which is a public repository of small open reading frames (sORFs) identified by ribosome profiling (RIBO-seq). After you get the file "sORFs.org.db.mouse.txt", you can proceed the netxt step. 

**2.2 download the Uniprot mouse reference protein database and supplment mouse sORFs proteins**

```python tofasta.py sORFs.org.db.mouse.txt sORFs.org.db.mouse.fasta

   cat sORFs.org.db.mouse.fasta uniprot.mouse.protein.fasta > uniprot.mouse.protein+sORFs.fasta
```
**2.3 create target and decoy protein database**,the decoy DB was produced by reversing protein sequences in the target database.

```python decoy.py uniprot.mouse.protein+sORFs.fasta uniprot.mouse.protein+sORFs.revCat.fasta```

**2.4 Index the target and decoy database**

```java -Xmx10000M -cp MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d uniprot.mouse.protein+sORFs.revCat.fasta -tda 0```


# Step 3 - Paralle database search using nextflow
**3.1 prepare a input, tab delimited file with two column**:  absolute file path of mzmL files and setname. The setnames are used to group mzmL. If you wish to get one FDR for MS files, use one setname for all files. If you wish to calculate separate FDR for a subset of MS files, use different setnames for different subgroups.

e.g. PXD001250_mzmlfiles_setnames.txt

/file_absolute_path/xxx.mzml	set1

/file_absolute_path/xxx.mzml	set1

….

**3.2 run nextflow to start database search**

```bash
#!/bin/bash


nextflow run main.nf --msgfjar MSGFPlus.jar \
--tdb uniprot.mouse.protein+sORFs.revCat.fasta  \ ## target decoy combined databases
--inst 3 \ #  0: low res LTQ, 1: Orbitrap/FTICR/Lumos, 2: TOF , 3: QE
--mods labelfree_Mods.txt \ #Modification file for MSGF+.
--mzmldef PXD001250_mzmlfiles_setnames.txt \# a tab deliminated text file with mzmlfilepath(absolute path) and setname 
--outdir results \ #Output directory
--activationFilter HCD \
--MS2error 0.02 \
--qval 0.01 \
-profile testing --tasks 4 --thread 1\
-resume \ # use it to resume the jobs from the last stopped process.
--PrecursorMassTolerance 10ppm \
--FragmentMethodID 3 \
--mzid2tsvConverter MzidToTsvConverter.exe
```
# Step 4 - Calculate class specific FDR for sORFs
```
sort -s -g -t$'\t' -k13,13 psmtable.txt > psmtable.sorted.txt
Python subgroupPepFDR.py --input psmtable.sorted.txt --output output_filename --knownproteins uniprot.fa --decoy_prefix XXX_  --psm_qval 0.01

```

# Step 5 - Blastp analysis of smORFs sequence

This step is to remomve smORF peptides that match to the sequences of known proteins. This is typically common even though smORF encoded peptides are included in known proteins, the search engine still assign smORF protein accession to these peptides due to known proteins don't have tryptic sites at the N or C-term of peptides. 

```
Python txt2fasta.py smORF_psmtable.qval0.01.txt SEP.fa

makeblastdb -in uniprot.mouse.protein.fa -dbtype prot

blastp -db uniprot.mouse.protein.fa -query SEP.fa -outfmt '6 qseqid sseqid pident qlen slen qstart qend sstart send mismatch positive gapopen gaps qseq sseq evalue bitscore' -num_threads 4 -max_target_seqs 1 -evalue 1000 -out blastp_out.txt

python parse_blastp_output.py --input smORF_psmtable.qval0.01.txt --blastp_result blastp_out.txt --fasta uniprot.mouse.protein.fa --output smORF_psmtable.qval0.01.blastp.txt
```


