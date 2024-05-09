setwd("~/Documents/smORF_project/PhyloCSF analysis")

library(readxl)
df.mouse = read_excel("mouse.combined.PhyloCSF.xlsx")
colnames(df.mouse)

df.mouse = df.mouse[df.mouse$NovelCodons!=0,]


df.mouse = df.mouse[df.mouse$RelBL_Unanno>=0.5,]

df.mouse = df.mouse[!is.na(df.mouse$PhyloCSFPsi_Unanno),]


df.human = read_excel("human.combined.PhyloCSF.xlsx")

df.human = df.human[df.human$NovelCodons!=0,]
df.human = df.human[df.human$RelBL_Unanno>=0.5,]

df.human = df.human[!is.na(df.human$PhyloCSFPsi_Unanno),]


df.cor = df.mouse[,c(4,6,58)]

df = df.cor

colnames(df) = c("CHR","BP","Prob.score")
# add additional column SNP, set it as row numbers
df$SNP  = rownames(df)
# reorder the data frame
df = df[,c(4,1,2,3)]
df$CHR = gsub("chrX","23",df$CHR)
df$CHR = gsub("chrY","24",df$CHR)
df$CHR = gsub("chrM","25",df$CHR)
df$CHR = gsub("chr","",df$CHR)
df$CHR = as.numeric(df$CHR)

library(qqman)
library(RColorBrewer)
pal <- brewer.pal(7,"Dark2")

pdf("mouse_smORF_phyloCSF_psi.score.pdf",width=5.5,height=4,useDingbats = FALSE)
manhattan(df,col=pal,
          suggestiveline = F,
          genomewideline = F,
          p = "Prob.score",logp = F,
          ylab="PhyloCSF_psi(decibans)",ylim=c(min(df$Prob.score),max(df$Prob.score)+10))
dev.off()

setwd("~/Documents/smORF_project/nextflow_work/smORF_table/")

df.psm = read.table("Mouse_tissue_all_novelsmORF.filter_psmtable.blastp.parsed.txt",
                    sep = "\t",stringsAsFactors = F, comment.char = "",quote = "",header = T)
df.psm.filter = df.psm[df.psm$blastp_category=="novelpep (map to known protein with more than 2 mismatched aa)" |
                      df.psm$blastp_category=="novelpep (no match to known protein found)",]
df.psm2 = read.table("Mouse_secretome_novelsmORF.filter_psmtable.blastp.parsed.txt",
                     sep = "\t",stringsAsFactors = F, comment.char = "",quote = "",header = T)
df.psm2.filter = df.psm2[df.psm2$blastp_category=="novelpep (map to known protein with more than 2 mismatched aa)" |
                         df.psm2$blastp_category=="novelpep (no match to known protein found)",]

knownpeps = read.table("/Users/yafengzhu/Documents/smORF_project/smORF_DB_mouse/Mouse_uniprot+gencode.v23.trypticpep.0-1miss.txt",stringsAsFactors = F)

df.psm.filter.merge = rbind(df.psm.filter[,c("Protein","SpecEValue")],df.psm2.filter[,c("Protein","SpecEValue")])

df.psm.filter = df.psm.filter.merge
df.psm.filter$clean_pep = apply(df.psm.filter,1,function(x) substr(x[10],3,nchar(x[10])-2))
df.psm.filter$Sequence=gsub("[^A-Z]","",df.psm.filter$clean_pep)
df.psm.filter = df.psm.filter[!df.psm.filter$Sequence %in% as.character(knownpeps$V1),]

library(plyr)
df.score = ddply(df.psm.filter, "Protein", summarise, best_score = min(SpecEValue),psm_count = length(Protein))
df.score = df.score[df.score$Protein%in%df.mouse$Protein,]
rownames(df.score) = df.score$Protein

df.cor = unique(df.mouse[,c("Protein","Chromosome","RegionStart","PhyloCSFPsi_Unanno")])
df.cor$score = -log10(df.score[df.cor$Protein,]$best_score)
df.cor$psm_count = df.score[df.cor$Protein,]$psm_count
df.cor$log2_count =  log2(df.cor$psm_count)
df.cor = na.omit(df.cor)
  
pdf("mouse.smORF.FACS.plot.pdf",width = 7, height = 5,useDingbats = F)
ggplot(data=df.cor, x=score, y=PhyloCSFPsi_Unanno )+ geom_point(data=df.cor, aes(x=score, y=PhyloCSFPsi_Unanno, fill=log2_count),color="black",pch=21) + 
  scale_fill_gradient(limits=c(2, 12), low= "yellow", high="red") + 
  labs(colour="log2(PSM count)", x="-log10(SpecEval)", y="PhyloCSFPsi(novel codons)")
dev.off()

setwd("~/Documents/smORF_project/nextflow_work_human/")
df1 = read.table("PXD000561_novpeps.blastp.filter.2mismatch.txt",header = T, sep = "\t",
                 comment.char = "",quote = "")
df2 = read.table("PXD000865_novpeps.blastp.filter.2mismatch.txt",header = T, sep = "\t",
                 comment.char = "",quote = "")
df3 = read.table("PXD010154_novpeps.blastp.filter.2mismatch.txt",header = T,sep = "\t",
                 comment.char = "",quote = "")
df.human.psm = rbind(df1,df2,df3)
unique(df.human.psm$blastp_category)

knownpeps2 = read.table("/Users/yafengzhu/Documents/smORF_project/nextflow_work_human/Human_uniprot+gencode.v32.trypticpep.0-1miss.txt",stringsAsFactors = F)
df.human.psm$clean_pep = apply(df.human.psm,1,function(x) substr(x[10],3,nchar(x[10])-2))
df.human.psm$Sequence=gsub("[^A-Z]","",df.human.psm$clean_pep)
df.human.psm = df.human.psm[!df.human.psm$Sequence %in% as.character(knownpeps2$V2),]

setwd("~/Documents/smORF_project/PhyloCSF analysis")
df.human = read_excel("human.combined.PhyloCSF.xlsx")
df.human = df.human[df.human$NovelCodons!=0,]

df.score = ddply(df.psm.filter, "Protein", summarise, best_score = min(SpecEValue),psm_count = length(Protein))
rownames(df.score) = df.score$Protein

df.cor = unique(df.mouse[,c("Protein","Chromosome","RegionStart","PhyloCSFPsi_Unanno")])
df.cor$score = -log10(df.score[df.cor$Protein,]$best_score)
df.cor$psm_count = df.score[df.cor$Protein,]$psm_count
df.cor$log2_count =  log2(df.cor$psm_count)
df.cor = na.omit(df.cor)

pdf("mouse.smORF.FACS.plot.pdf",width = 7, height = 5,useDingbats = F)
ggplot(data=df.cor, x=score, y=PhyloCSFPsi_Unanno )+ geom_point(data=df.cor, aes(x=score, y=PhyloCSFPsi_Unanno, fill=log2_count),color="black",pch=21) + 
  scale_fill_gradient(limits=c(2, 8), low= "yellow", high="red") + 
  labs(colour="log2(PSM count)", x="-log10(SpecEval)", y="PhyloCSFPsi(novel codons)")
dev.off()

write.table(df.cor,"mouse.smoRFs.PhyloCSFPsi+MSscore.txt",sep = "\t",row.names = F,quote = F)