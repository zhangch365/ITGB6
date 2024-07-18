setwd('/data/caihua2/ITGB6/上皮分析')
kang_epi<-readRDS('/data/caihua2/ITGB6/kang_epi.rds')
devtools::install_github("navinlabcode/copykat")
library(copykat)
library(Seurat)
# 读取数据
BiocManager::install("homologene")
library(homologene)
getwd()

# genelist=rownames(counts)
# genelist<-as.data.frame(genelist)
# genelist<-homologene(genelist$genelist, inTax = 10090, outTax = 9606)
# genelist=genelist[!duplicated(genelist[,2]),]
# 
# counts<-counts[genelist$`10090`,]
# rownames(counts)<-genelist$`9606`
# counts[1:4,1:4]

counts = as.matrix(ITGB6_4sample_tumor@assays$RNA@counts)
sc_cnv = copykat(rawmat = counts,
                 ngene.chr = 5,
                 sam.name = 'HNSC',
                 n.cores = 32,
                 genome = 'mm10')#小鼠是mm10,人是默认
pred<-as.data.frame(sc_cnv$prediction)









