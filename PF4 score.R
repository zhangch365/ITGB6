tcga_mrna_clin[1:3,1:35]
tcga_mrna<-tcga_mrna_clin[,-c(1:34)]
tcga_clin<-tcga_mrna_clin[,c(1:34)]
rownames(tcga_mrna)<-tcga_clin$sample_id
tcga_mrna[1:5,1:5]

load('/data/xlk/caihua/ITGB6/TCGA分析/TCGA_GTEx_pancancer_mrna_pheno.rdata')
tcga_gtex_mrna_pheno[1:3,1:35]
tcga_gtex_mrna<-tcga_gtex_mrna_pheno[,-c(1:4)]
tcga_gtex_clin<-tcga_gtex_mrna_pheno[,c(1:4)]

rownames(tcga_gtex_mrna)<-tcga_gtex_clin$sample_id
tcga_gtex_mrna[1:5,1:5]

table(tcga_clin$project)
table(tcga_gtex_clin$project)
table(tcga_gtex_clin$sample_type)


tcga_clin$tumor1<-ifelse(tcga_clin$project%in%c('ACC','DLBC','GBM','LAML','LGG','MESO','OV','TGCT','UCS','UVM'),'1','2')
tcga_clin1<-subset(tcga_clin,tumor1=='2')

library(genefilter)
library(GSVA)
library(Biobase)
library(tidyr)
library(tibble)
tcga_gtex_mrna[1:5,1:5]
gene_set <-read.csv('mmc3.csv')
bg_genes <- split(as.matrix(gene_set)[,1], gene_set[,2])
gsva_matrix<- gsva(as.matrix(tcga_gtex_mrna), bg_genes,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix<-as.data.frame(t(gsva_matrix))

#gsva_matrix<-gsva_matrix[tcga_gtex_clin$sample_id,]
identical(tcga_gtex_clin$sample_id,rownames(gsva_matrix))
gsva_matrix1<-gsva_matrix
gsva_matrix1$Cell_type<-tcga_gtex_clin$project

table(tcga_gtex_clin$sample_type)
gsva_matrix1$Group<-ifelse(tcga_gtex_clin$sample_type%in%c('TCGA_tumor'),'tumor','normal')
table(gsva_matrix1$Cell_type)
colnames(gsva_matrix1)[1]<-'score'
gsva_matrix1$pair<-ifelse(gsva_matrix1$Cell_type%in%c('MESO','UVM'),'1','2')
gsva_matrix2<-subset(gsva_matrix1,pair=='2')


cols = c("#00B2FF","orange")
p<-ggplot(gsva_matrix2,aes(Cell_type,score,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cancer Type", y = "Expression") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = cols)+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..))
p
ggsave("PF4 score pancancer+GTEX_2 .pdf", p, width = 12, height = 4.5,dpi = 600)








