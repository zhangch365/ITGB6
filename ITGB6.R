library("harmony")
ITGB6 <- NormalizeData(ITGB6, normalization.method = "LogNormalize", scale.factor = 10000)
ITGB6 <- FindVariableFeatures(ITGB6, selection.method = "vst", nfeatures = 2000)
ITGB6<-ScaleData(ITGB6)%>%RunPCA(verbose=FALSE)
system.time({ITGB6 <- RunHarmony(ITGB6, group.by.vars = "orig.ident")})


ITGB6<-FindNeighbors(ITGB6, reduction = "harmony",dims = 1:30)
ITGB6<-FindClusters(ITGB6, reduction = "harmony",resolution = 0.5)
ITGB6<-RunUMAP(ITGB6,dims = 1:30,reduction = "harmony")
#ITGB6<-RunTSNE(ITGB6,dims = 1:30,reduction = "harmony")



ITGB6$keep <- ifelse(ITGB6$seurat_clusters %in% c(19),"not","keep")
ITGB6<- subset(ITGB6, subset = keep == "keep")


Idents(ITGB6)<-ITGB6$new.lables
DimPlot(ITGB6, reduction = "umap", pt.size = 0.5,label = T)

Idents(ITGB6)<-ITGB6$new.lables
ITGB6.markers <- FindAllMarkers(ITGB6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ITGB6.markers,'ITGB6.markers.csv')


################分群#########
ITGB6@meta.data$new.lables<-"Fibroblast.cells"
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(11))]<-"CD8.T.cells"#
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(4,12))]<-"CD4.T.cells"
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(13))]<-"Neutrophil.cells"#
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(5))]<-"MoMfDC.cells"
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(2,3,6,10))]<-"Endothelial.cells"
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(17))]<-"Cycling.immune.cells"
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(7))]<-"γδ.T.cells"
ITGB6@meta.data$new.lables[which(ITGB6@meta.data$seurat_clusters %in% c(8))]<-"Tumor.cells"



Idents(ITGB6)<-ITGB6$new.lables
DimPlot(ITGB6,label = T)
table(ITGB6$orig.ident,ITGB6$new.lables)
table(ITGB6$new_group,ITGB6$new.lables)


bulk_dataset<-read.csv('/data/caihua2/ITGB6/上皮分析/HNSC_FPKM_mRNA.csv',row.names = 1)
bulk_survival<-read.csv('/data/caihua2/ITGB6/上皮分析/TCGA-HNSC.survival.csv',row.names = 1)
genelist<-intersect(colnames(bulk_dataset),rownames(bulk_survival))
genelist<-as.data.frame(genelist)
bulk_dataset<-bulk_dataset[,genelist$genelist]
bulk_survival<-bulk_survival[genelist$genelist,]
bulk_dataset[1:3,1:4]#行是基因，列代表细胞（这里不再是单个细胞，而是病人的bulk样本）
head(bulk_survival)
identical(rownames(bulk_survival),colnames(bulk_dataset))
table(bulk_survival)
phenotype <- bulk_survival
colnames(phenotype) <- c("time","status")
bulk_dataset<-as.matrix(bulk_dataset)



infos1 <- Scissor(bulk_dataset, sc_dataset_ITGB6, phenotype,tag =NULL, alpha = 0.05, 
                  family = "cox", #因为表型信息是临床生存信息，故用COX回归模型
                  Save_file = 'Scissor_ITGB6.RData')


#visualize the Scissor selected cells by adding a new annotation in the Seurat object
Scissor_select <- rep(0, ncol(sc_dataset_ITGB6))#创建一个列表，用来表示4102个细胞，
names(Scissor_select) <- colnames(sc_dataset_ITGB6)#给列表中每一个数赋予细胞编号
Scissor_select[infos1$Scissor_pos] <- 1#被选为Scissor+的细胞赋值为1
Scissor_select[infos1$Scissor_neg] <- 2#被选为Scissor-的细胞赋值为2
sc_dataset_ITGB6 <- AddMetaData(sc_dataset_ITGB6, metadata = Scissor_select, col.name = "scissor")#将表示4102个细胞分类的列表添加到sc_dataset_ITGB6这个Seurat对象中
DimPlot(sc_dataset_ITGB6, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))#可视化


infos2 <- Scissor(bulk_dataset, sc_dataset_ITGB6, phenotype, tag = NULL, alpha = 0.134, 
                  family = "cox", Load_file = "Scissor_ITGB6.RData")

Scissor_select <- rep(0, ncol(sc_dataset_ITGB6))
names(Scissor_select) <- colnames(sc_dataset_ITGB6)
Scissor_select[infos2$Scissor_pos] <- 1
Scissor_select[infos2$Scissor_neg] <- 2
sc_dataset_ITGB6 <- AddMetaData(sc_dataset_ITGB6, metadata = Scissor_select, col.name = "scissor")
p<-DimPlot(sc_dataset_ITGB6, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p
ggsave("kang_epi scissor umap.pdf", p, width = 6, height = 4.5,dpi = 600)
DimPlot(sc_dataset_ITGB6,label = T)
table(sc_dataset_ITGB6$seurat_clusters)



