####
DefaultAssay(Epithelial) <- "RNA"
emt_features <- c("SERPINE1","TGFBI","MMP10","LAMC2","P4HA2","PDPN","ITGA5","LAMA3","CDH13","TNC","MMP2","EMP3","INHBA","LAMB3","VIM","SEMA3C","PRKCDBP","ANXA5","DHRS7","ITGB1","ACTN1","CXCR7","ITGB6","IGFBP7","THBS1","PTHLH","TNFRSF6B","PDLIM7","CAV1","DKK3","COL17A1","LTBP1",
                  "COL5A2","COL1A1","FHL2","TIMP3","PLAU","LGALS1","PSMD2","CD63","HERPUD1","TPM1","SLC39A14","C1S","MMP1","EXT2","COL4A2","PRSS23","SLC7A8","SLC31A2","ARPC1B","APP","MFAP2","MPZL1","DFNA5","MT2A","MAGED2","ITGA6","FSTL1","TNFRSF12A","IL32","COPB2","PTK7",
                  "OCIAD2","TAX1BP3","SEC13","SERPINH1","TPM4","MYH9","ANXA8L1","PLOD2","GALNT2","LEPREL1","MAGED1","SLC38A5","FSTL3","CD99","F3","PSAP","NMRK1","FKBP9","DSG2","ECM1","HTRA1","SERINC1","CALU","TPST1","PLOD3","IGFBP3","FRMD6","CXCL14","SERPINE2",
                  "RABAC1","TMED9","NAGK","BMP1","ESYT1","STON2","TAGLN","GJA1")

library(Hmisc)
library(ggpubr)
m_emt_features <- capitalize(tolower(emt_features))
m_emt_features
m_emt_features <- list(m_emt_features)
write.csv(m_emt_features,'m_emt_features.csv')


Epithelial <- AddModuleScore(object = Epithelial, features = m_emt_features, name = 'm_emt_features')
colnames(Epithelial@meta.data)

col.num <- length(levels(Epithelial@active.ident))
table(Epithelial@meta.data$m_emt_features1)

my_comparisons <- list( c("0", "3"), c("2", "3"), c("1", "3"), c("4", "3"), c("5", "3"))
table(Epithelial@meta.data$seurat_clusters)

p13 <- ggviolin(Epithelial@meta.data, "seurat_clusters", "m_emt_features1", fill = "seurat_clusters",
                size = 0.176)+             #add = "boxplot", add.params = list(fill = "white",size=0.176,color="black",shape="dose"),
  geom_boxplot(width=0.2,position = position_identity(),fill="white",outlier.size = 0.001,size=0.176)+                    #,notch=T, outlier.colour ="black"
  stat_compare_means(comparisons = my_comparisons,tip.length=0.01,bracket.size=0.176,show.legend=F)+         #用于在ggplot图形中自动添加P值和显著性水平  #step.increase=0.05跟label.y=c(0.5,0.55,0.6,0.65,0.7)差不多   #,label.y=c(0.5,0.55,0.6,0.65,0.7)
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position ="none")+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_line(linetype=1,color="black",size=0.176),axis.line.y=element_line(linetype=1,color="black",size=0.176),axis.line.x=element_line(linetype=1,color="black",size=0.176))+
  #scale_y_continuous(limits = c(-0.15,0.5))+                 #控制 y轴长度
  #scale_fill_manual(values =c('#666699', '#bd88f5',"#CC6600",'#fa2c7b','#AB3282','#D6E7A3','#6666FF','#ffa235','#91d024','#2aaaef', '#FFFF00')) +
  #scale_fill_brewer(palette="Set1")+   ###改颜色的地方，具体查看这个网址https://www.jianshu.com/p/4af3b009df20###
  theme(
    #panel.border = element_rect(color = 'black', fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background   = element_blank(),
    axis.title.x       = element_blank(),
    axis.title.y       = element_blank(),
    axis.ticks         = element_blank(),
    axis.text.x= element_text(hjust=1),
    axis.text.y= element_text(hjust=1))
p13
ggsave("./Epithelial cells/epi_emtscore_seurat_clusters.png",plot = p13,width = 8, height = 6, dpi = 300)
ggsave("./Epithelial cells/epi_emtscore_seurat_clusters.pdf",plot = p13,width = 8, height = 6, dpi = 300)



DefaultAssay(Epithelial) <- "RNA"
cell_cycle_features <- c("CENPA","HIST1H2BC","CCNB2","HMGB2","LOCKD","HMGN2","H2AFZ","HSP90B1","NUCKS1","HMGB3","KNSTRN","HMGB1","FDPS","STMN1","UBE2S","PTMS","CDC20","H1FX","H2AFV","IDI1","CENPE","HNRNPA1",
                  "CDCA3","CDCA8","DYNLL1","EMP1","LSM5","EPCAM","HIST1H1C","TPX2","ANP32E","2700094K13RIK","PTMA","LSM4","CKS2","HNRNPA2B1","HNRNPA3","HMGN5","EIF5A","NDE1","DBI","PDIA6","TMSB10",
                  "HMMR","CALM3","GCAT","CDKN3","YWHAE","DDX39","HMGN1","INSIG1","GM10282","RANGAP1","CBX3","RAD21","KIF23","MARCKS","LBR","SFPQ","HDGF","BIRC5","GT(ROSA)26SOR","ARL2BP","RACGAP1","SNRPF",
                  "CALM1","TMA7","CANX","NPM3","GCSH","NUDCD2","PSMC1","CTNNB1","SAE1","KRT8","DAZAP1","HMGN3","ACAT2","LMNB1","PPIA","RAN","CEP89","PSRC1","HP1BP3","YBX1","IGFBPL1","HMGCS1","RBM3",
                  "AKR1B3","CENPW","PFDN4","CCDC34","GNB1","CBX1","RANGRF","HNRNPDL","SET","PPP1R14B","CENPF","ODC1",
                  "UBE2C","TOP2A","KPNA2","NUSAP1","CCNB1","CENPF","PRC1","CKS2","ARL6IP1","H2AFX","CDC20","HMGB2","CDK1","BIRC5","TUBB4B","HIST1H2AP","CDCA8","TUBA1B","CENPA","LOCKD","TUBA1C","PBK","SPC25",
                  "HIST1H2AE","TPX2","CCNA2","CDCA3","PLK1","MKI67","HN1","UBE2S","HMMR","CENPE","KIF22","KIF23","TUBB5","SMC4","AURKA","CCNB2","H2AFZ","RACGAP1","SPC24","BUB3","CKS1B","TACC3","FAM64A",
                  "TMPO","MIS18BP1","STMN1","KIF11","KNSTRN","KIF20B","CKAP2","CKAP2L","CASC5","CDCA2","HJURP","HIST1H1B","SHCBP1","SMC2","H2AFV","AURKB","KIF20A","H1FX","KIF2C","SGOL1","ECT2","HIST1H1E",
                  "HMGN2","SGOL2A","INCENP","NUCKS1","RANGAP1","CDKN2C","2810417H13RIK","DEPDC1A","ASPM","MXD3","RAN","G2E3","SAPCD2","NDE1","LMNB1","REEP4","CCDC34","FZR1","KIFC1","CDKN3","CKAP5",
                  "BUB1B","ARHGAP11A","TRIM59","DLGAP5","HMGB1","ANP32E","PSRC1","CALM2","NDC80","NUF2","CIT",
                  "2810417H13RIK","RRM2","HIST1H2AP","LIG1","TYMS","HELLS","HIST1H2AE","SLBP","PCNA","HIST1H1B","CCNE2","TUBA1B","MCM3","GMNN","RANBP1","TOP2A","MCM6","MCM5","DUT","MCM7",
                  "DEK","UNG","TUBB5","RRM1","TK1","NASP","SMC2","DTL","TIPIN","MCM4","DCTPP1","MCM2","PBK","ATAD2","PRIM1","RPA2","ESCO2","FEN1","HIST1H1E","HIST1H2AN","DNAJC9","DNMT1","CLSPN",
                  "LSM2","ZFP367","RFC5","RFC2","STMN1","H2AFZ","CDK1","SPC24","CHAF1B","CDCA7","DTYMK","ALYREF","SIVA1","ANP32B","CBX5","UHRF1","RAD51","CENPH","ARL6IP6","SYCE2","EZH2","GINS2",
                  "HAT1","CENPK","SRSF7","HMGB2","NAP1L1","FAM111A","USP1","CHAF1A","NRM","ORC6","CDT1","CDK2","HN1L","SSRP1","RFC3","SMC6","TMPO","RFC4","PTGES3","PSIP1","POLA1","NCL","APITD1",
                  "CKS1B","SHMT1","PAICS","SUPT16","RAN","SMC4","RPA3","HJURP","NSG1","HSPD1","NOP56","CENPM"
)
library(Hmisc)
library(ggpubr)
m_cell_cycle_features <- capitalize(tolower(cell_cycle_features))
m_cell_cycle_features
m_cell_cycle_features <- list(m_cell_cycle_features)

Epithelial <- AddModuleScore(object = Epithelial, features = m_cell_cycle_features, name = 'm_cell_cycle_features')
colnames(Epithelial@meta.data)

col.num <- length(levels(Epithelial@active.ident))
table(Epithelial@meta.data$m_cell_cycle_features1)

my_comparisons <- list( c("0", "1"), c("2", "1"), c("1", "3"), c("4", "1"), c("5", "1"))
table(Epithelial@meta.data$seurat_clusters)

p13 <- ggviolin(Epithelial@meta.data, "seurat_clusters", "m_cell_cycle_features1", fill = "seurat_clusters",
                size = 0.176)+             #add = "boxplot", add.params = list(fill = "white",size=0.176,color="black",shape="dose"),
  geom_boxplot(width=0.15,position = position_identity(),fill="white",outlier.size = 0.001,size=0.176)+                    #,notch=T, outlier.colour ="black"
  stat_compare_means(comparisons = my_comparisons,tip.length=0.01,bracket.size=0.176,show.legend=F)+         #用于在ggplot图形中自动添加P值和显著性水平  #step.increase=0.05跟label.y=c(0.5,0.55,0.6,0.65,0.7)差不多   #,label.y=c(0.5,0.55,0.6,0.65,0.7)
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position ="none")+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_line(linetype=1,color="black",size=0.176),axis.line.y=element_line(linetype=1,color="black",size=0.176),axis.line.x=element_line(linetype=1,color="black",size=0.176))+
  #scale_y_continuous(limits = c(-0.15,0.5))+                 #控制 y轴长度
  #scale_fill_manual(values =c('#666699', '#bd88f5',"#CC6600",'#fa2c7b','#AB3282','#D6E7A3','#6666FF','#ffa235','#91d024','#2aaaef', '#FFFF00')) +
  #scale_fill_brewer(palette="Set1")+   ###改颜色的地方，具体查看这个网址https://www.jianshu.com/p/4af3b009df20###
  theme(
    #panel.border = element_rect(color = 'black', fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background   = element_blank(),
    axis.title.x       = element_blank(),
    axis.title.y       = element_blank(),
    axis.ticks         = element_blank(),
    axis.text.x= element_text(hjust=1),
    axis.text.y= element_text(hjust=1))
p13
ggsave("./Epithelial cells/epi_cell_cyclescore_seurat_clusters.png",plot = p13,width = 8, height = 6, dpi = 300)
ggsave("./Epithelial cells/epi_cell_cyclescore_seurat_clusters.pdf",plot = p13,width = 8, height = 6, dpi = 300)

