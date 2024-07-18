rownames(tcga_clin)<-tcga_clin$sample_id
tcga_clin$group<-ifelse(as.numeric(substr(rownames(tcga_clin),14,15)) < 10,'tumor','normal')
table(tcga_clin$group)
tcga_clin_tumor<-subset(tcga_clin,group=='tumor')
gsva_matrix_TCGA_tumor<-gsva_matrix[rownames(tcga_clin_tumor),]
gsva_matrix_TCGA_tumor<-as.data.frame(gsva_matrix_TCGA_tumor)
rownames(gsva_matrix_TCGA_tumor)<-rownames(tcga_clin_tumor)
tcga_clin_tumor$score<-gsva_matrix_TCGA_tumor$gsva_matrix_TCGA_tumor
tcga_clin_tumor$new_group<-ifelse(tcga_clin_tumor$score>median(tcga_clin_tumor$score),'high','low')
table(tcga_clin_tumor$new_group2)



library(coin)
library(survminer)
library(survival)
fit <- survfit(Surv(OS.time,OS) ~ tcga_clin_tumor_early$new_group, data = tcga_clin_tumor_early)
d <- ggsurvplot_list(fit,data = tcga_clin_tumor_early,pval = T,##是否添加P值
                     conf.int = F,### 是否添加置信区间
                     risk.table.col = "strata", 
                     ###linetype = "strata",
                     surv.median.line = "hv", # 是否添加中位生存线
                     risk.table.y.text.col = F,risk.table.y.text = FALSE,
                     ggtheme = theme_bw()+theme(legend.text = element_text(colour = c("red", "blue")))
                     +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                     +theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),axis.title.y.left = element_text(size = 16,face = "bold",vjust = 1),axis.title.x.bottom = element_text(size = 16,face = "bold",vjust = 0))
                     +theme(axis.text.x.bottom = element_text(size = 12,face = "bold",vjust = -0.8,colour = "black"))
                     +theme(axis.text.y.left = element_text(size = 12,face = "bold",vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                     +theme(legend.title = element_text(face = "bold",family = "Times",colour = "black",size = 12))
                     +theme(legend.text = element_text(face = "bold",family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                     palette = c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                     xlab = "Days")##随访的时间时天，就写Days，月份就写Months

d









