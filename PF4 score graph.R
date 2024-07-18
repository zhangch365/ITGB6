BiocManager::install("tinyarray")
BiocManager::install("factoextra")
BiocManager::install("factoextra")
library(tinyarray)
library(factoextra)
library(FactoMineR)
library(survival)
setwd('/data/caihua4/TCGA')



df_risk <- GEO_expr[c("LGALS3BP",
                      "ADGRE1",
                      "CX3CR1",
                      "PF4",
                      "C3AR1",
                      "MRC1",
                      "TREM2",
                      "SLAMF9",
                      "CLEC4A"),]
df_risk <- as.data.frame(t(df_risk))
identical(rownames(df_risk),rownames(TCGA_scissor_sur))
dat_risk <- cbind(df_risk, TCGA_scissor_sur)

s=Surv(time, event) ~LGALS3BP+
  ADGRE1+
  CX3CR1+
  PF4+
  C3AR1+
  MRC1+
  TREM2+
  SLAMF9+
  CLEC4A
write.csv(dat_risk,'dat_risk.csv')
dat_risk<-read.csv('dat_risk.csv',row.names = 1)
model <- coxph(s, data = dat_risk )
summary(model,data=dat)


RiskScore<-predict(model,type = "risk")
names(RiskScore) = rownames(dat_risk)


fp <- RiskScore
phe<-dat_risk
fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
fp_dat$riskgroup= ifelse(fp_dat$fp>=median(fp_dat$fp),'high','low')

sur_dat=data.frame(patientid=1:length(fp),time=phe[names(sort(fp)),'time'],event=phe[names(sort(fp )),'event']) 
sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
exp_dat=dat_risk[names(sort(fp)),(ncol(dat_risk)-7):ncol(dat_risk)]
fp_dat$fp <- log(fp_dat$fp)

library(ggplot2)
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
  scale_colour_manual(values = c("red","green"))+
  theme_bw()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
  geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
p1

p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("red","green"))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
p2

fp_dat$riskgroup
dim(df_risk)
an = data.frame(group = fp_dat$riskgroup,
                row.names = rownames(df_risk))
library(pheatmap)
mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)
tmp=t(scale(df_risk))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
p3 <- pheatmap(tmp,col= mycolors,annotation_col = an,show_colnames = F,cluster_cols = F)
p3
dev.new()

library(ggplotify)
plots = list(p1,p2,as.ggplot(as.grob(p3)))
library(gridExtra)
lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) #布局矩阵
grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
dev.new()

############


risk_sv <- as.data.frame(RiskScore)
risk_sv$barcode <- rownames(risk_sv)

identical(rownames(TCGA_scissor_sur),rownames(risk_sv))
library(survminer)
#####
# write.csv(TCGA_surv,'TCGA_surv.csv')
TCGA_surv <- cbind(TCGA_scissor_sur,RiskScore)
TCGA_surv$time <- as.numeric(TCGA_surv$time)
my.surv <-Surv(as.numeric(TCGA_surv$time), TCGA_surv$event)
TCGA_surv$group=ifelse(TCGA_surv$RiskScore > median(TCGA_surv$RiskScore),"high","low")
kmfit1 <- survfit(my.surv~group, data = TCGA_surv)
ggsurvplot(kmfit1, conf.int = F, pval = T, risk.table = T, ncesor.plot = T)
p<-ggsurvplot(kmfit1,palette = c("#FF3030", "#0000FF"),
              linetype ="solid",##linetype = "solid"
              xlab ="Time in Day", #x轴标题
              surv.median.line = "hv", #增加中位生存时间
              #risk.table = TRUE,# 添加风险表
              title = "Overall survival",
              pval = TRUE, # 添加P值
              pval.coord=c(15,0.25),
              legend.title = "",# 设置图例标题，这里设置不显示标题，用空格替代
              legend = c(0.8,0.75), # 指定图例位置
              legend.labs = c("Pf4_high", "Pf4_low") ,# 指定图例分组标签
              ggtheme =theme_survminer())###theme_test()
p



library(tinyarray)

risk_sv$riskgroup= ifelse(risk_sv$RiskScore>=median(risk_sv$RiskScore),'high','low')
group_list = risk_sv$riskgroup #定义组别
group_list = factor(group_list,levels = c("low","high")) # 转换为因子
table(group_list) #查看各自的数量
write.csv(df_risk,'df_risk.csv')
df_risk<-read.csv('df_risk.csv',row.names = 1)
p<-draw_boxplot(t(scale(df_risk)),group_list,color = c('#00B2FF','orange'),width = 0.7)
p
ggsave("TCGA pancancer pf4signature .pdf", p, width = 6, height = 4.5,dpi = 600)

risk_sv$response<-GEO_clin$BOR_binary
write.csv(risk_sv,'risk_sv.csv')

gene_set <-read.csv('mmc1.csv')
bg_genes <- split(as.matrix(gene_set)[,1], gene_set[,2])
gsva_matrix_HNSC<- gsva(as.matrix(GEO_expr), bg_genes,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix_HNSC<-as.data.frame(t(gsva_matrix_HNSC))
risk_sv$score<-gsva_matrix_HNSC$`Pf4+Macrophages`
risk_sv$riskgroup= ifelse(risk_sv$score>=median(risk_sv$score),'high','low')
write.csv(risk_sv,'risk_sv.csv')


TCGA_tRNA_risk <- cbind(dat_risk, risk_sv)
write.csv(TCGA_tRNA_risk,'TCGA_score.csv')






TCGA_surv <- cbind(TCGA_scissor_sur,RiskScore)
TCGA_surv$time <- as.numeric(TCGA_surv$time)
my.surv <-Surv(as.numeric(TCGA_surv$time), TCGA_surv$event)
identical(rownames(TCGA_surv),rownames(gsva_matrix_HNSC))
TCGA_surv$score<-gsva_matrix_HNSC$`Pf4+Macrophages`
TCGA_surv$group=ifelse(TCGA_surv$score > median(TCGA_surv$score),"high","low")
kmfit1 <- survfit(my.surv~group, data = TCGA_surv)
ggsurvplot(kmfit1, conf.int = F, pval = T, risk.table = T, ncesor.plot = T)
p<-ggsurvplot(kmfit1,palette = c("#FF3030", "#0000FF"),
              linetype ="solid",##linetype = "solid"
              xlab ="Time in Day", #x轴标题
              surv.median.line = "hv", #增加中位生存时间
              #risk.table = TRUE,# 添加风险表
              title = "Overall survival",
              pval = TRUE, # 添加P值
              pval.coord=c(15,0.25),
              legend.title = "",# 设置图例标题，这里设置不显示标题，用空格替代
              legend = c(0.8,0.75), # 指定图例位置
              legend.labs = c("Pf4_high", "Pf4_low") ,# 指定图例分组标签
              ggtheme =theme_survminer())###theme_test()
p









