######################################################

valcano_data <- data.frame(genes=sc_dataset_kang_scissor.markers$gene, 
                           logFC=sc_dataset_kang_scissor.markers$avg_log2FC, 
                           group=rep("NotSignificant", 
                                     nrow(sc_dataset_kang_scissor.markers)),
                           FDR=sc_dataset_kang_scissor.markers$p_val,
                           stringsAsFactors = F,
                           new_group=sc_dataset_kang_scissor.markers$cluster)
valcano_data$logFC<-ifelse(sc_dataset_kang_scissor.markers$cluster=='0',-(valcano_data$logFC),valcano_data$logFC)
valcano_data[which(valcano_data['FDR'] < 4.1242282940505E-87 & 
                     valcano_data['logFC'] > 0.585),"group"] <- "Increased"
valcano_data[which(valcano_data['FDR'] < 4.1242282940505E-87 &
                     valcano_data['logFC'] < -0.585),"group"] <- "Decreased"

cols = c("darkgrey","#00B2FF","orange")
names(cols) = c("NotSignificant","Increased","Decreased")

valcano_data<-read.csv('valcano_data.csv')
library(ggplot2)
vol1 <- ggplot(valcano_data, aes(x = logFC, y = -log10(FDR), color = group))+
  scale_colour_manual(values = cols) +
  ggtitle(label = "Volcano Plot 1", subtitle = "Scissor1 volcano plot") +
  geom_point(size = 2.5, alpha = 1, na.rm = T) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "right") + 
  xlab(expression(log[2]("logFC"))) + 
  ylab(expression(-log[10]("FDR"))) +
  geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 0.2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -0.2, colour="#990000", linetype="dashed")+ 
  scale_y_continuous(trans = "log1p")

vol1

ggsave("ITGB6上皮scenic火山图 .pdf", vol1, width = 6, height = 4.5,dpi = 600)
ggsave("Scissor 1 火山图 .pdf", vol1, width = 6, height = 4.5,dpi = 600)

