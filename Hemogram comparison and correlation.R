rm(list=ls())

library(vroom)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(Hmisc)

MC_hemo<- read_csv('data/hemogram plot.csv')%>% gather(`NEUT`: `NLR`, key = "Hemo_feature", value = "Value")

MC_hemo


ggboxplot(MC_hemo, x = "Category", y = "Value", add = "jitter",color = "Category", palette = c("#ee0000", "#3b4992"))+
  facet_wrap(~Hemo_feature, nrow=2, scales="free")+
  background_grid()+labs(title="", x="", y="Propotion")+
  stat_compare_means(method = "wilcox.test", paired=TRUE)



# Correlation #----
MC_hemo2<- read_csv('data/correlation plot.csv')%>% gather(`QMG`: `QOL_15`, key = "Clinical_score", value = "Value")

library("ggpubr")
ggscatter(MC_hemo2, x = "NLR", y = "QMG", color = '#3b4992',
          add = "reg.line", add.params = list(color = "#ee0000", fill = "lightgray"),
          conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NLR", ylab = "MMT score")

## facet
ggscatter(MC_hemo2, x = "NLR", y = "Value", 
          color = "#3b4992",
          facet.by = "Clinical_score", scales = "free_x",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#ee0000", fill = "lightgray"))+
  stat_cor(r.accuracy = 0.001, method = "spearman", label.y = 75)

# TEST #--
MC_hemo3<- read_csv('data/correlation plot NEU and LYM.csv')%>% gather(`QMG`: `QOL_15`, key = "Clinical_score", value = "Value")

ggscatter(MC_hemo3, x = "LYM", y = "Value", 
          color = "#3b4992",
          facet.by = "Clinical_score", scales = "free_x",
          add = "reg.line", conf.int = TRUE, 
          add.params = list(color = "#ee0000", fill = "lightgray"))+
  stat_cor(r.accuracy = 0.001, method = "spearman", label.y = 75)


