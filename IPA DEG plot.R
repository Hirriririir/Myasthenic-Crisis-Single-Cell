library("readxl")
library(tidyverse)

## All DEG (MC vs. After MC) 
go_enrich <- read_excel('./data/MC_DEG/All DEG.xls') %>% head(20)
go_enrich <- go_enrich[order(go_enrich$Logp, decreasing = FALSE), ]
go_enrich$Pathway <- factor(go_enrich$Pathway, levels = go_enrich$Pathway)


#柱形图，纵坐标是 GO Term，横坐标是各 GO Term 的富集得分（Enrichment_score），颜色按 p 值着色
ggplot(go_enrich, aes(Pathway, Ratio)) +
  geom_col(aes(fill = Zscore), width = 0.5) +
  scale_fill_gradient(high = 'red', low = '#5d7ca8') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Gene ratio (IPA)')+
  ## add percentage labels
  geom_text(aes(label = Logp), colour = "white" ,hjust = 1)



## M crisis  #----
M_DEG_Crisis <- read_excel('./data/MC_DEG/M DEG.xls') %>% head(20)

M_DEG_Crisis <- M_DEG_Crisis[order(M_DEG_Crisis$Logp, decreasing = FALSE), ]
M_DEG_Crisis$Pathway <- factor(M_DEG_Crisis$Pathway, levels = M_DEG_Crisis$Pathway)


ggplot(M_DEG_Crisis, aes(Pathway, Ratio)) +
  geom_col(aes(fill = Zscore), width = 0.5) +
  scale_fill_gradient(high = 'red', low = '#5d7ca8') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Gene ratio (IPA)')+
  ## add percentage labels
  geom_text(aes(label = Logp), colour = "white" ,hjust = 1)


## M3 crisis  #----
M3_DEG_Crisis <- read_excel('./data/MC_DEG/M3 DEG.xls') %>% head(20)

M3_DEG_Crisis <- M3_DEG_Crisis[order(M3_DEG_Crisis$Logp, decreasing = FALSE), ]
M3_DEG_Crisis$Pathway <- factor(M3_DEG_Crisis$Pathway, levels = M3_DEG_Crisis$Pathway)

ggplot(M3_DEG_Crisis, aes(Pathway, Ratio)) +
  geom_col(aes(fill = Zscore), width = 0.5) +
  scale_fill_gradient(high = 'red', low = '#5d7ca8') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Gene ratio (IPA)')+
  ## add percentage labels
  geom_text(aes(label = Logp), colour = "white" ,hjust = 1)


## M3 vs. rest M crisis  #----
M3MR_DEG <- read_excel('./data/MC_DEG/M3 Mrest DEG.xls') %>% head(20)

M3MR_DEG <- M3MR_DEG[order(M3MR_DEG$Logp, decreasing = FALSE), ]
M3MR_DEG$Pathway <- factor(M3MR_DEG$Pathway, levels = M3MR_DEG$Pathway)

ggplot(M3MR_DEG, aes(Pathway, Ratio)) +
  geom_col(aes(fill = Zscore), width = 0.5) +
  scale_fill_gradient(high = 'red', low = '#5d7ca8') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Gene ratio (IPA)')+
  ## add percentage labels
  geom_text(aes(label = Logp), colour = "white" ,hjust = 1)



## T crisis #----
T_DEG_Crisis <- read_excel('./data/MC_DEG/T DEG.xls') %>% head(20)

T_DEG_Crisis <- T_DEG_Crisis[order(T_DEG_Crisis$Logp, decreasing = FALSE), ]
T_DEG_Crisis$Pathway <- factor(T_DEG_Crisis$Pathway, levels = T_DEG_Crisis$Pathway)

ggplot(T_DEG_Crisis, aes(Pathway, Ratio)) +
  geom_col(aes(fill = Zscore), width = 0.5) +
  scale_fill_gradient(high = 'red', low = '#5d7ca8') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Gene ratio (IPA)')+
  ## add percentage labels
  geom_text(aes(label = Logp), colour = "white" ,hjust = 1)

## B crisis #----

B_DEG_Crisis <- read_excel('./data/MC_DEG/B DEG.xls') %>% head(20)

B_DEG_Crisis <- B_DEG_Crisis[order(B_DEG_Crisis$Logp, decreasing = FALSE), ]
B_DEG_Crisis$Pathway <- factor(B_DEG_Crisis$Pathway, levels = B_DEG_Crisis$Pathway)

ggplot(B_DEG_Crisis, aes(Pathway, Ratio)) +
  geom_col(aes(fill = Zscore), width = 0.5) +
  scale_fill_gradient(high = 'red', low = '#5d7ca8') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Gene ratio (IPA)')+
  ## add percentage labels
  geom_text(aes(label = Logp), colour = "white" ,hjust = 1)