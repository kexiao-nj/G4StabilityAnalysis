library(ggplot2)
library(cowplot)
library(here)
library(ggsci)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(RColorBrewer)
library(ppcor)
setwd(here::here("analysis"))

# read prepared datasets
cols <- c("index", "chrom", "start", "end", "pqs.position", "eG4Sig", "Stability", "phyloP", "ChromState", "ATACSig", "Type", "TFnum")
in_df <- data.frame(read.csv('../prepared_datasets/K562.tsv', sep = "\t"))
final_df <- in_df[,cols]
in_df <- data.frame(read.csv('../prepared_datasets/HepG2.tsv', sep = "\t"))
final_df <- rbind(final_df, in_df[,cols])
in_df <- data.frame(read.csv('../prepared_datasets/293T.tsv', sep = "\t"))
final_df <- rbind(final_df, in_df[,cols])

final_df$Type <- factor(final_df$Type, levels =c("K562", "HepG2", "293T"))
final_df$ChromState <- factor(final_df$ChromState, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))

Types <- c("K562", "HepG2", "293T")
items <- c('Stability', 'eG4Sig', 'ATACSig', 'TFnum', 'phyloP')

# calculate corr
intens_df <- data.frame('x'=NA, 'y'=NA, 'cor'=NA, 'p'=NA, 'star'=NA, 'Type'=NA)

p_to_star <- function(p) {
  cut(p,
      breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
      labels = c("***", "**", "*", ""))
}

for (t in Types) {
  for (i in 1:length(items)) {
    for (j in 1:i) {
      x <- final_df[(final_df$Type == t), items[i]]
      y <- final_df[(final_df$Type == t), items[j]]
      corr <- cor.test(x,y, method='spearman', exact=FALSE)
      cat(t, ": ", items[i], ", ", items[j], " corr is", corr$estimate, " p value is", corr$p.value,"\n")
      intens_df <- intens_df %>% add_row(x=items[i], y=items[j], cor=corr$estimate, p=corr$p.value, star=p_to_star(corr$p.value), Type=t)
    }
  } 
}


intens_df <- intens_df[!is.na(intens_df$x),]
intens_df[intens_df == 'TFnum'] <- '#TF'
intens_df[intens_df$x == intens_df$y, c('cor', 'p', 'star')] <- NA

items_4cor <- c('Stability', 'eG4Sig', 'ATACSig', '#TF', 'phyloP')
intens_df$x <- factor(intens_df$x, levels = items_4cor)
intens_df$y <- factor(intens_df$y, levels = items_4cor)
intens_df$Type <- factor(intens_df$Type, levels = Types)

g <- ggplot(intens_df, aes(x, y, color = cor)) +
  geom_point(aes(size= abs(cor))) +
  geom_text(aes(label=round(cor,2), x=x, y=y), col ="black", size = 6) +
  # geom_text(aes(label=round(cor,2), x=x, y=y), col ="black", size = 6, angle = 30) +
  geom_text(aes(label=as.character(star), x=y, y=x), col ="black", size = 5) +
  geom_abline(slope = -1, intercept = 6, color = 'gray80') +
  scale_y_discrete(limits = rev(items_4cor)) +
  scale_color_gradient(low = "white", high = "dodgerblue", name = "Corr") +
  theme_bw() +
  scale_size_continuous(range = c(1, 25), guide = FALSE) +
  theme(axis.title = element_blank(), axis.text = element_text(size=15),
    axis.text.x = element_text(angle = 30, vjust = 0.5),
    legend.title = element_text(size = 18, margin = margin(b = 10)),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.key.height = unit(4.5, "line")) +
  facet_grid(~Type)
g
# ggsave("sig-state-cor/all.corr.png", g, height = 5, width = 12)
ggsave(here::here("output-fig/fig2-corr.pdf"), g, height = 5.5, width = 16)
