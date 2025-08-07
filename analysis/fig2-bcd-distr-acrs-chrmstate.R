library(ggplot2)
library(cowplot)
library(here)
library(ggsci)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ppcor)
setwd(here::here("analysis"))

# emissions_18_core_K27ac.txt is from the "extended 18-state" model which was trained by ChromHMM, and could be downloaded from the Roadmap epigenomics project.
infile <- "../prepared_datasets/emissions_18_core_K27ac.txt"
emission_df <- data.frame(read.csv(infile, sep = "\t"))

colnames(emission_df)[1] <- str_split_fixed(colnames(emission_df)[1], '[.]', n = 2)[,1]
id <- colnames(emission_df)[1]
vars <- colnames(emission_df)[-1]

new_df <- melt(emission_df, id.vars=c(id), measure.vars=vars, variable.name='Mark', value.name='Ratio', na.rm = T)
new_df$state <- factor(new_df$state, levels = c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk", "7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts", "13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies"))

# preparing heatmap for emission matrix 
g_emission <- ggplot(new_df, aes(Mark, state, fill=Ratio)) +
  geom_tile() +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  ggtitle("Emissions") + 
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(angle=90,vjust = 0.5),
        axis.title.x = element_blank(), 
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(r = 0)) +
  scale_fill_distiller(palette = "Blues", direction = 1)
g_emission

# preparing plot for state names 
name_df <- data.frame(x = rep(1, 18), state = c("Active TSS", "Flanking TSS", "Flanking TSS Upstream", "Flanking TSS Downstream", "Strong transcription", "Weak transcription", "Genic enhancer1", "Genic enhancer2", "Active Enhancer 1", "Active Enhancer 2", "Weak Enhancer", "ZNF genes & repeats", "Heterochromatin", "Bivalent/Poised TSS", "Bivalent Enhancer", "Repressed PolyComb", "Weak Repressed PolyComb", "Quiescent/Low"), color = c("red", "coral1", "coral1", "coral1", "chartreuse4", "darkgreen", "greenyellow", "greenyellow", "orange", "orange", "yellow", "mediumaquamarine", "paleturquoise", "indianred", "darkkhaki", "gray60", "gray80", "white"))
name_df$state <- factor(name_df$state, levels = name_df$state)
g_name <- ggplot(name_df, aes(x, state, fill=state, color=state)) +
  geom_tile() +
  geom_text(aes(label=state), col ="black", size = 3) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text=element_blank(),
        axis.title = element_blank(), 
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = NA),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0) # Left margin
                             ) +
  scale_fill_manual(values=name_df$color) +
  scale_color_manual(values=name_df$color)
g_name


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


# The distribution plot for K562
typelist <- c('Stability', 'eG4Sig', 'ATACSig', 'phyloP', 'TFnum')
celltype <- 'K562'
median_df <- data.frame(Type=typelist)
median_df$medscore <- 0
for (r in typelist){
  median_df[median_df$Type == r, 'medscore'] <- median(final_df[final_df$Type == celltype & final_df$ChromState == "1_TssA", r])
}

g_Stability<-ggplot(final_df[final_df$Type == celltype,], aes(Stability, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.6))) +
  scale_y_discrete(drop = FALSE) +
  xlab('Stability Level (MM %)') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'Stability', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_Stability

g_score<-ggplot(final_df[final_df$Type == celltype,], aes(eG4Sig, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.7))) +
  scale_y_discrete(drop = FALSE) +
  xlab('eG4 Signal Sig. Score') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'eG4Sig', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_score

g_atac<-ggplot(final_df[final_df$Type == celltype,], aes(ATACSig, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  # scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.5))) +
  scale_y_discrete(drop = FALSE) +
  xlab('ATAC-Seq Signal Intensity') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'ATACSig', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_atac

g_phylop<-ggplot(final_df[final_df$Type == celltype,], aes(phyloP, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  # scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(drop = FALSE) +
  xlab('phyloP Score') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'phyloP', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_phylop

g_tfno<-ggplot(final_df[final_df$Type == celltype,], aes(TFnum, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.3))) +
  scale_y_discrete(drop = FALSE) +
  # scale_y_reverse() +
  xlab('# of Colocalized TFs (Norm.)') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'TFnum', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_tfno

g_all <- plot_grid(g_emission, g_name, g_Stability, g_score, g_atac, g_tfno, g_phylop,
                   nrow = 1,
                   # labels = "AUTO",
                   axis = 'tb',
                   rel_widths = c(0.5, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                   align = "h"
)
g_all
ggsave(here::here("output-fig/fig2-b-sigdistr-K562.pdf"), g_all, height = 4.5, width = 16)



# The distribution plot for HepG2
typelist <- c('Stability', 'eG4Sig', 'ATACSig', 'phyloP', 'TFnum')
celltype <- 'HepG2'
median_df <- data.frame(Type=typelist)
median_df$medscore <- 0
for (r in typelist){
  median_df[median_df$Type == r, 'medscore'] <- median(final_df[final_df$Type == celltype & final_df$ChromState == "1_TssA", r])
}

g_Stability<-ggplot(final_df[final_df$Type == celltype,], aes(Stability, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.6))) +
  scale_y_discrete(drop = FALSE) +
  xlab('Stability Level (MM %)') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'Stability', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_Stability

g_score<-ggplot(final_df[final_df$Type == celltype,], aes(eG4Sig, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.7))) +
  scale_y_discrete(drop = FALSE) +
  xlab('eG4 Signal Sig. Score') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'eG4Sig', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_score

g_atac<-ggplot(final_df[final_df$Type == celltype,], aes(ATACSig, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  # scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.5))) +
  scale_y_discrete(drop = FALSE) +
  xlab('ATAC-Seq Signal Intensity') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'ATACSig', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_atac

g_phylop<-ggplot(final_df[final_df$Type == celltype,], aes(phyloP, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  # scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(drop = FALSE) +
  xlab('phyloP Score') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'phyloP', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_phylop

g_tfno<-ggplot(final_df[final_df$Type == celltype,], aes(TFnum, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.3))) +
  scale_y_discrete(drop = FALSE) +
  # scale_y_reverse() +
  xlab('# of Colocalized TFs (Norm.)') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'TFnum', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_tfno

g_all <- plot_grid(g_emission, g_name, g_Stability, g_score, g_atac, g_tfno, g_phylop,
                   nrow = 1,
                   # labels = "AUTO",
                   axis = 'tb',
                   rel_widths = c(0.5, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                   align = "h"
)
g_all
ggsave(here::here("output-fig/fig2-c-sigdistr-HepG2.pdf"), g_all, height = 4.5, width = 16)



# The distribution plot for 293T
typelist <- c('Stability', 'eG4Sig', 'ATACSig', 'phyloP', 'TFnum')
celltype <- '293T'
median_df <- data.frame(Type=typelist)
median_df$medscore <- 0
for (r in typelist){
  median_df[median_df$Type == r, 'medscore'] <- median(final_df[final_df$Type == celltype & final_df$ChromState == "1_TssA", r])
}

g_Stability<-ggplot(final_df[final_df$Type == celltype,], aes(Stability, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.6))) +
  scale_y_discrete(drop = FALSE) +
  xlab('Stability Level (MM %)') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'Stability', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_Stability

g_score<-ggplot(final_df[final_df$Type == celltype,], aes(eG4Sig, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.7))) +
  scale_y_discrete(drop = FALSE) +
  xlab('eG4 Signal Sig. Score') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'eG4Sig', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_score

g_atac<-ggplot(final_df[final_df$Type == celltype,], aes(ATACSig, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  # scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.5))) +
  scale_y_discrete(drop = FALSE) +
  xlab('ATAC-Seq Signal Intensity') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'ATACSig', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_atac

g_phylop<-ggplot(final_df[final_df$Type == celltype,], aes(phyloP, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  # scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(drop = FALSE) +
  xlab('phyloP Score') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'phyloP', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_phylop

g_tfno<-ggplot(final_df[final_df$Type == celltype,], aes(TFnum, ChromState)) +
  geom_jitter(alpha=0.5, size=0.2, color='gray') +
  geom_boxplot(aes(color=ChromState, fill = ChromState), outlier.size = 0.1,
               lwd=0.3, alpha=0.5) +
  scale_x_continuous(expand = expansion(mult = c(0, -0.3))) +
  scale_y_discrete(drop = FALSE) +
  # scale_y_reverse() +
  xlab('# of Colocalized TFs (Norm.)') +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text.x=element_text(vjust = 0.5),
        axis.text.y=element_blank(), 
        legend.position = 'none',
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 0) # Left margin
                             ) +
  geom_vline(aes(xintercept = median_df[median_df$Type == 'TFnum', 'medscore']), color='lightblue4', linetype=2) +
  scale_color_manual(values=name_df$color, drop = FALSE) +
  scale_fill_manual(values=name_df$color, drop = FALSE)
g_tfno

g_all <- plot_grid(g_emission, g_name, g_Stability, g_score, g_atac, g_tfno, g_phylop,
                   nrow = 1,
                   # labels = "AUTO",
                   axis = 'tb',
                   rel_widths = c(0.5, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
                   align = "h"
)
g_all
ggsave(here::here("output-fig/fig2-d-sigdistr-293T.pdf"), g_all, height = 4.5, width = 16)
