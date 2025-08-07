library(ggplot2)
library(cowplot)
library(ggsci)
library(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ExperimentHub)
library(clusterProfiler)
setwd(here::here("analysis"))


# read prepared datasets
cols <- c("index", "chrom", "start", "end", "pqs.position", "eG4Sig", "Stability", "Type", "annotation", "geneId", "transcriptId", "distanceToTSS")
in_df <- data.frame(read.csv('../prepared_datasets/K562.tsv', sep = "\t"))
in_df <- in_df[,cols]
in_df$mmpctlvl <- 'Low'
in_df[in_df$Stability > median(in_df$Stability), 'mmpctlvl'] <- 'High'
in_df$mmpctlvl <- factor(in_df$mmpctlvl, levels = c('Low', 'High'))
final_df <- in_df[,c(cols, 'mmpctlvl')]

in_df <- data.frame(read.csv('../prepared_datasets/HepG2.tsv', sep = "\t"))
in_df$mmpctlvl <- 'Low'
in_df[in_df$Stability > median(in_df$Stability), 'mmpctlvl'] <- 'High'
in_df$mmpctlvl <- factor(in_df$mmpctlvl, levels = c('Low', 'High'))
final_df <- rbind(final_df, in_df[,c(cols, 'mmpctlvl')])

in_df <- data.frame(read.csv('../prepared_datasets/293T.tsv', sep = "\t"))
in_df$mmpctlvl <- 'Low'
in_df[in_df$Stability > median(in_df$Stability), 'mmpctlvl'] <- 'High'
in_df$mmpctlvl <- factor(in_df$mmpctlvl, levels = c('Low', 'High'))
final_df <- rbind(final_df, in_df[,c(cols, 'mmpctlvl')])

final_df$Type <- factor(final_df$Type, levels =c("K562", "HepG2", "293T"))



high_vec_K562_vinc <- final_df[(final_df$distanceToTSS > -200) & (final_df$distanceToTSS < 200) & final_df$Type=='K562' & final_df$mmpctlvl=='High', 'geneId']
low_vec_K562_vinc <- final_df[(final_df$distanceToTSS > -200) & (final_df$distanceToTSS < 200) & final_df$Type=='K562' & final_df$mmpctlvl=='Low', 'geneId']

high_vec_HepG2_vinc <- final_df[(final_df$distanceToTSS > -200) & (final_df$distanceToTSS < 200) & final_df$Type=='HepG2' & final_df$mmpctlvl=='High', 'geneId']
low_vec_HepG2_vinc <- final_df[(final_df$distanceToTSS > -200) & (final_df$distanceToTSS < 200) & final_df$Type=='HepG2' & final_df$mmpctlvl=='Low', 'geneId']

high_vec_293T_vinc <- final_df[(final_df$distanceToTSS > -200) & (final_df$distanceToTSS < 200) & final_df$Type=='293T' & final_df$mmpctlvl=='High', 'geneId']
low_vec_293T_vinc <- final_df[(final_df$distanceToTSS > -200) & (final_df$distanceToTSS < 200) & final_df$Type=='293T' & final_df$mmpctlvl=='Low', 'geneId']


geneClsts<-list(low_vec_K562_vinc, low_vec_HepG2_vinc, low_vec_293T_vinc, 
                high_vec_K562_vinc, high_vec_HepG2_vinc, high_vec_293T_vinc)
names(geneClsts)<-c('K562 Low', 'HepG2 Low', '293T Low', 
                    'K562 High', 'HepG2 High', '293T High')

cmpClstGO <- compareCluster(geneClusters = geneClsts, fun = enrichGO,
                            ont = "BP",OrgDb='org.Hs.eg.db')
p1 <- dotplot(cmpClstGO, showCategory = 10, label_format=80, font.size=14) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 15, vjust = 1),
        axis.text.y = element_text(size = 10)) +
  scale_color_material("red", reverse = TRUE)

p1
save_plot(here::here("output-fig/fig5-tss-vincity-GOenrich-4stability.pdf"), p1, base_height=12, base_asp=0.8)


