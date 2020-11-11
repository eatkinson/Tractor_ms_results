setwd('/Users/elizabeth/Dropbox/Projects/MGH-Broad/Tractor/AdmixedAfUKBB/SelectionCriteriaImages/')

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
#install.packages("mixtools")
library(mixtools)
library(sp)
library(plyr)

AfEurPCs <- read.table('pcs.withinonly.UKBB_AfEur_QCed_lipids.txt', header=T)
#4576 people

plot_pca <- function(dataset, first_pc, second_pc) {
  pc_biplot <- ggplot(dataset, aes_string(x=first_pc, y=second_pc, color="darkblue")) +
    geom_point(alpha=0.2) +
    guides(color=guide_legend(override.aes = list(alpha=1))) +
    theme_classic() +
  return(pc_biplot)
}

p <- ggplot(AfEurPCs, aes_string(x='PC1', y='PC2')) +
  geom_point(alpha=0.2) +
  guides(color=guide_legend(override.aes = list(alpha=1))) +
  theme_classic() + ggtitle("Re-generated PCA for admixed Afr-Eur UKBB participants") 

p <- plot_pca(AfEurPCs, 'PC1', 'PC2')  # -> pc1_2
ggsave('AfEur_ukbb_pca_self.png', p, width=6, height=5)


######
#plot again the outcome of PCA with the 1kG reference

setwd('/Users/elizabeth/Dropbox/Projects/MGH-Broad/Tractor/AdmixedAfUKBB/SelectionCriteriaImages')

library(plyr)
library(ggplot2)
read.table('AFR_EUR.1kg_order.indivs.pops', header = T) -> AFR_EUR.1kg_order.indivs.pops
read.table('JoinedIDs.txt', header=T) -> JoinedIDs
join(AFR_EUR.1kg_order.indivs.pops, JoinedIDs, by="ID", type="right") -> joined
joined[,c(1,3)] -> joined1
write.table(joined1, 'joined1.txt', quote=F, row.names=F, col.names=F)

read.table('UKB_admixed.pong_ind2popdata.txt', header = T) -> UKB_admixed.pong_ind2popdata.txt

#plot up the PCA
pcs <- read.table('pcs.UBK_AfrEur_1kg', head=T)

#join with Pop fields
join(pcs, UKB_admixed.pong_ind2popdata.txt, by="IID", type="inner") -> joined

pca_plot1 <- function(first_pc, second_pc, filename) {
  p_pca <- ggplot(filename, aes_string(x=first_pc, y=second_pc, color='Pop', fill='Pop')) +
    geom_point(alpha=0.6) +
    scale_color_discrete(name = "Pop", breaks="Pop") +
    scale_fill_discrete(name = "Pop", breaks="Pop") + guides(color=guide_legend(override.aes = list(alpha=1))) +
    ggtitle("PCA for admixed Afr-Eur UKBB participants with 1000G AFR and EUR") +
    theme_classic()
  return(p_pca)
}

pca_plot1('PC1', 'PC2', joined)

p <- ggplot(joined, aes_string(x='PC1', y='PC2', color='Pop', fill='Pop')) +
    geom_point(alpha=0.6) +  ggtitle("PCA for admixed Afr-Eur UKBB participants with 1000G AFR and EUR") +
    theme_classic()
save_pdf('PCA_UKB_1kg.pdf', p)

