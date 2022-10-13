#!/usr/bin/Rscript
#$ -j y
#$ -cwd

library(qqman)
phenoname = 'TC'

plot_man <- function(phenoname) {
  phenoname.man <- manhattan(na.omit(read.table(paste(phenoname, '.CI.assoc.linear', sep=""), header=T)), chr = "CHR", bp = "BP", p = "P", snp = "SNP",
                         col = c("gray10", "gray60"), chrlabs = NULL,
                         suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
                         highlight = NULL, logp = TRUE, annotatePval = NULL,
                         annotateTop = TRUE, main = paste(phenoname, "GWAS"))
  return(phenoname.man)
}


plot_qq <- function(phenoname) {
  phenoname.qq <- qq((na.omit(read.table(paste(phenoname, '.CI.assoc.linear', sep=""), header=T)))$P, main = paste(phenoname, 'Q-Q plot'))
  return(phenoname.qq)
}

phenoname.man <- plot_man(phenotype)
phenoname.qq <- plot_qq(phenotype)
save_plot(paste(phenoname.man, '.man.pdf', sep=""), phenoname.man, base_height=7, base_width=14)
save_plot(paste(phenoname.qq, '.qq.pdf', sep=""), phenoname.qq, base_height=7, base_width=7)

plots <- plot_grid(phenoname.man, phenoname.qq, ncol=2, labels='AUTO', align='h')
save_plot(paste(phenoname, '.manqq.pdf', sep=""), plots, base_height=7, base_width=20)
