####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.S2 Figures for the comparisons of qMiSeq and qPCR
####

# names for the output folder and workspace file
qmiseq.out.04 <- "04_QMiSeq_CompFigsOut"
ws.out.04 <- file.path(qmiseq.out.04, "04_QMiSeq_CompFigsOut.RData")

# load workspace
load(ws.out.04)

# load library
library(ggplot2)
library(cowplot)
library(reshape2)
source("functions/HelperFuncs.R")

# comparisons of qMiSeq and qPCR
# total eDNA, all samples
fs3a <- ggplot(d.comp, aes(x = total_reads, y = total_qpcr, colour = slope))
fs3a <- PlotStyle3(fs3a, lm2.1)
fs3a <- PlotStyle(fs3a) + theme(legend.position = "none",
                                axis.text = element_text(size=9)) + ggtitle("Total fish eDNA")
fs3a <- fs3a + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor2.1^2))
# total eDNA, exclude ourliers
fs3b <- ggplot(d.comp, aes(x = total_reads, y = total_qpcr, colour = slope))
fs3b <- PlotStyle3(fs3b, lm4.1) + ylim(0, 3000)
fs3b <- PlotStyle(fs3b) + theme(legend.position = "none",
                                axis.text = element_text(size=9)) + ggtitle("Total fish eDNA")
fs3b <- fs3b + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor4.1^2))

# Engraulis japonicus, all samples
fs3c <- ggplot(d.comp, aes(x = eng_reads, y = eng_qpcr, colour = slope))
fs3c <- PlotStyle3(fs3c, lm2.2)
fs3c <- PlotStyle(fs3c) + theme(legend.position = "none",
                                axis.text = element_text(size=9)) + ggtitle("Japanese anchovy")
fs3c <- fs3c + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                       label = sprintf("R^2 == %0.3f", cor2.2^2))
# Engraulis japonicus, exclude ourliers
fs3d <- ggplot(d.comp, aes(x = eng_reads, y = eng_qpcr, colour = slope))
fs3d <- PlotStyle3(fs3d, lm4.2) + xlim(0, 3000) + ylim(0, 150)
fs3d <- PlotStyle(fs3d) + theme(legend.position = "none",
                                axis.text = element_text(size=9)) + ggtitle("Japanese anchovy")
fs3d <- fs3d + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor4.2^2))

# Trachurus japonicus, all samples
fs3e <- ggplot(d.comp, aes(x = tra_reads, y = tra_qpcr, colour = slope))
fs3e <- PlotStyle3(fs3e, lm2.3)
fs3e <- PlotStyle(fs3e) + theme(legend.position = "none",
                                axis.text = element_text(size=9)) + ggtitle("Japanese jack mackerel")
fs3e <- fs3e + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor2.3^2))
# Trachurus japonicus, exclude ourliers
fs3f <- ggplot(d.comp, aes(x = tra_reads, y = tra_qpcr, colour = slope))
fs3f <- PlotStyle3(fs3f, lm4.3) + xlim(0, 350) + ylim(0, 20)
fs3f <- PlotStyle(fs3f) + theme(legend.position = "none",
                                axis.text = element_text(size=9)) + ggtitle("Japanese jack mackerel")
fs3f <- fs3f + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor4.3^2))

# output figure
dev.off()
quartz(width = 6.5, height = 8.5) # quartz function is only for Mac
plot_grid(fs3a, fs3b, fs3c, fs3d, fs3e, fs3f, ncol=2, align="hv", labels=c("a","b","c","d","e","f"))
figS3.name <- file.path(qmiseq.out.00, "FigureS3.png")
quartz.save(figS3.name)
