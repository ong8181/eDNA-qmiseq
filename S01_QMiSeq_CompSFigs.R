####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.S1 Figures for the comparisons of qMiSeq and qPCR
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
source("functions/S1_HelperFuncs.R")

# comparisons of qMiSeq and qPCR
# total eDNA, all samples
fS1a <- ggplot(d.comp, aes(x = total_reads, y = total_qpcr, colour = slope))
fS1a <- PlotStyle3(fS1a, lm2.1)
fS1a <- PlotStyle(fS1a) + theme(legend.position = "none")
fS1a <- fS1a + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor2.1^2))
# total eDNA, exclude ourliers
fS1b <- ggplot(d.comp, aes(x = total_reads, y = total_qpcr, colour = slope))
fS1b <- PlotStyle3(fS1b, lm4.1) + ylim(0, 3000)
fS1b <- PlotStyle(fS1b) + theme(legend.position = "none")
fS1b <- fS1b + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor4.1^2))

# Engraulis japonicus, all samples
fS1c <- ggplot(d.comp, aes(x = eng_reads, y = eng_qpcr, colour = slope))
fS1c <- PlotStyle3(fS1c, lm2.2)
fS1c <- PlotStyle(fS1c) + theme(legend.position = "none")
fS1c <- fS1c + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor2.2^2))
# Engraulis japonicus, exclude ourliers
fS1d <- ggplot(d.comp, aes(x = eng_reads, y = eng_qpcr, colour = slope))
fS1d <- PlotStyle3(fS1d, lm4.2) + xlim(0, 3000) + ylim(0, 150)
fS1d <- PlotStyle(fS1d) + theme(legend.position = "none")
fS1d <- fS1d + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor4.2^2))

# Trachurus japonicus, all samples
fS1e <- ggplot(d.comp, aes(x = tra_reads, y = tra_qpcr, colour = slope))
fS1e <- PlotStyle3(fS1e, lm2.3)
fS1e <- PlotStyle(fS1e) + theme(legend.position = "none")
fS1e <- fS1e + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor2.3^2))
# Trachurus japonicus, exclude ourliers
fS1f <- ggplot(d.comp, aes(x = tra_reads, y = tra_qpcr, colour = slope))
fS1f <- PlotStyle3(fS1f, lm4.3) + xlim(0, 350) + ylim(0, 20)
fS1f <- PlotStyle(fS1f) + theme(legend.position = "none")
fS1f <- fS1f + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                        label = sprintf("R^2 == %0.3f", cor4.3^2))

# output figure
dev.off()
quartz(width = 7.5, height = 9.5) # quartz function is only for Mac
plot_grid(fS1a, fS1b, fS1c, fS1d, fS1e, fS1f, ncol=2, align="hv", labels=c("a","b","c","d","e","f"))
figS1.name <- file.path(qmiseq.out.00, "FigureS1.png")
quartz.save(figS1.name)
