####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.4 Figures for the comparisons of qMiSeq and qPCR
####

# names for the output folder and workspace file
qmiseq.out.03 <- "03_QMiSeq_ComparisonOut"
qmiseq.out.04 <- "04_QMiSeq_CompFigsOut"
dir.create(qmiseq.out.04, showWarnings = FALSE)
ws.out.03 <- file.path(qmiseq.out.03, "03_QMiSeq_ComparisonOut.RData")
ws.out.04 <- file.path(qmiseq.out.04, "04_QMiSeq_CompFigsOut.RData")

# load workspace
load(ws.out.03)

# load library
library(ggplot2)
library(cowplot)
library(reshape2)

# comparisons of qMiSeq and qPCR
# total eDNA, all samples
f3a <- ggplot(d.comp, aes(x = total_qmiseq, y = total_qpcr, colour = slope))
f3a <- PlotStyle2(f3a, lm1.1)
f3a <- PlotStyle(f3a) + theme(legend.position = "none")
f3a <- f3a + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                      label = sprintf("R^2 == %0.3f", cor1.1^2))
# total eDNA, exclude ourliers
f3b <- ggplot(d.comp, aes(x = total_qmiseq, y = total_qpcr, colour = slope))
f3b <- PlotStyle2(f3b, lm3.1) + xlim(0, 500) + ylim(0, 3000)
f3b <- PlotStyle(f3b) + theme(legend.position = "none")
f3b <- f3b + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                      label = sprintf("R^2 == %0.3f", cor3.1^2))

# Engraulis japonicus, all samples
f3c <- ggplot(d.comp, aes(x = eng_qmiseq, y = eng_qpcr, colour = slope))
f3c <- PlotStyle2(f3c, lm1.2) + geom_abline(slope = 1, intercept = 0, linetype = 2)
f3c <- PlotStyle(f3c) + theme(legend.position = "none")
f3c <- f3c + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                      label = sprintf("R^2 == %0.3f", cor1.2^2))
# Engraulis japonicus, exclude ourliers
f3d <- ggplot(d.comp, aes(x = eng_qmiseq, y = eng_qpcr, colour = slope))
f3d <- PlotStyle2(f3d, lm3.2) + xlim(0, 150) + ylim(0, 150) + geom_abline(slope = 1, intercept = 0, linetype = 2)
f3d <- PlotStyle(f3d) + theme(legend.position = "none")
f3d <- f3d + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                      label = sprintf("R^2 == %0.3f", cor3.2^2))

# Trachurus japonicus, all samples
f3e <- ggplot(d.comp, aes(x = tra_qmiseq, y = tra_qpcr, colour = slope))
f3e <- PlotStyle2(f3e, lm1.3) + geom_abline(slope = 1, intercept = 0, linetype = 2)
f3e <- PlotStyle(f3e) + theme(legend.position = "none")
f3e <- f3e + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                      label = sprintf("R^2 == %0.3f", cor1.3^2))
# Trachurus japonicus, exclude ourliers
f3f <- ggplot(d.comp, aes(x = tra_qmiseq, y = tra_qpcr, colour = slope))
f3f <- PlotStyle2(f3f, lm3.3) + xlim(0, 20) + ylim(0, 20) + geom_abline(slope = 1, intercept = 0, linetype = 2)
f3f <- PlotStyle(f3f) + theme(legend.position = "none")
f3f <- f3f + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, parse = T,
                      label = sprintf("R^2 == %0.3f", cor3.3^2))

# output figure
dev.off()
quartz(width = 7.5, height = 9.5) # quartz function is only for Mac
plot_grid(f3a, f3b, f3c, f3d, f3e, f3f, ncol=2, align="hv", labels=c("a","b","c","d","e","f"))
fig3.name <- file.path(qmiseq.out.00, "Figure3.png")
quartz.save(fig3.name)

# save workspace
save.image(ws.out.04)
