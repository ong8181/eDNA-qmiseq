####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.S3 Bland-Altman plot
####

# names for the output folder and workspace file
qmiseq.out.04 <- "04_QMiSeq_CompFigsOut"
ws.out.04 <- file.path(qmiseq.out.04, "04_QMiSeq_CompFigsOut.RData")

# load workspace
load(ws.out.04)
source("functions/HelperFuncs.R")

# load library
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(BlandAltmanLeh)

# add date information
d.comp$date <- as.POSIXct(as.character(na.omit(cor.slope.d)$date), format="%Y%m%d")
# simple time series
# total eDNA (qMiSeq and qPCR)
d.comp$norm_total_qmiseq <- as.numeric(scale(d.comp$total_qmiseq))
d.comp$norm_total_qpcr <- as.numeric(scale(d.comp$total_qpcr))
d.comp <- d.comp[d.comp$slope >= 10, ] # include samples w/o PCR inhibitio

#### Bland-Altman plot (difference plot)
## for raw values
# total eDNA
ba1 <- bland.altman.stats(d.comp$norm_total_qmiseq, d.comp$norm_total_qpcr)
# Enguraulis japonicus eDNA
ba2 <- bland.altman.stats(d.comp$eng_qmiseq, d.comp$eng_qpcr)
# Trachurus japonicus eDNA
ba3 <- bland.altman.stats(d.comp$tra_qmiseq, d.comp$tra_qpcr)

## visualize results
# raw values
fs4a <- bland.altman.plot(d.comp$norm_total_qmiseq, d.comp$norm_total_qpcr, graph.sys = "ggplot2")
fs4a <- fs4a + labs(title = "Total fish DNA") + theme_bw()
fs4b <- bland.altman.plot(d.comp$eng_qmiseq, d.comp$eng_qpcr, graph.sys = "ggplot2")
fs4b <- fs4b + labs(title = "Japanese anchovy") + theme_bw()
fs4c <- bland.altman.plot(d.comp$tra_qmiseq, d.comp$tra_qpcr, graph.sys = "ggplot2")
fs4c <- fs4c + labs(title = "Japanese Jack meckerel") + theme_bw()

# output figure
dev.off()
quartz(width = 8, height = 3) # quartz function is only for Mac
plot_grid(fs4a, fs4b, fs4c, ncol=3, align="hv", labels=c("a","b","c"))
figS4.name <- file.path(qmiseq.out.00, "FigureS4.png")
quartz.save(figS4.name)
