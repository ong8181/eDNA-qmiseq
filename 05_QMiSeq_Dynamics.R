####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.5 Illustrations of temporal dynamics of eDNA
####

# names for the output folder and workspace file
qmiseq.out.04 <- "04_QMiSeq_CompFigsOut"
qmiseq.out.05 <- "05_QMiSeq_DynamicsOut"
dir.create(qmiseq.out.05, showWarnings = FALSE)
ws.out.04 <- file.path(qmiseq.out.04, "04_QMiSeq_CompFigsOut.RData")
ws.out.05 <- file.path(qmiseq.out.05, "05_QMiSeq_DynamicsOut.RData")

# load workspace
load(ws.out.04)

# load library
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)

# add date information
d.comp$date <- as.POSIXct(as.character(na.omit(cor.slope.d)$date), format="%Y%m%d")

# simple time series
# total eDNA (qMiSeq and qPCR)
d.comp$norm_total_qmiseq <- as.numeric(scale(d.comp$total_qmiseq))
d.comp$norm_total_qpcr <- as.numeric(scale(d.comp$total_qpcr))
d.comp.ts1 <- melt(d.comp[,c("date", "norm_total_qmiseq", "norm_total_qpcr")],
                   id.var = "date")
f4a <- ggplot(d.comp.ts1, aes(x = date, y = value, group = variable))
f4a <- f4a + geom_line(aes(linetype = variable), size = 0.5)
f4a <- PlotStyle(f4a) + scale_linetype_discrete(name = NULL, labels = c("qMiSeq", "qPCR"))
f4a <- f4a + theme(legend.position = c(0.1, 0.8), legend.text = element_text(size = 12))
f4a <- f4a + ylab(expression(paste("Normalized eDNA (copies ", {µl}^-1, ")")))
f4a <- f4a + xlab("Sampling date") + scale_x_datetime(labels = date_format('%Y-%m'))

# Engraulis japonicus (qMiSeq and qPCR)
d.comp.ts2 <- melt(d.comp[,c("date", "eng_qmiseq", "eng_qpcr")],
                   id.var = "date")
f4b <- ggplot(d.comp.ts2, aes(x = date, y = value, group = variable))
f4b <- f4b + geom_line(aes(linetype = variable), size = 0.5)
f4b <- PlotStyle(f4b) + scale_linetype_discrete(name = NULL, labels = c("qMiSeq", "qPCR"))
f4b <- f4b + theme(legend.position = c(0.1, 0.8), legend.text = element_text(size = 12))
f4b <- f4b + ylab(expression(paste("eDNA (copies ", {µl}^-1, ")")))
f4b <- f4b + xlab("Sampling date") + scale_x_datetime(labels = date_format('%Y-%m'))

# Trachurus japonicus (qMiSeq and qPCR)
d.comp.ts3 <- melt(d.comp[,c("date", "tra_qmiseq", "tra_qpcr")],
                   id.var = "date")
f4c <- ggplot(d.comp.ts3, aes(x = date, y = value, group = variable))
f4c <- f4c + geom_line(aes(linetype = variable), size = 0.5)
f4c <- PlotStyle(f4c) + scale_linetype_discrete(name = NULL, labels = c("qMiSeq", "qPCR"))
f4c <- f4c + theme(legend.position = c(0.1, 0.8), legend.text = element_text(size = 12))
f4c <- f4c + ylab(expression(paste("eDNA (copies ", {µl}^-1, ")")))
f4c <- f4c + xlab("Sampling date") + scale_x_datetime(labels = date_format('%Y-%m'))

# output figure
dev.off()
quartz(width = 6, height = 7) # quartz function is only for Mac
plot_grid(f4a, f4b, f4c, ncol=1, align="hv", labels=c("a","b","c"))
fig4.name <- file.path(qmiseq.out.00, "Figure4.png")
quartz.save(fig4.name)

# save workspace
save.image(ws.out.05)
