####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.5 Illustrations of temporal dynamics of eDNA
####

# names for the output folder and workspace file
qmiseq.out.05 <- "05_QMiSeq_LogBAPlotOut"
qmiseq.out.06 <- "06_QMiSeq_DynamicsOut"
dir.create(qmiseq.out.06, showWarnings = FALSE)
ws.out.05 <- file.path(qmiseq.out.05, "05_QMiSeq_LogBAPlotOut.RData")
ws.out.06 <- file.path(qmiseq.out.06, "06_QMiSeq_DynamicsOut.RData")

# load workspace
load(ws.out.05)

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
f5a <- ggplot(d.comp.ts1, aes(x = date, y = value, group = variable))
f5a <- f5a + geom_line(aes(linetype = variable), size = 0.5)
f5a <- PlotStyle(f5a) + scale_linetype_discrete(name = NULL, labels = c("qMiSeq", "qPCR"))
f5a <- f5a + theme(legend.position = c(0.1, 0.7), legend.text = element_text(size = 12))
f5a <- f5a + ylab(expression(paste("Normalized eDNA (copies ", mu, l^-1, ")")))
f5a <- f5a + xlab("Sampling date") + scale_x_datetime(labels = date_format('%Y-%m'))


# Engraulis japonicus (qMiSeq and qPCR)
d.comp.ts2 <- melt(d.comp[,c("date", "eng_qmiseq", "eng_qpcr")],
                   id.var = "date")
f5b <- ggplot(d.comp.ts2, aes(x = date, y = value, group = variable))
f5b <- f5b + geom_line(aes(linetype = variable), size = 0.5)
f5b <- PlotStyle(f5b) + scale_linetype_discrete(name = NULL, labels = c("qMiSeq", "qPCR"))
f5b <- f5b + theme(legend.position = c(0.1, 0.7), legend.text = element_text(size = 12))
f5b <- f5b + ylab(expression(paste("eDNA (copies ", mu, l^-1, ")")))
f5b <- f5b + xlab("Sampling date") + scale_x_datetime(labels = date_format('%Y-%m'))

# Trachurus japonicus (qMiSeq and qPCR)
d.comp.ts3 <- melt(d.comp[,c("date", "tra_qmiseq", "tra_qpcr")],
                   id.var = "date")
f5c <- ggplot(d.comp.ts3, aes(x = date, y = value, group = variable))
f5c <- f5c + geom_line(aes(linetype = variable), size = 0.5)
f5c <- PlotStyle(f5c) + scale_linetype_discrete(name = NULL, labels = c("qMiSeq", "qPCR"))
f5c <- f5c + theme(legend.position = c(0.1, 0.7), legend.text = element_text(size = 12))
f5c <- f5c + ylab(expression(paste("eDNA (copies ", mu, l^-1, ")")))
f5c <- f5c + xlab("Sampling date") + scale_x_datetime(labels = date_format('%Y-%m'))

# output figure
dev.off()
quartz(width = 6, height = 7) # quartz function is only for Mac
plot_grid(f5a, f5b, f5c, ncol=1, align="hv", labels=c("a","b","c"))
fig5.name <- file.path(qmiseq.out.00, "Figure5.png")
quartz.save(fig5.name)

# save workspace
save.image(ws.out.06)
