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

#### Bland-Altman plot (difference plot)
## for raw values
# total eDNA
ba1 <- bland.altman.stats(d.comp$norm_total_qmiseq, d.comp$norm_total_qpcr)
# Enguraulis japonicus eDNA
ba2 <- bland.altman.stats(d.comp$eng_qmiseq, d.comp$eng_qpcr)
# Trachurus japonicus eDNA
ba3 <- bland.altman.stats(d.comp$tra_qmiseq, d.comp$tra_qpcr)

## for log-transformed values
log.d.comp <- as.data.frame(apply(d.comp[,1:9], 2, function(x) log10(x+0.5)))
log.d.comp$slope <- d.comp$slope
# linear regressions
summary(lm(log.d.comp$total_qmiseq ~ log.d.comp$total_qpcr))
cor.s1 <- cor(log.d.comp$total_qmiseq, log.d.comp$total_qpcr)
summary(lm(log.d.comp$eng_qmiseq ~ log.d.comp$eng_qpcr))
cor.s2 <- cor(log.d.comp$eng_qmiseq, log.d.comp$eng_qpcr)
summary(lm(log.d.comp$tra_qmiseq ~ log.d.comp$tra_qpcr))
cor.s3 <- cor(log.d.comp$tra_qmiseq, log.d.comp$tra_qpcr)
# total eDNA
lba1 <- bland.altman.stats(log.d.comp$total_qmiseq, d.comp$norm_total_qpcr)
# Enguraulis japonicus eDNA
lba2 <- bland.altman.stats(log.d.comp$eng_qmiseq, log.d.comp$eng_qpcr)
# Trachurus japonicus eDNA
lba3 <- bland.altman.stats(log.d.comp$tra_qmiseq, log.d.comp$tra_qpcr)

## visualize results
# raw values
fs4a <- bland.altman.plot(d.comp$norm_total_qmiseq, d.comp$norm_total_qpcr, graph.sys = "ggplot2")
fs4a <- fs4a + labs(title = "Total fish DNA") + theme_bw()
fs4b <- bland.altman.plot(d.comp$eng_qmiseq, d.comp$eng_qpcr, graph.sys = "ggplot2")
fs4b <- fs4b + labs(title = "Japanese anchovy") + theme_bw()
fs4c <- bland.altman.plot(d.comp$tra_qmiseq, d.comp$tra_qpcr, graph.sys = "ggplot2")
fs4c <- fs4c + labs(title = "Japanese jack mackerel") + theme_bw()

# log-transformed values 
fs5a <- ggplot(log.d.comp, aes(y = total_qmiseq, x = total_qpcr, colour = slope))
fs5a <- fs5a + geom_point() + labs(title = "Total fish DNA")
fs5a <- PlotStyle4(fs5a) + geom_smooth(method = "lm", color = "black", se = F, size = 0.5) 
fs5a <- PlotStyle(fs5a) + theme(legend.position = "none")
fs5a <- fs5a + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.8, parse = T,
                        label = sprintf("R^2 == %0.3f", cor.s1^2))

fs5b <-ggplot(log.d.comp, aes(y = eng_qmiseq, x = eng_qpcr, color = slope))
fs5b <- fs5b + geom_point() + labs(title = "Japanese anchovy")
fs5b <- PlotStyle4(fs5b) + geom_smooth(method = "lm", color = "black", se = F, size = 0.5) 
fs5b <- PlotStyle(fs5b) + theme(legend.position = "none")
fs5b <- fs5b + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.8, parse = T,
                        label = sprintf("R^2 == %0.3f", cor.s2^2))

fs5c <- ggplot(log.d.comp, aes(y = tra_qmiseq, x = tra_qpcr, color = slope))
fs5c <- fs5c + geom_point() + labs(title = "Japanese jack mackerel")
fs5c <- PlotStyle4(fs5c) + geom_smooth(method = "lm", color = "black", se = F, size = 0.5) 
fs5c <- PlotStyle(fs5c) + theme(legend.position = "none")
fs5c <- fs5c + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.8, parse = T,
                        label = sprintf("R^2 == %0.3f", cor.s3^2))

fs5d <- bland.altman.plot(log.d.comp$total_qmiseq, log.d.comp$total_qpcr, graph.sys = "ggplot2")
fs5d <- fs5d + labs(title = "Total fish DNA") + theme_bw()
fs5e <- bland.altman.plot(log.d.comp$eng_qmiseq, log.d.comp$eng_qpcr, graph.sys = "ggplot2")
fs5e <- fs5e + labs(title = "Japanese anchovy") + theme_bw()
fs5f <- bland.altman.plot(log.d.comp$tra_qmiseq, log.d.comp$tra_qpcr, graph.sys = "ggplot2")
fs5f <- fs5f + labs(title = "Japanese jack mackerel") + theme_bw()



# output figure
dev.off()
quartz(width = 8, height = 3) # quartz function is only for Mac
plot_grid(fs4a, fs4b, fs4c, ncol=3, align="hv", labels=c("a","b","c"))
figS2.name <- file.path(qmiseq.out.00, "FigureS4.png")
quartz.save(figS2.name)

dev.off()
quartz(width = 8.5, height = 6) # quartz function is only for Mac
plot_grid(fs5a, fs5b, fs5c,
          fs5d, fs5e, fs5f, ncol=3, align="hv", labels=c("a","b","c","d","e","f"))
figS3.name <- file.path(qmiseq.out.00, "FigureS5.png")
quartz.save(figS3.name)
