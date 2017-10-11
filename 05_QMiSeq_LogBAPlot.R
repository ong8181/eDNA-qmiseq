####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.5 Log-transformed Bland-Altman plot
####

# names for the output folder and workspace file
qmiseq.out.04 <- "04_QMiSeq_CompFigsOut"
qmiseq.out.05 <- "05_QMiSeq_LogBAPlotOut"
dir.create(qmiseq.out.05, showWarnings = FALSE)
ws.out.04 <- file.path(qmiseq.out.04, "04_QMiSeq_CompFigsOut.RData")
ws.out.05 <- file.path(qmiseq.out.05, "05_QMiSeq_LogBAPlotOut.RData")

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
## for log-transformed values
log.d.comp <- as.data.frame(apply(d.comp[,1:9], 2, function(x) log2(x+0.5)))
log.d.comp$slope <- d.comp$slope
log.d.comp <- log.d.comp[log.d.comp$slope >= 10, ] # include samples w/o PCR inhibition
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
# log-transformed values 
f4a <- ggplot(log.d.comp, aes(y = total_qmiseq, x = total_qpcr, colour = slope))
f4a <- f4a + geom_point() + labs(title = "Total fish DNA")
f4a <- PlotStyle4(f4a) + geom_smooth(method = "lm", color = "black", se = F, size = 0.5) 
f4a <- PlotStyle(f4a) + theme(legend.position = "none")
f4a <- f4a + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.8, parse = T,
                        label = sprintf("R^2 == %0.3f", cor.s1^2))

f4b <-ggplot(log.d.comp, aes(y = eng_qmiseq, x = eng_qpcr, color = slope))
f4b <- f4b + geom_point() + labs(title = "Japanese anchovy")
f4b <- PlotStyle4(f4b) + geom_smooth(method = "lm", color = "black", se = F, size = 0.5) 
f4b <- PlotStyle(f4b) + theme(legend.position = "none")
f4b <- f4b + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.8, parse = T,
                        label = sprintf("R^2 == %0.3f", cor.s2^2))

f4c <- ggplot(log.d.comp, aes(y = tra_qmiseq, x = tra_qpcr, color = slope))
f4c <- f4c + geom_point() + labs(title = "Japanese Jack meckerel")
f4c <- PlotStyle4(f4c) + geom_smooth(method = "lm", color = "black", se = F, size = 0.5) 
f4c <- PlotStyle(f4c) + theme(legend.position = "none")
f4c <- f4c + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.8, parse = T,
                        label = sprintf("R^2 == %0.3f", cor.s3^2))

f4d <- bland.altman.plot(log.d.comp$total_qmiseq, log.d.comp$total_qpcr, graph.sys = "ggplot2")
f4d <- f4d + labs(title = "Total fish DNA") + theme_bw()
f4e <- bland.altman.plot(log.d.comp$eng_qmiseq, log.d.comp$eng_qpcr, graph.sys = "ggplot2")
f4e <- f4e + labs(title = "Japanese anchovy") + theme_bw()
f4f <- bland.altman.plot(log.d.comp$tra_qmiseq, log.d.comp$tra_qpcr, graph.sys = "ggplot2")
f4f <- f4f + labs(title = "Japanese Jack meckerel") + theme_bw()



# output figure
dev.off()
quartz(width = 8.8, height = 6) # quartz function is only for Mac
plot_grid(f4a, f4b, f4c,
          f4d, f4e, f4f, ncol=3, align="hv", labels=c("a","b","c","d","e","f"))
fig4.name <- file.path(qmiseq.out.00, "Figure4.png")
quartz.save(fig4.name)

# save workspace
save.image(ws.out.05)
