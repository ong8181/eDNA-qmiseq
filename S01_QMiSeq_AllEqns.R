####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.S1 Correction equations for all samples
####

# names for the output folder and workspace file
qmiseq.out.00 <- "00_QMiSeq_Figures"
#dir.create(qmiseq.out.00, showWarnings = FALSE)

# load data
d.sample <- read.csv("data/01_sample_data.csv")
d.reads.all <- read.csv("data/02_miseq_reads.csv")

# scientific names and copy numbers of standard DNA
std.name <- c("Saurogobio.immaculatus",
              "Elopichthys.bambusa",
              "Carassioides.acuminatus",
              "Labeo.coubie",
              "Acanthopsoides.gracilentus")
std.copy.1ul <- c(500,250,100,50,25)

# combine standard DNA data
std.reads <- d.reads.all[,std.name]
std.all <- cbind(d.sample, std.reads)
std.sample <- std.all[std.all$treatment == "sample",]

# show the correction equations for all samples
# (copy numbers of standard DNA v.s. reads of standard DNA)
# load library
library(ggplot2)
library(cowplot)
library(reshape2)

# calculate adjusted R2 of the correction equations
d.a1 <- melt(std.sample[,c(2,4:8)], id.vars = "date", variable.name = "STD.sp")
colnames(d.a1)[3] <- "reads"
d.a1$copy <- NaN
d.a1[d.a1[,"STD.sp"] == std.name[1],"copy"] <- std.copy.1ul[1]
d.a1[d.a1[,"STD.sp"] == std.name[2],"copy"] <- std.copy.1ul[2]
d.a1[d.a1[,"STD.sp"] == std.name[3],"copy"] <- std.copy.1ul[3]
d.a1[d.a1[,"STD.sp"] == std.name[4],"copy"] <- std.copy.1ul[4]
d.a1[d.a1[,"STD.sp"] == std.name[5],"copy"] <- std.copy.1ul[5]

# calculate residual structures
resid.raw <- apply((std.sample[,4:8]), 1, function(x) summary(lm(as.numeric(x)~std.copy.1ul+0))$resid)
resid.df <- as.data.frame(resid.raw)
colnames(resid.df) <- as.character(std.sample$date)
resid.df$copy <- std.copy.1ul

d.a2 <- melt(resid.df, id.vars = "copy")
colnames(d.a2) <- c("copy", "date", "residual")

# visualize all regressions
fs1 <- ggplot(d.a1, aes(x = copy, y = reads, group = date))
fs1 <- fs1 + geom_point() + facet_wrap(~date)
fs1 <- fs1 + geom_smooth(method = "lm", se = F, size = 0.8, formula = y ~ x + 0)
fs1 <- fs1 + xlab(expression(paste("The copy number of standard DNA (copies ", mu, l^-1, ")")))
fs1 <- fs1 + ylab("Sequence reads of standard DNA")
fs1 <- fs1 + theme_bw()

# visualize all residuals
fs2 <- ggplot(d.a2, aes(x = copy, y = residual, group = date))
fs2 <- fs2 + geom_point() + facet_wrap(~date) + geom_hline(yintercept = 0, linetype = 2)
fs2 <- fs2 + xlab(expression(paste("The copy number of standard DNA (copies ", mu, l^-1, ")")))
fs2 <- fs2 + ylab("Residuals of the regression")
fs2 <- fs2 + theme_bw()


# output figure
dev.off()
quartz(height=10, width=10)
fs1 # regression line
figS1.name <- file.path(qmiseq.out.00, "FigureS1.png")
quartz.save(figS1.name)

dev.off()
quartz(height=10, width=10)
fs2 # residual structure
figS2.name <- file.path(qmiseq.out.00, "FigureS2.png")
quartz.save(figS2.name)

