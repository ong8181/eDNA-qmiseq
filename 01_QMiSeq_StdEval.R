####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.1 Evaluations of the quality of correction equations
####

# names for the output folder and workspace file
qmiseq.out.00 <- "00_QMiSeq_Figures"
qmiseq.out.01 <- "01_QMiSeq_StdEvalOut"
dir.create(qmiseq.out.00, showWarnings = FALSE)
dir.create(qmiseq.out.01, showWarnings = FALSE)
ws.out.01 <- file.path(qmiseq.out.01, "01_QMiSeq_StdEvalOut.RData")

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

# calculate adjusted R2 of the correction equations
cor.raw <- apply((std.sample[,4:8]), 1, function(x) summary(lm(as.numeric(x)~std.copy.1ul+0))$r.squared)
slope.raw <- apply((std.sample[,4:8]), 1, function(x) summary(lm(as.numeric(x)~std.copy.1ul+0))$coefficients[1])
cor.slope.d <- data.frame(date = std.sample$date,
                          treatment=std.sample$treatment,
                          r2 = cor.raw,
                          slope = slope.raw)


# visualize the results
# load library
library(ggplot2)
library(cowplot)
library(reshape2)

# load helper functions
source("functions/HelperFuncs.R")

# make ggplot objects
# Figure 2a (max, median, and minimum slope regressions)
max.slope <- as.numeric(c(std.sample[which.max(slope.raw), 4:8], 0))
med.slope <- as.numeric(c(std.sample[which.min(abs(slope.raw - median(slope.raw))), 4:8], 0))
min.slope <- as.numeric(c(std.sample[which.min(slope.raw), 4:8], 0))

d.f2a <- melt(data.frame(copy = c(std.copy.1ul, 0),
                    Max.slope = max.slope,
                    Med.slope = med.slope,
                    Min.slope = min.slope),
              id.vars = "copy")

f2a <- ggplot(d.f2a, aes(x = copy, y = value, group = variable, label = variable, colour = variable))
f2a <- f2a + geom_point(size = 2) + xlim(-20,520) + scale_color_manual(name = "Regression slope", values = c("red3", "darkred", "black"))
f2a <- f2a + geom_smooth(method = "lm", size = 0.5, se = F)
f2a <- f2a + xlab(expression(paste("Copy numbers of standard DNA (", {µl}^-1, ")")))
f2a <- f2a + ylab("Sequence reads")
f2a <- f2a + scale_y_continuous(labels = scales::comma)
f2a <- PlotStyle(f2a) + theme(legend.position = c(0.25, 0.75))

# Figure 2b
read.std <- melt(data.frame(std.sample[,c(2,4:8)],
                       slope = cor.slope.d$slope),
                 id.var = c("date", "slope"))
read.std$copy <- NaN
read.std[read.std$variable == std.name[1], "copy"] <- std.copy.1ul[1]
read.std[read.std$variable == std.name[2], "copy"] <- std.copy.1ul[2]
read.std[read.std$variable == std.name[3], "copy"] <- std.copy.1ul[3]
read.std[read.std$variable == std.name[4], "copy"] <- std.copy.1ul[4]
read.std[read.std$variable == std.name[5], "copy"] <- std.copy.1ul[5]
read.std$copy_fac <- factor(read.std$copy, levels=c("25","50","100","250","500"))

f2b <- ggplot(read.std, aes(x = copy_fac, y = value, color = slope))
f2b <- f2b + geom_jitter(shape=16, position = position_jitter(0.2), size = 0.7, alpha=0.5)
f2b <- PlotStyle(f2b) +  scale_colour_gradient2(low = "black", mid = "red3", high = "red3", midpoint = 30)
f2b <- f2b  + xlab(expression(paste("Copy numbers of standard DNA (", {µl}^-1, ")")))
f2b <- f2b + ylab("Sequence reads") + theme(legend.position = c(0.2, 0.65))
f2b <- f2b + scale_y_continuous(labels = scales::comma)

# Figure 2c
d.f2c <- data.frame(r2 = cor.slope.d$r2)
f2c <- ggplot(d.f2c, aes(x=r2))
f2c <- f2c + geom_histogram(stat = "bin", binwidth = 0.04, colour = "black", fill = "gray")
f2c <- PlotStyle(f2c) + xlim(0,1) + geom_hline(yintercept = 0)
f2c <- f2c + xlab(expression(paste("Adjusted ",{R}^2))) + ylab("Count")

# Figure 2d
d.f2d <- data.frame(slope = cor.slope.d$slope)
f2d <- ggplot(d.f2d, aes(x=slope))
f2d <- f2d + geom_histogram(stat = "bin", binwidth = 4, colour = "black", fill = "gray")
f2d <- PlotStyle(f2d) + xlim(-3,60) + geom_hline(yintercept = 0)
f2d <- f2d + xlab("Slope") + ylab("Count")

# output figures
quartz(width = 8.5, height = 6) # quartz function is only for Mac
plot_grid(f2a, f2b, f2c, f2d, ncol=2, align="hv", labels=c("a","b","c","d"))
fig2.name <- file.path(qmiseq.out.00, "Figure2.png")
quartz.save(fig2.name)

# save results
save.image(ws.out.01)
