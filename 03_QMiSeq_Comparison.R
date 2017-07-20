####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.3 Comparisons of qMiSeq and qPCR
####

# names for the output folder and workspace file
qmiseq.out.02 <- "02_QMiSeq_ConversionOut"
qmiseq.out.03 <- "03_QMiSeq_ComparisonOut"
dir.create(qmiseq.out.03, showWarnings = FALSE)
ws.out.02 <- file.path(qmiseq.out.02, "02_QMiSeq_ConversionOut.RData")
ws.out.03 <- file.path(qmiseq.out.03, "03_QMiSeq_ComparisonOut.RData")

# load workspace
load(ws.out.02)

# load data
d.calc.read <- read.csv(converted.name2)
d.qpcr <- read.csv("data/03_qPCR_copy.csv")
raw.reads.sample <- d.reads.all[d.sample$treatment == "sample",]

# preparations for comparisons
d.comp <- data.frame(total_qmiseq = rowSums(d.calc.read),
                     total_reads = rowSums(raw.reads.sample),
                     total_qpcr = d.qpcr$MiFish.SYBR,
                     eng_qmiseq = d.calc.read[,"Engraulis.japonicus"],
                     eng_reads = raw.reads.sample[,"Engraulis.japonicus"],
                     eng_qpcr = d.qpcr[,"Engraulis.japonicus"],
                     tra_qmiseq = d.calc.read[,"Trachurus.japonicus"],
                     tra_reads = raw.reads.sample[,"Trachurus.japonicus"],
                     tra_qpcr = d.qpcr[,"Trachurus.japonicus"],
                     slope = cor.slope.d$slope,
                     r2 = cor.slope.d$r2)
d.comp <- na.omit(d.comp) # exclude one NA sample

# simple linear regressions
# qMiSeq v.s. qPCR, all values
cor1.1 <- cor(d.comp$total_qmiseq, d.comp$total_qpcr)
cor1.2 <- cor(d.comp$eng_qmiseq, d.comp$eng_qpcr)
cor1.3 <- cor(d.comp$tra_qmiseq, d.comp$tra_qpcr)

# reads v.s. qPCR, all values
cor2.1 <- cor(d.comp$total_reads, d.comp$total_qpcr)
cor2.2 <- cor(d.comp$eng_reads, d.comp$eng_qpcr)
cor2.3 <- cor(d.comp$tra_reads, d.comp$tra_qpcr)

# qMiSeq v.s. qPCR, exclude outliers
# (apparent outliers were visually judged and excluded)
d.comp.e1 <- d.comp[-c(20, 33),] # exclude outliers
d.comp.e1.0 <- d.comp.e1[d.comp.e1$slope < 10,]
d.comp.e1 <- d.comp.e1[d.comp.e1$slope >= 10,] # include samples without PCR inhibition
cor3.1 <- cor(d.comp.e1$total_qmiseq, d.comp.e1$total_qpcr)

d.comp.e2 <- d.comp[-c(33),] # exclude outliers
d.comp.e2.0 <- d.comp.e2[d.comp.e2$slope < 10,]
d.comp.e2 <- d.comp.e2[d.comp.e2$slope >= 10,] # include samples without PCR inhibition
cor3.2 <- cor(d.comp.e2$eng_qmiseq, d.comp.e2$eng_qpcr)

d.comp.e3 <- d.comp[-c(20, 46),] # exclude outliers
d.comp.e3.0 <- d.comp.e3[d.comp.e3$slope < 10,] # include samples without PCR inhibition
d.comp.e3 <- d.comp.e3[d.comp.e3$slope >= 10,] # include samples without PCR inhibition
cor3.3 <- cor(d.comp.e3$tra_qmiseq, d.comp.e3$tra_qpcr)

# reads v.s. qPCR, exclude outliers
# (the same samples were excluded as outliers for fair comparisons)
cor4.1 <- cor(d.comp.e1$total_reads, d.comp.e1$total_qpcr)
cor4.2 <- cor(d.comp.e2$eng_reads, d.comp.e2$eng_qpcr)
cor4.3 <- cor(d.comp.e3$tra_reads, d.comp.e3$tra_qpcr)

# preparations for visualizations
lm1.1 <- lm(d.comp$total_qpcr ~ d.comp$total_qmiseq)
lm1.2 <- lm(d.comp$eng_qpcr ~ d.comp$eng_qmiseq)
lm1.3 <- lm(d.comp$tra_qpcr ~ d.comp$tra_qmiseq)

lm2.1 <- lm(d.comp$total_qpcr ~ d.comp$total_reads)
lm2.2 <- lm(d.comp$eng_qpcr ~ d.comp$eng_reads)
lm2.3 <- lm(d.comp$tra_qpcr ~ d.comp$tra_reads)

lm3.1 <- lm(d.comp.e1$total_qpcr ~ d.comp.e1$total_qmiseq)
lm3.2 <- lm(d.comp.e2$eng_qpcr ~ d.comp.e2$eng_qmiseq)
lm3.3 <- lm(d.comp.e3$tra_qpcr ~ d.comp.e3$tra_qmiseq)

lm4.1 <- lm(d.comp.e1$total_qpcr ~ d.comp.e1$total_reads)
lm4.2 <- lm(d.comp.e2$eng_qpcr ~ d.comp.e2$eng_reads)
lm4.3 <- lm(d.comp.e3$tra_qpcr ~ d.comp.e3$tra_reads)

# save workspace
save.image(ws.out.03)
