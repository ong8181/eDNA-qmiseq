####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.2 The conversion of sequence reads to calculated copy numbers
####

# names for the output folder and workspace file
qmiseq.out.01 <- "01_QMiSeq_StdEvalOut"
qmiseq.out.02 <- "02_QMiSeq_ConversionOut"
dir.create(qmiseq.out.02, showWarnings = FALSE)
ws.out.01 <- file.path(qmiseq.out.01, "01_QMiSeq_StdEvalOut.RData")
ws.out.02 <- file.path(qmiseq.out.02, "02_QMiSeq_ConversionOut.RData")

# load workspace
load(ws.out.01)

# load functions
source("functions/HelperFuncs.R")

# convert reads 
converted.reads <- data.frame()
for(i in 1:nrow(d.sample)){
  converted.reads <- rbind(converted.reads, ConvertReads(sample.row = i))
}

# replace NaN with NA
converted.reads[is.na(converted.reads)] <- NA

# exclude negative control samples
converted.reads.sample <- converted.reads[d.sample$treatment == "sample",]
std.name.cols <- match(std.name, colnames(converted.reads.sample))
converted.reads.sample <- converted.reads.sample[,-std.name.cols]

# file name for converted reads
converted.name1 <- file.path(qmiseq.out.02, "converted_reads_all.csv")
converted.name2 <- file.path(qmiseq.out.02, "converted_reads_sample.csv")

write.csv(converted.reads, converted.name1, row.names = F)
write.csv(converted.reads.sample, converted.name2, row.names = F)

# save workspace
save.image(ws.out.02)
