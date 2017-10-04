####
#### R code for Ushio et al.
#### "Quantitative monitoring of multispecies fish environmental DNA using high-throughput sequencing"
#### No.S1 Helper functions
####

# for ggplot2
PlotStyle <-  function(ggobject){
  return(ggobject + theme_bw() + theme(axis.text.x = element_text(angle=0),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.text = element_text(size=12),
                                       axis.title = element_text(size=12),
                                       panel.background=element_rect(colour="black", fill=NA, size=0.8)))
  #return(ggobject)
}

# for conversion from sequence reads to copy numbers
ConvertReads <- function(sample.row = 1){
  x <- std.all[sample.row, 4:8]
  y <- d.reads.all[sample.row,]
  slope <- lm(as.numeric(x) ~ std.copy.1ul + 0)$coefficients
  calc.copy <- y/slope
  return(calc.copy)
}

# for scattered plots
PlotStyle2 <-  function(ggobject, lm.result){
  return(ggobject + geom_point(size = 2) +
           scale_colour_gradient2(low = "black", mid = "red3", high = "red3", midpoint = 30) +
           xlab(expression(paste("qMiSeq (copies ", mu, l^-1, ")"))) +
           ylab(expression(paste("qPCR (copies ", mu, l^-1, ")"))) +
           geom_abline(intercept = lm.result$coefficients[1], slope = lm.result$coefficients[2])
  )
}

# for supplementary scattered plots
PlotStyle3 <-  function(ggobject, lm.result){
  return(ggobject + geom_point(size = 2) +
           scale_colour_gradient2(low = "black", mid = "red3", high = "red3", midpoint = 30) +
           xlab("Sequence reads") +
           ylab(expression(paste("qPCR (copies ", mu, l^-1, ")"))) +
           geom_abline(intercept = lm.result$coefficients[1], slope = lm.result$coefficients[2])
  )
}

# for supplementary scattered plots
PlotStyle4 <-  function(ggobject){
  return(ggobject + geom_point(size = 2) +
           scale_colour_gradient2(low = "black", mid = "red3", high = "red3", midpoint = 30) +
           xlab(expression(paste("Log(qMiSeq [copies ", mu, l^-1, "])"))) +
           ylab(expression(paste("Log(qPCR [copies ", mu, l^-1, "])")))
  )
}
