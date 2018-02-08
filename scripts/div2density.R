library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

setwd(args[1])

divergencesFl <- read.csv(args[2], header=TRUE, sep='\t')
label <- c(rep('Uncorrected', nrow(divergencesFl)), rep('Corrected', nrow(divergencesFl)))
values <- c(divergencesFl$divergence, divergencesFl$correctedDivergence)
divergences <- data.frame(divergences=values, type=label)
most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence))




pdf(file=paste(args[1], '/dC_vs_d.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergences, color=label, fill=label)) +
  geom_density(alpha=0.5) +
  theme_minimal() +
    xlab("LTR divergence") +
  theme(legend.justification=c('right', 'top'),
        legend.key.size=unit(40,"pt"),
        legend.position=c(0.76, .76),
        legend.text=element_text(size=20, face="italic"),
        legend.title = element_blank(),
        axis.text.x = element_text(size=16, face="bold"),
        axis.text.y = element_text(size=16, face="bold"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16)) + 
  	scale_x_continuous(limits=c(0, most))

dev.off()



divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$divergence)
most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence))

pdf(file=paste(args[1], '/divByClass_uncorrected.density.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergence, color=classification, fill=classification)) +
  geom_density(alpha=0.5) +
  theme_minimal() +
    xlab("LTR divergence") +
  theme(legend.justification=c('right', 'top'),
        legend.key.size=unit(40,"pt"),
        legend.position=c(0.76, .76),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        axis.text.x = element_text(size=16, face="bold"),
        axis.text.y = element_text(size=16, face="bold"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16)) + 
  	scale_x_continuous(limits=c(0, most))


dev.off()

pdf(file=paste(args[1], '/divByClass_uncorrected.hist.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergence, color=classification, fill=classification)) +
    geom_histogram(bins = 100, alpha=0.5) +
    theme_minimal() + 
    xlab("LTR divergence") +
  
    theme(legend.justification=c('right', 'top'),
          legend.key.size=unit(40,"pt"),
        legend.position=c(0.76, .76),
                  #legend.position="none",
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          axis.text.x = element_text(size=16, face="bold"),
          axis.text.y = element_text(size=16, face="bold"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16)) + 
  	scale_x_continuous(limits=c(0, most))


dev.off()



divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$correctedDivergence)
most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence))

pdf(file=paste(args[1], '/divByClass_geneconvCorrected.density.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergence, color=classification, fill=classification)) +
  geom_density(alpha=0.5) +
  theme_minimal() +
    xlab("LTR divergence  (corrected for gene conversion)") +
#  ylab("density") +
  theme(legend.justification=c('right', 'top'),
        legend.key.size=unit(40,"pt"),
        legend.position=c(0.76, .76),
        legend.text=element_text(size=20, face="italic"),
        legend.title = element_blank(),
        axis.text.x = element_text(size=16, face="bold"),
        axis.text.y = element_text(size=16, face="bold"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16)) + 
  	scale_x_continuous(limits=c(0, most))


dev.off()

pdf(file=paste(args[1], '/divByClass_geneconvCorrected.hist.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergence, color=classification, fill=classification)) +
    geom_histogram(bins = 100, alpha=0.5) +
    theme_minimal() + 
    xlab("LTR divergence  (corrected for gene conversion)") +
  
    theme(legend.justification=c('right', 'top'),
          legend.key.size=unit(40,"pt"),
        legend.position=c(0.76, .76),
                  #legend.position="none",
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          axis.text.x = element_text(size=16, face="bold"),
          axis.text.y = element_text(size=16, face="bold"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16)) + 
  	scale_x_continuous(limits=c(0, most))


dev.off()

