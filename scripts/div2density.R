library(ggplot2)
library(plyr) # For count()

# Inputs:
# --------
# Arg1: wd (for output)
# Arg2: LTR divergences


args <- commandArgs(trailingOnly=TRUE)

# Output files go here
setwd(args[1]) # UNCOMMENT
# setwd('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir') # COMMENT
# Second argument is the LTR divergence summary tab file.
divergencesFl <- read.csv(args[2], header=TRUE, sep='\t') # UNCOMMENT
# divergencesFl <- read.csv('LTR_divergences.tab', header=TRUE, sep='\t') # COMMENT


#--------------------------------------------------------------------------------------------
# All elements' divergence estimates: 1 density plot, overlayed corrected vs. uncorrected
#--------------------------------------------------------------------------------------------

label <- c(rep('Uncorrected', nrow(divergencesFl)), rep('Corrected', nrow(divergencesFl)))
values <- c(divergencesFl$divergence, divergencesFl$correctedDivergence)
divergences <- data.frame(divergences=values, type=label)
most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence)) # Get largest divergence value for setting the x axis range in the plots

# First plot: divergences for all elements
pdf(file=paste(args[1], '/All.raw_vs_corrected_divergences.density.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergences, color=label, fill=label)) +
  geom_density(alpha=0.5) +
  theme_minimal() +
  
    xlab("LTR divergence") +
  theme(panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        legend.justification=c('right', 'top'),
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


#--------------------------------------------------------------------------------------------
# Histogram, density, and violin plots for each classif, for corrected and uncorrected divergences
#--------------------------------------------------------------------------------------------
# All superfamilies on same plot: density plot

divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$divergence)
most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence)) # Get largest divergence value for setting the x axis range in the plots

pdf(file=paste(args[1], '/Superfamily.raw.density.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergence, color=classification, fill=classification)) +
  geom_density(alpha=0.5) +
  theme_minimal() +
    xlab("LTR divergence") +
  theme(panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        legend.justification=c('right', 'top'),
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


#--------------------------------------------------------------------------------------------
# All superfamilies on same plot: histogram (uncorrected)

pdf(file=paste(args[1], '/Superfamily.raw_divergences.hist.pdf', sep=''), height=5)

ggplot(divergences, aes(x=divergence, color=classification, fill=classification)) +
    geom_histogram(bins = 100, alpha=0.5) +
    theme_minimal() + 
    xlab("LTR divergence") +
  
    theme(panel.grid.minor = element_blank(),  
          panel.grid.major = element_blank(),
          legend.justification=c('right', 'top'),
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


#--------------------------------------------------------------------------------------------
# Each superfamily on different plot: density plot (uncorrected)
superfamilies <- as.character(unique(divergences$classification))
divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$divergence)
most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence))
# superfam_ct <- count(divergences$classification)


for (i in 1:length(superfamilies)){
  current_recs <- divergences[which(divergences$classification == superfamilies[i]),]
  pdf(file=paste(args[1], '/', superfamilies[i], '.density.pdf', sep=''), height=5)
  
  # pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '.raw_divergences.density.pdf', sep=''), height=5)

  p <- ggplot(current_recs, aes(x=divergence)) +
    geom_density(alpha=0.5) +
    theme_minimal() +
      xlab("LTR divergence (raw)") +
  #  ylab("density") +
    theme(panel.grid.minor = element_blank(),  
          panel.grid.major = element_blank(),
          legend.justification=c('right', 'top'),
          legend.key.size=unit(40,"pt"),
          legend.position="none",
          legend.text=element_text(size=20, face="italic"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=16, face="bold"),
          axis.text.y = element_text(size=16, face="bold"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16)) + 
    	scale_x_continuous(limits=c(0, most))
  
  print(p)
  
  dev.off()

}
#--------------------------------------------------------------------------------------------
# Each superfamily on different plot: histogram (uncorrected)
for (i in 1:length(superfamilies)){
  current_recs <- divergences[which(divergences$classification == superfamilies[i]),]
  pdf(file=paste(args[1], '/', superfamilies[i], '.histogram.pdf', sep=''), height=5)
  
  # pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '.raw_divergences.histogram.pdf', sep=''), height=5)

  p<- ggplot(current_recs, aes(x=divergence)) +
    geom_histogram(bins = 100, alpha=0.5) +
    theme_minimal() + 
    xlab("LTR divergence (raw") +
  
    theme(panel.grid.minor = element_blank(),  
          panel.grid.major = element_blank(),
          legend.justification=c('right', 'top'),
          legend.key.size=unit(40,"pt"),
          legend.position="non",
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          axis.text.x = element_text(size=16, face="bold"),
          axis.text.y = element_text(size=16, face="bold"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16)) + 
  	scale_x_continuous(limits=c(0, most))

  print(p)

  dev.off()
}



#--------------------------------------------------------------------------------------------
# Each superfamily on different plot: density plot (corrected)
superfamilies <- as.character(unique(divergences$classification))
divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$correctedDivergence)
most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence))
# superfam_ct <- count(divergences$classification)


for (i in 1:length(superfamilies)){
  current_recs <- divergences[which(divergences$classification == superfamilies[i]),]
  pdf(file=paste(args[1], '/', superfamilies[i], '.density.pdf', sep=''), height=5)
  
  # pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '.corrected_divergences.density.pdf', sep=''), height=5)

  p <- ggplot(current_recs, aes(x=divergence)) +
    geom_density(alpha=0.5) +
    theme_minimal() +
      xlab("LTR divergence (raw)") +
  #  ylab("density") +
    theme(panel.grid.minor = element_blank(),  
          panel.grid.major = element_blank(),
          legend.justification=c('right', 'top'),
          legend.key.size=unit(40,"pt"),
          legend.position="none",
          legend.text=element_text(size=20, face="italic"),
          legend.title = element_blank(),
          axis.text.x = element_text(size=16, face="bold"),
          axis.text.y = element_text(size=16, face="bold"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16)) + 
    	scale_x_continuous(limits=c(0, most))
  
  print(p)
  
  dev.off()

}
#--------------------------------------------------------------------------------------------
# Each superfamily on different plot: histogram (corrected)
for (i in 1:length(superfamilies)){
  current_recs <- divergences[which(divergences$classification == superfamilies[i]),]
  pdf(file=paste(args[1], '/', superfamilies[i], '.histogram.pdf', sep=''), height=5)
  
  # pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '.corrected_divergences.histogram.pdf', sep=''), height=5)

  p<- ggplot(current_recs, aes(x=divergence)) +
    geom_histogram(bins = 100, alpha=0.5) +
    theme_minimal() + 
    xlab("LTR divergence (raw") +
  
    theme(panel.grid.minor = element_blank(),  
          panel.grid.major = element_blank(),
          legend.justification=c('right', 'top'),
          legend.key.size=unit(40,"pt"),
          legend.position="non",
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          axis.text.x = element_text(size=16, face="bold"),
          axis.text.y = element_text(size=16, face="bold"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16)) + 
  	scale_x_continuous(limits=c(0, most))

  print(p)

  dev.off()
}



#--------------------------------------------------------------------------------------------
# Box/Violins, by-superfamily & by-cluster
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# All superfamilies on the same plot: violin plots (uncorrected)
superfam_ct <- count(divergences$classification)
divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$divergence)

most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence)) # Get largest divergence value for setting the x axis range in the plots

pdf(file=paste(args[1], '/Superfamilies.raw_divergences_boxplots.pdf', sep=''), height=4, width=16)
 
# pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '_LTRdiv_boxplots.pdf', sep=''), height=4, width=16)


ggplot(divergences, aes(x=classification, y=divergence)) +
  theme_bw() +
  theme(legend.position='none', 
        #panel.grid.minor = element_blank(),  
        #panel.grid.major = element_blank(), 
        plot.title=element_blank(), 
        axis.text.y=element_text(size=16, face="bold"), 
        #axis.title.y=element_text(size=18, face="bold"), 
        axis.text.x=element_text(size=16, face="bold")) +
        theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5)) +
  xlab("") +
  ylim(c(0,1)) +
  ylab("LTR divergence") +
  geom_violin(scale='count') +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.9)
  # geom_boxplot()
  # geom_violin(scale='count', draw_quantiles=c(0.25,0.5,0.75))
  #scale_y_log10(breaks=c(10, 100, 1000, 10000, 100000, 1000000))


#--------------------------------------------------------------------------------------------
# All superfamilies on the same plot: violin plots (corrected)
superfam_ct <- count(divergences$classification)
divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$correctedDivergence)

most <- max(c(divergencesFl$divergence, divergencesFl$correctedDivergence)) # Get largest divergence value for setting the x axis range in the plots

pdf(file=paste(args[1], '/Superfamilies.corrected_divergences.boxplots.pdf', sep=''), height=4, width=16)
 
# pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '_LTRdiv_boxplots.pdf', sep=''), height=4, width=16)


ggplot(divergences, aes(x=classification, y=divergence)) +
  theme_bw() +
  theme(legend.position='none', 
        #panel.grid.minor = element_blank(),  
        #panel.grid.major = element_blank(), 
        plot.title=element_blank(), 
        axis.text.y=element_text(size=16, face="bold"), 
        #axis.title.y=element_text(size=18, face="bold"), 
        axis.text.x=element_text(size=16, face="bold")) +
        theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5)) +
  xlab("") +
  ylim(c(0,1)) +
  ylab("LTR divergence") +
  geom_violin(scale='count') +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.9)
  # geom_boxplot()
  # geom_violin(scale='count', draw_quantiles=c(0.25,0.5,0.75))
  #scale_y_log10(breaks=c(10, 100, 1000, 10000, 100000, 1000000))





#--------------------------------------------------------------------------------------------
# All clusters on the same plot: violin plots (uncorrected) (within each superfamily)

# divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$divergence, cluster=as.factor(divergencesFl$cluster), clusterSize=as.factor(divergencesFl$clusterSize))
divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$divergence, cluster=as.factor(paste(as.character(divergencesFl$cluster), "\n(", as.character(divergencesFl$clusterSize), ")", sep="")))
superfamilies <- as.character(unique(divergences$classification))

for (i in 1:length(superfamilies)){
  current_recs <- divergences[which(divergences$classification == superfamilies[i]),]
  # Remove singleton dopelton clusters
  non_singleton_clusters <- as.character(count(current_recs$cluster)$x[which(count(current_recs$cluster)$freq >2 )])
  current_recs <- current_recs[current_recs$cluster %in% non_singleton_clusters,]
  
  pdf(file=paste(args[1], '/', superfamilies[i], '.raw_divergences.boxplots.pdf', sep=''), height=4, width=16)
 
  # pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '_LTRdiv_boxplots.pdf', sep=''), height=4, width=16)
  
  p <- ggplot(current_recs, aes(x=cluster, y=divergence)) +
    theme_bw() +
    theme(legend.position='none', 
          #panel.grid.minor = element_blank(),  
          #panel.grid.major = element_blank(), 
          plot.title=element_blank(), 
          axis.text.y=element_text(size=16, face="bold"), 
          axis.title.y=element_text(size=14, face="bold"),
          axis.title.x=element_text(size=14, face="bold"),
          axis.text.x=element_text(size=8, face="bold")) +
          theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5)) +
    xlab(superfamilies[i]) +
    ylim(c(0,1)) +
    ylab("LTR divergence") +
    # geom_violin(scale='count', draw_quantiles=c(0.25,0.5,0.75))
    geom_boxplot()
  
  print(p)
  
  dev.off()
}


#--------------------------------------------------------------------------------------------
# All clusters on the same plot: violin plots (corrected) (within each superfamily)

# divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$divergence, cluster=as.factor(divergencesFl$cluster), clusterSize=as.factor(divergencesFl$clusterSize))
divergences <- data.frame(classification=divergencesFl$classification, divergence=divergencesFl$correctedDivergence, cluster=as.factor(paste(as.character(divergencesFl$cluster), "\n(", as.character(divergencesFl$clusterSize), ")", sep="")))
superfamilies <- as.character(unique(divergences$classification))

for (i in 1:length(superfamilies)){
  current_recs <- divergences[which(divergences$classification == superfamilies[i]),]
  # Remove singleton dopelton clusters
  non_singleton_clusters <- as.character(count(current_recs$cluster)$x[which(count(current_recs$cluster)$freq >2 )])
  current_recs <- current_recs[current_recs$cluster %in% non_singleton_clusters,]
  
  pdf(file=paste(args[1], '/', superfamilies[i], '.corrected_divergences.boxplots.pdf', sep=''), height=4, width=16)
 
  # pdf(file=paste('/Users/mathewsimenc/dropbox_csuf/LTRAn/R_scripts/workingdir', '/', superfamilies[i], '_LTRdiv_boxplots.pdf', sep=''), height=4, width=16)
  
  p <- ggplot(current_recs, aes(x=cluster, y=divergence)) +
    theme_bw() +
    theme(legend.position='none', 
          #panel.grid.minor = element_blank(),  
          #panel.grid.major = element_blank(), 
          plot.title=element_blank(), 
          axis.text.y=element_text(size=16, face="bold"), 
          axis.title.y=element_text(size=14, face="bold"),
          axis.title.x=element_text(size=14, face="bold"),
          axis.text.x=element_text(size=8, face="bold")) +
          theme(axis.text.x=element_text(angle=0, vjust=0, hjust=0.5)) +
    xlab(superfamilies[i]) +
    ylim(c(0,1)) +
    ylab("LTR divergence") +
    # geom_violin(scale='count', draw_quantiles=c(0.25,0.5,0.75))
    geom_boxplot()
  
  print(p)
  
  dev.off()
}
