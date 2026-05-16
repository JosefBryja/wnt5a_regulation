#================================================
# Visualize TFs motifs on wnt5a promotor region
#================================================

###################################
### Libraries and data loading ####
###################################

# Libraries
library(ggplot2)

# Load the data
motifs <- read.csv("Analysis/high_confidence_hits.txt", sep = "\t")


###################################################
### Visualization of localization of homer hits ###
###################################################

ggplot(motifs, aes(x = Offset)) +
  geom_point(aes(y = 1, col = MotifScore), size = 8)
