##---------------------------------------------------------------------#
##
## Title: plot signal intensities coloured by if they were positioned on the
##        outside edge of the plate or not
##
##                    
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#


# The lid lifted off the plate during the lab process. 
# We would hypothesise that this would be more likely
# to affect outer samples.


#----------------------------------------------------------------------#
# LOAD PACKAGES and DATA
#----------------------------------------------------------------------#

library(ggplot2)

load(file = file.path(QCDir, "QCmetrics.rdat"))



#----------------------------------------------------------------------#
# CREATE PLOT
#----------------------------------------------------------------------#

# create object with sample positions for outer samples to colour the plot by
outer <- unique(c(paste0(toupper(letters[1:8]), 1),
           paste0(toupper(letters[1:8]), 6),
           paste0("A", 1:12),
           paste0("H", 1:12)))

QCmetrics$platePosition <- rep("Inner")
QCmetrics$platePosition[QCmetrics$Location %in% outer] <- "Outer" 


ggplot(QCmetrics, aes(x = M.median, y = U.median, colour = platePosition)) +
  geom_point() +
  xlab("Median M intensity") +
  ylab("Median U intensity") +
  ggtitle("Signal Intensities by Sample Position")
