##---------------------------------------------------------------------#
##
## Title: Explore bisulphite conversion probes 
##
## Purpose of script: look more closely at bisulphite conversion probes
##                    to see which to include in QCmetrics
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

#----------------------------------------------------------------------#
# LOAD PACKAGES etc.
#----------------------------------------------------------------------#

library("dplyr")
library("ggplot2")

source("config.r")

#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#


load(file = file.path(QCDir, "bsCon.rdat"))
BSconAll <- as.data.frame(BSconAll)
BSconAll$Address <- row.names(BSconAll)


load(file = file.path(QCDir, "rgSet.rdat"))


#----------------------------------------------------------------------#
# EXTRACT CTRL PROBE DATA
#----------------------------------------------------------------------#


ctrls <- metadata(rgSet)$ictrl
bs.both <- as.data.frame(rbind(ctrls[ctrls$Type == "BISULFITE CONVERSION I",], ctrls[ctrls$Type == "BISULFITE CONVERSION II",]))
bs.both$CorU <- rep("typeII")
bs.both$CorU[1:10] <- substr(bs.both$ExtendedType, nchar(bs.both$ExtendedType)-1, nchar(bs.both$ExtendedType)-1)[1:10] # only type 1 probes have C or U


#----------------------------------------------------------------------#
# PLOT BS PROBES
#----------------------------------------------------------------------#

# create object containing all info needed for plot
BSconAll <- left_join(BSconAll, as.data.frame(bs.both %>% dplyr::select(Address, CorU)))
plotdf <- reshape2::melt(BSconAll)

# plot
ggplot(plotdf, aes(x=probe, y=value, fill=CorU))+
  geom_boxplot()+
  ylab("Methylation")+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Bisulphite Conversion Probes on MSA Array")



#----------------------------------------------------------------------#
# CALCULATE BSCON EFFICIENCY EXCLUDING U PROBES
#----------------------------------------------------------------------#


keepProbes <- bs.both$Address[bs.both$CorU != "U"] # get probes which are not U

# get bscon probes
bs.type1 <- ctrls$Address[ctrls$Type == "BISULFITE CONVERSION I"]
bs.type2 <- ctrls$Address[ctrls$Type == "BISULFITE CONVERSION II"]

bs.green.type1 <- assays(rgSet)$Green[bs.type1,]
bs.green.type2 <- assays(rgSet)$Green[bs.type2,]
bs.red.type1 <- assays(rgSet)$Red[bs.type1,]
bs.red.type2 <- assays(rgSet)$Red[bs.type2,]

red.med <- apply(bs.red.type1, 1, median)
green.med <- apply(bs.green.type1, 1, median)

redGreater<-rowMeans(bs.red.type1) > rowMeans(bs.green.type1)

BScon1<-rbind(
  bs.red.type1[which(redGreater),] / ( bs.red.type1[which(redGreater),] + bs.green.type1[which(redGreater),] ), 
  bs.green.type1[which(!redGreater),] / ( bs.red.type1[which(!redGreater),] + bs.green.type1[which(!redGreater),] )
)

BScon2 <- bs.red.type2 / ( bs.red.type2 + bs.green.type2 )

BSconAll<-rbind(BScon1, BScon2)


# remove U probes

BSconAll <- BSconAll[keepProbes,]

BScon.med<-apply(BSconAll, 2, median)
BSconAll<-as.data.frame(rbind(BSconAll, BScon.med)*100)



#----------------------------------------------------------------------#
# PLOT NEW  HISTOGRAM
#----------------------------------------------------------------------#

histDF <- as.data.frame(t(BSconAll))

ggplot(histDF, aes(x = BScon.med))+
  geom_histogram()+
  geom_vline(xintercept=bsConThresh, colour = "red", linetype="dashed")+
  geom_vline(xintercept=80, colour = "red", linetype="dashed")+
  
  ggtitle("BSconversion exclusing U probes")
S







