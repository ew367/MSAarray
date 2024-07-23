##---------------------------------------------------------------------#
##
## Title: correlate MSA test data to data on same samples ran on
##        450k array
##
##                    
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# 11 sammples on MSA array don't appear to be in UCL dataset - looking
# at pool ID in SCZ and Individual_ID in MSA

# [1] "C2436"  "00S127" "0S2282" "C2434"  "C2478"  "0S2297" "00S222" "00S069" "C2549"  "C2324"  "C0254" 


# there are 121277 probes on both arrays

#----------------------------------------------------------------------#
# LOAD PACKAGES and DATA
#----------------------------------------------------------------------#


library(ENmix)
library(ggplot2)
library(dplyr)
library(minfi)

source("config.r")


# load manifest
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)


# get 450k probes
msa450k <- man[!(man$Methyl450_Locus_Match == ""), c("IlmnID", "Name", "Methyl450_Locus_Match")] #

# load MSA data
sampleSheet <- read.csv(pheno, stringsAsFactors = F)
sampleSheet <- sampleSheet[sampleSheet$Study == "UCL SZ",]

load(file = file.path(QCDir, "mraw.rdat"))
MSAbetas <- as.data.frame(getB(mraw))


#----------------------------------------------------------------------#
#  Import UCL SCZ (450k) data
#----------------------------------------------------------------------#

load("/lustre/projects/Research_Project-MRC190311/DNAm/mrcSCZBlood/UCL/2_normalised/normalised.rdata")
# betas.ucl and pheno.ucl

# subset data to samples there is matched data for
pheno <- pheno.ucl[pheno.ucl$Pool_ID %in% sampleSheet$Individual_ID,]
pheno <- pheno[order(pheno$Pool_ID),]

betas <- as.data.frame(betas.ucl[,pheno$Basename])
betas <- betas[msa450k$Methyl450_Locus_Match,]

# match back as there are some missing (see notes)
sampleSheet <- sampleSheet[sampleSheet$Individual_ID %in% pheno$Pool_ID,]
sampleSheet <- sampleSheet[order(sampleSheet$Individual_ID),]
MSAbetas <- MSAbetas[,sampleSheet$Basename]

# check that matched samples are in the same order
identical(sampleSheet$Individual_ID, pheno$Pool_ID)
# TRUE


#----------------------------------------------------------------------#
# Join betas from 450k and MSA data into single object
#----------------------------------------------------------------------#

# join betas from epicv1 and MSA
plotBetas <- betas
plotBetas$Methyl450_Locus_Match <- row.names(plotBetas)
plotBetas <- left_join(plotBetas, msa450k)
MSAbetas$IlmnID <- row.names(MSAbetas)
plotBetas <- left_join(plotBetas, MSAbetas)
plotBetas <- plotBetas[complete.cases(plotBetas),]

# create data frame for colouring distribution plot by array version
distGroup <- bind_rows(
  sampleSheet %>% 
    dplyr::select(Basename) %>%
    mutate(tech = "MSA"),
  
  pheno %>%
    dplyr::select(Basename) %>%
    mutate(tech = "450k"))


#----------------------------------------------------------------------#
# Plot cor stat accross all samples
#----------------------------------------------------------------------#

# define function to apply to all probes
probeWiseCor <- function(x){
  corStat <- cor.test(as.numeric(x[pheno$Basename]), 
                      as.numeric(x[sampleSheet$Basename]))$estimate
  return(corStat)
}


# apply to each row (probe) of the matrix
start.time <- Sys.time()
corStat <- apply(plotBetas, 1, probeWiseCor)
end.time <- Sys.time()
round(end.time - start.time,2) #1.94 mins


## plot the distribution
corStat <- as.data.frame(corStat)

ggplot(corStat, aes(x=corStat))+
  geom_density()+
  ggtitle("Probewise R values for 450k vs MSA array probes")


#----------------------------------------------------------------------#
# Plot density, correlations and violin plot
#----------------------------------------------------------------------#


pdf("plots/450KvsMSAarray.pdf")

# plot distribution coloured by MSA/V1

plotBetas <- plotBetas[,distGroup$Basename]
ArrayType <- distGroup$tech
  
# denisty plot
minfi::densityPlot(as.matrix(plotBetas), sampGroups = ArrayType, main = "Whole Blood")
  
  
  # get mean methylation for UCL SZ (450k) data
  meanBetas <- as.data.frame(rowMeans(betas[,pheno$Basename]))
  colnames(meanBetas) <- "mean450k"
  meanBetas$Methyl450_Locus_Match <- row.names(meanBetas)
  
  # get mean methylation for MSA data
  meanMSA <- as.data.frame(rowMeans(MSAbetas[,sampleSheet$Basename]))
  colnames(meanMSA) <- "meanMSA"
  meanMSA$IlmnID <- row.names(meanMSA)
  
  # join together
  corPlotDF <- left_join(meanBetas, msa450k)
  corPlotDF <- left_join(corPlotDF, meanMSA)
  corPlotDF <- corPlotDF[complete.cases(corPlotDF),]
  
  p <- ggplot(corPlotDF, aes(mean450k, meanMSA))+
    geom_point()+
    ggtitle(paste0("Whole Blood: R = ", signif(cor.test(corPlotDF$mean450k, corPlotDF$meanMSA)$estimate, 2)))+
    geom_abline(colour="red", linetype="dashed")
  
  print(p)
  
  
# violin plot
  
violinDF <- melt(corPlotDF)
  
p2 <- ggplot(violinDF, aes(x=variable, y=value))+
    geom_violin()+
  xlab("Array")+
  ylab("Methylation")+
  ggtitle("Whole Blood")
  
  print(p2)


dev.off()

