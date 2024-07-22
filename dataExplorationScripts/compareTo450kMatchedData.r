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
load(file = file.path(QCDir, "mraw.rdat"))
MSAbetas <- as.data.frame(getB(mraw))
sampleSheet <- read.csv(pheno, stringsAsFactors = F)


# load SCZ blood data
load("/lustre/projects/Research_Project-MRC190311/DNAm/mrcSCZBlood/UCL/2_normalised/normalised.rdata")
# betas.ucl and pheno.ucl


#----------------------------------------------------------------------#
# UCL SCZ (450k)
#----------------------------------------------------------------------#

# subset data to samples there is matched data for
pheno <- pheno.ucl[pheno.ucl$Pool_ID %in% sampleSheet$Individual_ID,]

betas <- as.data.frame(betas.ucl[,pheno$Basename])
betas <- betas[msa450k$Methyl450_Locus_Match,]

sampleSheet <- sampleSheet[sampleSheet$Individual_ID %in% pheno$Pool_ID,]
MSAbetas <- MSAbetas[,sampleSheet$Basename]


# create data frame for colouring distribution plot by array version
distGroup <- bind_rows(
  sampleSheet %>% 
    dplyr::select(Basename) %>%
    mutate(tech = "MSA"),
  
  pheno %>%
    dplyr::select(Basename) %>%
    mutate(tech = "450k"))

# join betas from epicv1 and MSA
plotBetas <- betas
plotBetas$Methyl450_Locus_Match <- row.names(plotBetas)
plotBetas <- left_join(plotBetas, msa450k)
MSAbetas$IlmnID <- row.names(MSAbetas)
plotBetas <- left_join(plotBetas, MSAbetas)
plotBetas <- plotBetas[complete.cases(plotBetas),]


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








