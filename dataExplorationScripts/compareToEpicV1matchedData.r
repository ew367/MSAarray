##---------------------------------------------------------------------#
##
##    Title:        Correlate MSA test data to epic V1 data
##
##    This script:  plots distribution of probe level cor statistics  
##                  creates density plots (one for each cell type)
##                  correlation plots for each cell type
##                  violin plot for each cell type
##                  
##                  
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

#https://emea.illumina.com/products/by-type/microarray-kits/infinium-mouse-methylation.html

# there are 278 probes which are in both MSA array and epic V1 but have
# different cg IDs

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

# get epicv1 probes
msaV1 <- man[!(man$EPICv1_Locus_Match == ""), c("IlmnID", "Name", "EPICv1_Locus_Match")] #

# load MSA data and subset to epicv1 data 
sampleSheet <- read.csv(pheno, stringsAsFactors = F)
sampleSheet <- sampleSheet[sampleSheet$Study == "A-Risk",]
sampleSheet$Sample_ID <- gsub(" ", "_", sampleSheet$Sample_ID)
sampleSheet <- sampleSheet[order(sampleSheet$Sample_ID),]


load(file = file.path(QCDir, "mraw.rdat"))
MSAbetas <- as.data.frame(getB(mraw))
MSAbetas <- MSAbetas[,sampleSheet$Basename]



#----------------------------------------------------------------------#
#  Import A-Risk (epic v1) data
#----------------------------------------------------------------------#


# load A risk data
load("/lustre/projects/Research_Project-MRC190311/DNAm/Arisk/2_processed/CellSortedAsthmaERisk_Normalised.rdat")
#betas and pheno

# subset Arisk data to samples there is MSA samples/probe data for
pheno <- pheno[pheno$PatientID %in% sampleSheet$Individual_ID,]
pheno <- pheno[-which(pheno$Sample.Type == "Buccal" | pheno$Sample.Type == "Nasal"),] #remove Buccal and Nasal cells
pheno$Sample_ID <- paste0(pheno2$PatientID, "_", pheno$Sample.Type)
pheno$Sample_ID <- gsub(" ", "_", pheno$Sample_ID)
pheno$Sample_ID <- gsub("w", "W", pheno$Sample_ID)
pheno$Sample_ID <- gsub("b", "B", pheno$Sample_ID)
pheno <- pheno[order(pheno$Sample_ID),]

#subset betas to samples and CpgS needed
betas <- as.data.frame(betas[,pheno$Basename])
betas <- betas[msaV1$EPICv1_Locus_Match,]


# check samples from the datasaets are now in the same order
identical(pheno$Sample_ID, sampleSheet$Sample_ID)

#----------------------------------------------------------------------#
# Join betas from epicv1 and MSA data into single object
#----------------------------------------------------------------------#

plotBetas <- betas
plotBetas$EPICv1_Locus_Match <- row.names(plotBetas)
plotBetas <- left_join(plotBetas, msaV1)
MSAbetas$IlmnID <- row.names(MSAbetas)
plotBetas <- left_join(plotBetas, MSAbetas)
plotBetas <- plotBetas[complete.cases(plotBetas),]


# create data frame for colouring distribution plot by array version
distGroup <- bind_rows(
  sampleSheet %>% 
    dplyr::select(Basename, Cell_Type) %>%
    mutate(tech = "MSA"),
  
  pheno %>%
    dplyr::select(Basename, Sample.Type) %>%
    mutate(tech = "epicV1") %>%
    rename(Cell_Type = Sample.Type))


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
  ggtitle("Probewise R values for EpicV1 vs MSA array probes")



#----------------------------------------------------------------------#
# Within Cell Type Plots
#----------------------------------------------------------------------#

#  B-cells  CD4 T-cells  CD8 T-cells Granulocytes Monocytes, whole blood

# for loop for each cell type: density plot, correlation plot

pdf("plots/EpicV1vsMSAarrayByCellType.pdf")

for(cell in unique(sampleSheet$Cell_Type)){

  # plot distribution, each cell type on seperate plot, coloured by MSA/V1

  # subset to just cell type betas
  cellDistBetas <- plotBetas[,distGroup$Basename[distGroup$Cell_Type == cell]]
  ArrayType <- distGroup$tech[distGroup$Cell_Type == cell]

  # denisty plot
  minfi::densityPlot(as.matrix(cellDistBetas), sampGroups = ArrayType, main=cell)
  


  # get mean methylation for the cell type from A risk data
  cellBetas <- as.data.frame(rowMeans(betas[,pheno$Basename[pheno$Sample.Type == cell]]))
  colnames(cellBetas) <- "meanEpicV1"
  cellBetas$EPICv1_Locus_Match <- row.names(cellBetas)

  # get mean methylation for the cell type from MSA data
  meanMSA <- as.data.frame(rowMeans(MSAbetas[,sampleSheet$Basename[sampleSheet$Cell_Type == cell]]))
  colnames(meanMSA) <- "meanMSA"
  meanMSA$IlmnID <- row.names(meanMSA)

  # join together
  corPlotDF <- left_join(cellBetas, msaV1)
  corPlotDF <- left_join(corPlotDF, meanMSA)
  corPlotDF <- corPlotDF[complete.cases(corPlotDF),]

  p <- ggplot(corPlotDF, aes(meanEpicV1, meanMSA))+
    geom_point()+
    ggtitle(paste0(cell, " - R = ", signif(cor.test(corPlotDF$meanEpicV1, corPlotDF$meanMSA)$estimate, 2)))+
    geom_abline(colour="red", linetype="dashed")
  
  print(p)
  
}

dev.off()

#----------------------------------------------------------------------#
# Violin plot - all cell types on one 
#----------------------------------------------------------------------#

# need a dataframe with cols for methylation value, cell type, array

t <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(t)= c("Methylation", "array", "cell")

for(cell in unique(sampleSheet$Cell_Type)){
  
  # get mean methylation for the cell type from A risk data
  cellBetas <- as.data.frame(rowMeans(betas[,pheno$Basename[pheno$Sample.Type == cell]]))
  colnames(cellBetas) <- "Methylation"
  cellBetas$array <- "epicV1"
  cellBetas$cell <- cell
  
# get mean methylation for the cell type from MSA data
  meanMSA <- as.data.frame(rowMeans(MSAbetas[,sampleSheet$Basename[sampleSheet$Cell_Type == cell]]))
  colnames(meanMSA) <- "Methylation"
  meanMSA$array <- "MSA"
  meanMSA$cell <- cell
  
  # join together
  
  all <- rbind(cellBetas, meanMSA)
  t <- rbind(t, all)

}


pdf("plots/EpicV1vsMSAarrayViolinPlot.pdf")
p2 <- ggplot(t, aes(x=cell, y=Methylation, colour=array))+
  geom_violin()+
  #xlab("Array")+
  ylab("Methylation")

print(p2)


dev.off()




  
