##---------------------------------------------------------------------#
##
##    Title:        Correlate MSA test data to epic V1 data
##
##    This script:  creates density plots (one for each cell type)
##                  correlation plots for each cell type
##                  violin plot for each cell type
##                  plots distribution of probe level cor statistics
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
#sampleSheet <- sampleSheet[sort(sampleSheet$Individual_ID),]

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
pheno2 <- pheno[pheno$PatientID %in% sampleSheet$Individual_ID,]
 #<- pheno[match(pheno$PatientID, sampleSheet$Individual_ID),]
#pheno <- pheno[sort(pheno$PatientID),]
pheno2 <- pheno2[-which(pheno2$Sample.Type == "Buccal" | pheno2$Sample.Type == "Nasal"),] #remove Buccal and Nasal cellsm

betas <- as.data.frame(betas[,pheno$Basename])
betas <- betas[msaV1$EPICv1_Locus_Match,]


#----------------------------------------------------------------------#
# Plot cor stat accross all samples
#----------------------------------------------------------------------#

# join betas from epicv1 and MSA
plotBetas <- betas
plotBetas$EPICv1_Locus_Match <- row.names(plotBetas)
plotBetas <- left_join(plotBetas, msaV1)
MSAbetas$IlmnID <- row.names(MSAbetas)
plotBetas <- left_join(plotBetas, MSAbetas)
plotBetas <- plotBetas[complete.cases(plotBetas),]


#define function to calcualte correlation stat for each probe independently
probeWiseCor <- function(x){
  corStat <- cor.test(as.numeric(x[names(x) %in% pheno$Basename]), 
                      as.numeric(x[names(x) %in% sampleSheet$Basename]))$estimate
  return(corStat)
}

# and by cell





# apply to each row (probe) of the matrix
start.time <- Sys.time()
t <- apply(plotBetas, 1, probeWiseCor)
end.time <- Sys.time()
round(end.time - start.time,2) #1.94 mins


## plot distribution
t <- as.data.frame(t)

ggplot(t, aes(x=t))+
  geom_density()



#----------------------------------------------------------------------#
# Plot cor stat accross all samples
#----------------------------------------------------------------------#

# plot for each cell type - B-cells  CD4 T-cells  CD8 T-cells Granulocytes    Monocytes 

# create data frame for colouring distribution plot by array version
distGroup <- bind_rows(
  sampleSheet %>% 
    filter(Study == "A-Risk") %>%
    dplyr::select(Basename, Cell_Type) %>%
    mutate(tech = "MSA"),
  
  pheno %>%
    dplyr::select(Basename, Sample.Type) %>%
    mutate(tech = "epicV1") %>%
    rename(Cell_Type = Sample.Type))


pdf("plots/EpicV1vsMSAarrayByCellType.pdf")
for(cell in unique(sampleSheet$Cell_Type)[2:6]){

  # plot distribution, each cell type on seperate plot, coloured by MSA/V1

  # subset to just cell type betas
  cellDistBetas <- plotBetas[,distGroup$Basename[distGroup$Cell_Type == cell]]
  ArrayType <- distGroup$tech[distGroup$Cell_Type == cell]

  # denisty plot
  minfi::densityPlot(as.matrix(cellDistBetas), sampGroups = ArrayType, main=cell)
  
  
  # probe wise cor plot
  
  cellProbeCorStat <- apply(plotBetas, 1, probeWiseCor)


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
  
  
  # plot distribution of probe level statistics
  
  # by cell type
  
  
  # all samples
  
  
}

dev.off()

# violin plot

# need a dataframe with cols for methylation value, cell type, array

t <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(t)= c("Methylation", "array", "cell")

for(cell in unique(sampleSheet$Cell_Type)[2:6]){
  
  # get mean methylation for the cell type from A risk data
  cellBetas <- as.data.frame(rowMeans(betas[,pheno$Basename[pheno$Sample.Type == cell]]))
  colnames(cellBetas) <- "Methylation"
  #cellBetas$EPICv1_Locus_Match <- row.names(cellBetas)
  cellBetas$array <- "epicV1"
  cellBetas$cell <- cell
  
# get mean methylation for the cell type from MSA data
  meanMSA <- as.data.frame(rowMeans(MSAbetas[,sampleSheet$Basename[sampleSheet$Cell_Type == cell]]))
  colnames(meanMSA) <- "Methylation"
  #meanMSA$IlmnID <- row.names(meanMSA)
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




  
