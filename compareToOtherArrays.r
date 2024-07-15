##---------------------------------------------------------------------#
##
## Title: correlate MSA test data to data on same samples but different
##        arrays
##
##                    
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

#----------------------------------------------------------------------#
# LOAD PACKAGES and DATA
#----------------------------------------------------------------------#


library(ENmix)
library(ggplot2)

source("config.r")


# load manifest
man <- fread(manifest, skip=7, fill=TRUE, data.table=F)

# load MSA data
load(file = file.path(QCDir, "mraw.rdat"))
MSAbetas <- getB(mraw)
sampleSheet <- read.csv(pheno, stringsAsFactors = F)


# load A risk data
load("/lustre/projects/Research_Project-MRC190311/DNAm/Arisk/2_processed/CellSortedAsthmaERisk_Normalised.rdat")
#betas and pheno



#----------------------------------------------------------------------#
# A-Risk (epic v1)
#----------------------------------------------------------------------#

#subset to overlapping samples
AriskPheno <- pheno[pheno$PatientID %in% sampleSheet$Individual_ID,]s
AriskPheno <- AriskPheno[-which(AriskPheno$Sample.Type == "Buccal" | AriskPheno$Sample.Type == "Nasal"),] #remove Buccal and Nasal cells

betas <- as.data.frame(betas)

AriskBetas <- betas[,AriskPheno$Basename]


# subset to overlapping probes

msaV1 <- man[!(man$EPICv1_Locus_Match == ""),]

AriskBetas <- AriskBetas[msaV1$EPICv1_Locus_Match,] # 161737 

AriskBetas$ProbeID <- row.names(Arisk)

# one plot for each cell type - B-cells  CD4 T-cells  CD8 T-cells Granulocytes    Monocytes 

# take mean methylation accross cell types
# also plot 8 individual samples within each cell type



# try on one cell

cell <- unique(sampleSheet$Cell_Type)[6]



MS







  
