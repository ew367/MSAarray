##---------------------------------------------------------------------#
##
## Title: MSA array
##
## Purpose of script: 1. distribution accross sites
##                    2. DNAmAge
##                    3. PhenoAge
##                    4. DunedinPace
##                    
##
## Dataset: Extend (epic v1)
##          or Sarah's EpicV2 blood data?
##
## Author: Emma Walker
##
## Date Created: 11-03-24
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Set up
#----------------------------------------------------------------------#

library(CETYGO)
library(stringr)

dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/Extend/"
  
refDir <- "/lustre/projects/Research_Project-MRC190311/references"
v1 <- file.path(refDir, "EPICArray/MethylationEPIC_v-1-0_B4.csv")
v2 <- file.path(refDir, "EPICArray/epicV2manifest_pidsley.csv")
msa <- file.path(refDir, "MSA/MSA-48v1-0_A1.csv")


# load in Extend data
load(file.path(dataDir, "EXTEND_batches_1_2_normalised_together.rdat"))
pheno <- read.csv(file.path(dataDir,"EXTEND_pheno_batches_1_2.csv"))

# load manifests
v1Man <- read.csv(v1, stringsAsFactors = F, skip = 7)
msaMan <- read.csv(msa, stringsAsFactors = F, skip = 7)

setwd("/lustre/projects/Research_Project-MRC190311/DNAm/MSA/")


#----------------------------------------------------------------------#
# Summary and match to epic v1
#----------------------------------------------------------------------#

## useful cols = EPICv1_Locus_Match (cgs) and EPICv1_ProbeSeq_Match (T/F)
# sites in MSA that are in v1 
msaV1 <- msaMan[!(msaMan$EPICv1_Locus_Match == ""),]

# look at weird probes
w <- msaV1[which(msaV1$EPICv1_ProbeSeq_Match != TRUE & msaV1$EPICv1_ProbeSeq_Match != FALSE),]
#write.csv(w, "MSAmultipleProbeSeqMatchInEpicV1.csv")

# get betas

MSAv1Betas <- betas[row.names(betas) %in% unlist(str_split(msaV1$EPICv1_Locus_Match, ";", simplify = F)),]

densityPlot(betas)
densityPlot(MSAv1Betas)


#----------------------------------------------------------------------#
# DNAmAge
#----------------------------------------------------------------------#

data(coef)
DNAmAge<-agep(betas, coef=coef) # all betas

# count how many clock CPGs are available betas
sum(names(coef) %in% row.names(betas))
# [1] 335/354

MSAage <- agep(MSAv1Betas, coef=coef) # MSA only betas
sum(names(coef) %in% row.names(MSAv1Betas))
# 332/354


# plot
plotdf <- as.data.frame(cbind(DNAmAge, MSAage, pheno$Age))
colnames(plotdf) <- c("allBetasAgePred", "MSABetasAgePred", "ReportedAge")

cor.test(plotdf$allBetasAgePred, plotdf$MSABetasAgePred)

ggplot(plotdf,aes(x=allBetasAgePred, y=MSABetasAgePred))+
  geom_point()+
  ggtitle("DNAmAge Extend epicV1, R=0.9999")


#----------------------------------------------------------------------#
# phenoAge
#----------------------------------------------------------------------#

imputation.knn <- function(x){
  # assumes you are only going to impute a few hundred sites. 
  nsamples = ncol(x)
  datMethUsed = t(x)
  dimnames1 = dimnames(datMethUsed)
  noMissingPerSample = rowSums(is.na(datMethUsed))
  
  if(nsamples > 1 & max(noMissingPerSample) > 0){
    datMethUsed = as.data.frame(t(impute.knn(datMethUsed)$data))
    colnames(datMethUsed) = dimnames1[[1]]
    return(datMethUsed)
  } else {
    return(x)
  }
}

# Calculate PhenoAge
pheno.age <- function(beta, PhenoAgeCoeff){
  intercept <- PhenoAgeCoeff$Weight[1]
  coeffs <- PhenoAgeCoeff$Weight[-1]
  probeOverlap <- match(PhenoAgeCoeff$CpG[-1], rownames(beta))
  coeffs <- coeffs[!is.na(probeOverlap)]
  probeOverlap <- probeOverlap[!is.na(probeOverlap)]
  beta <- beta[probeOverlap,]
  imputed_beta <- imputation.knn(beta)
  imputed_dat <- as.matrix(imputed_beta)
  predAge <- intercept+coeffs %*% imputed_dat
  return(t(predAge))
}

# read in phenoAge coeffs
PhenoAgeCoeff <- read.csv("/lustre/home/ew367/SiyiAgeClock/demoAgeClock/PhenoAgeCoeff.csv", stringsAsFactors = F, fill = T)

# for all betas
V1phenoAge <- pheno.age(betas, PhenoAgeCoeff)
sum(PhenoAgeCoeff$CpG %in% row.names(betas))
#512/513


#MSAbetas
MSAphenoAge <- pheno.age(MSAv1Betas, PhenoAgeCoeff)
sum(PhenoAgeCoeff$CpG %in% row.names(MSAv1Betas))
#508



plotdf <- as.data.frame(cbind(V1phenoAge, MSAphenoAge, pheno$Age))
colnames(plotdf) <- c("allBetasAgePred", "MSABetasAgePred", "ReportedAge")

cor.test(plotdf$allBetasAgePred, plotdf$MSABetasAgePred)

ggplot(plotdf,aes(x=allBetasAgePred, y=MSABetasAgePred))+
  geom_point()+
  ggtitle("phenoAge Extend epicV1, R=0.9998")



#----------------------------------------------------------------------#
# DunedinPace
#----------------------------------------------------------------------#

# nb this needs to be done on R v4.

library(DunedinPACE)

# V1
v1PredAgeDun <- try(PACEProjector(betas, proportionOfProbesRequired=0.7))
v1paceresult <- data.frame("IID" = names(v1PredAgeDun$DunedinPACE), "v1PredAge" = v1PredAgeDun$DunedinPACE)

# MSA
msaPredAgeDun <- try(PACEProjector(MSAv1Betas, proportionOfProbesRequired=0.4)) # won't work higher than this
MSApaceresult <- data.frame("IID" = names(msaPredAgeDun$DunedinPACE), "msaPredAge" = msaPredAgeDun$DunedinPACE)


#mean(dun$PredAge)

#get CpGs for clock
dunPaceProbes <- getRequiredProbes(backgroundList = F)
sum(unlist(dunPaceProbes) %in% row.names(betas))
#159/173
sum(unlist(dunPaceProbes) %in% row.names(MSAv1Betas))
#158/173


# including background probes
# betas - 19797/20000 (99%)
# MSAbetas - 9569/20000 (48%)

load("POAv1Andmsa.rdat")
colnames(POAres)[2:3] <- c("allBetasPOA", "MSABetasPOA")

# plot

cor.test(POAres$allBetasPOA, POAres$MSABetasPOA)


ggplot(POAres,aes(x=allBetasPOA, y=MSABetasPOA))+
  geom_point()+
  xlim(0.5, 1.5)+
  ylim(0.5, 1.5)+
  geom_abline(linetype="dashed", colour="red")+
  ggtitle("Dunedin Pace of Aging, Extend epicV1, R=0.995")






#----------------------------------------------------------------------#
# CETYGO
#----------------------------------------------------------------------#

# all extend betas
rowIndex<-rownames(betas)[rownames(betas) %in% rownames(modelBloodCoef)]
predProp<-projectCellTypeWithError(betas, modelBloodCoef[rowIndex,])


# MSA probe betas only
rowIndexMSA<-rownames(MSAv1Betas)[rownames(MSAv1Betas) %in% rownames(modelBloodCoef)]
predPropMSA<-projectCellTypeWithError(MSAv1Betas, modelBloodCoef[rowIndex,])

# cor and plot
cor.test(predProp[,7], predPropMSA[,7])

cetygoDF <- as.data.frame(cbind(predProp[,7], predPropMSA[,7]))
colnames(cetygoDF) <- c("allBetas", "MSABetas")

ggplot(cetygoDF,aes(x=allBetas, y=MSABetas))+
  geom_point()+
  geom_abline(linetype="dashed", colour="red")+
  ggtitle("CETYGO score Extend epicV1, R=0.975")


#----------------------------------------------------------------------#
# Hypo/hemi/hypermethlyation plots
#----------------------------------------------------------------------#

## MSA

MSAmeanBetas <- as.data.frame(rowMeans(MSAv1Betas))
colnames(MSAmeanBetas) <- "MSAmean"
sum(MSAmeanBetas$MSAmean < 0.2)
# [1] 37288
sum(MSAmeanBetas$MSAmean > 0.8)
# [1] 31526
sum(MSAmeanBetas$MSAmean >= 0.2 & MSAmeanBetas$MSAmean <= 0.8)
#[1] 77043

MSAmeanBetas$ID <- row.names(MSAmeanBetas)


## all extend

v1meanBetas <- as.data.frame(rowMeans(betas))
colnames(v1meanBetas) <- "v1mean"
sum(v1meanBetas$v1mean < 0.2)
# [1] 215647
sum(v1meanBetas$v1mean > 0.8)
# [1] 247908
sum(v1meanBetas$v1mean >= 0.2 & v1meanBetas$v1mean <= 0.8)
#[1] 345943

# define hypo/hemi/hyper state of probes
v1meanBetas$state <- rep(NA)
v1meanBetas$state[v1meanBetas$v1mean < 0.2] <- "hypo"
v1meanBetas$state[v1meanBetas$v1mean > 0.8] <- "hyper"
v1meanBetas$state[v1meanBetas$v1mean >= 0.2 & v1meanBetas$v1mean <= 0.8] <- "hemi"

# create col for array type
v1meanBetas$ID <- row.names(v1meanBetas)
v1meanBetas$array[v1meanBetas$ID %in% MSAmeanBetas$ID] <- "MSA"
v1meanBetas$array[is.na(v1meanBetas$array)] <- "V1"


##plot

ggplot(v1meanBetas, aes(x=state, y=v1mean, colour = array))+
  geom_boxplot()


ggplot(v1meanBetas, aes(x=state, fill = array))+
geom_bar(position = "fill")


## denisty plot

densityPlot(betas, sampGroups = v1meanBetas$array)


#----------------------------------------------------------------------#
# RS probes
#----------------------------------------------------------------------#

# same 59 probes, but with different ilmIDs

rsProbes <- msaV1[msaV1$Probe_Type == "rs",]

length(unique(rsProbes$EPICv1_Locus_Match))
#[1] 59
length(unique(rsProbes$Name))
#[1] 59
length(unique(rsProbes$IlmnID))
#[1] 374

#rsbetas<-betas[grep("rs", rownames(betas)),]
#snpCor<-cor(rsbetas, use = "pairwise.complete.obs")


#----------------------------------------------------------------------#
# Variance
#----------------------------------------------------------------------#

varV1 <- apply(betas, 1, var)
varV1 <- as.data.frame((varV1))
varV1 <- cbind(varV1, v1meanBetas$state, v1meanBetas$array)
colnames(varV1) <- c("var", "state", "array")


ggplot(varV1, aes(x=state, y=log10(var), colour = array))+
  geom_boxplot()


ggplot(varV1, aes(x=array, y=log10(var), colour=array))+
  geom_boxplot()



