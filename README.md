# Methylation Screening Array
QC of cell sorted methylation array data run on the Ilumina Infinium Methylation Screening BeadChip \
([https://emea.illumina.com/products/by-type/microarray-kits/infinium-mouse-methylation.html](https://emea.illumina.com/products/by-type/microarray-kits/infinium-methylation-screening-array.html))
\
\
This repository contains analysis pipelines for the pre-processing and analysis of the data generated as part of this project. They are developed for use with a HPC system with the SLURM job scheduler.


# Guidance for using the pipeline scripts

PREQUISITES: 
* A project folder containing the following sub directories:
  * 0_metadata
  * 1_raw
  * 2_normalised
  * 3_analysis
  * logFiles
  
* A config file in the project folder containing file paths and thresholds etc specific to their project (see exampleData/config.r for an example)

* A file named sampleSheet.csv in the 0_metadata folder which contains the following columns, and any other associated phenotype data:
    "Basename","Chip_ID","Chip_Position","Plate","Batch", "Individual_ID","Sample_ID","Cell_Type","N_Nuclei","Group","Age","Sex". (see exampleData/sampleSheet.csv)
    
* raw idat files in the 1_raw folder



# QUALITY CONTROL STAGE

OUTPUT: 
* Quality control objects, metrics and html reports are located in 2_normalised/QCmetrics


The QCjobSubmission.sh script automates the quality control pipeline and can be run on the command line with:
sbatch QCjobSubmission.sh <filepath/to/projectFolder>

This script will execute:
* calcQCmetrics.r
* Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args $1
* cellTypeChecks.r (tba)
* Rscript -e "rmarkdown::render('cellTypeQC.rmd', output_file='cellTypeQC.html')" (tba)
* normalisation.r (tba)



