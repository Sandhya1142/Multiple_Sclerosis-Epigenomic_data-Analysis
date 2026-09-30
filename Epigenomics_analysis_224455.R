#Epigenomic analysis of dataset GSE224455

#Loading the libraries
library(GEOquery)
library(BSgenome)
library(ChAMP)
library(FDb.InfiniumMethylation.hg19)
library(org.Hs.eg.db)
library(geneLenDataBase)
library(here)
library(limma)
library(ChAMPdata)


#################################################
#STEP 1 - Data Extraction and pre-processing
#################################################

options(timeout = 3600)
gse <- getGEO("GSE224455", getGPL=FALSE, GSEMatrix = TRUE)
gse <- gse[[1]]
pheno <- pData(gse)
pheno


#Downloading supplementary files
getGEOSuppFiles("GSE224455", baseDir=here())
list.files(here("GSE224455"))

csv <- data.frame(Sample_Name = pheno$geo_accession,
                  Sentrix_ID = pheno$`sentrixid:ch1`, 
                  Sentrix_Position =pheno$`sentrixposition:ch1`,
                  Diagnosis=pheno$title,
                  Condition = pheno$source_name_ch1,
                  Age=pheno$`age:ch1`,
                  Sex= pheno$`seks:ch1`,
                  PMI=pheno$characteristics_ch1.2)
csv$Condition <- ifelse(
  csv$Condition == "Brain lesion",
  "Brain_lesion",
  "Brain_NAWM"
)
csv$Diagnosis <- sub(" ","_",csv$Diagnosis)
csv$PMI <- sub("pmi \\(minutes\\):", "", csv$PMI)
untar(here("GSE224455", "GSE224455_RAW.tar"), exdir = here("idat"))


#Mapping Sentrix ID and position to file names
idat_original <- here("idat")
dir.create(here("idat_champ"))
csv[,c("Sample_Name", "Sentrix_ID", "Sentrix_Position", "Condition")]
args(list.files)


for (i in 1:nrow(csv)){
  gsm <- csv$Sample_Name[i]
  sentrix_id <- csv$Sentrix_ID[i]
  sentrix_position <- csv$Sentrix_Position[i]
  
  green_file <- list.files(idat_original,
                           pattern=paste0("^",gsm,"_.*_Grn\\.idat\\.gz$"),
                           full.names = TRUE)
  red_file <- list.files(idat_original,
                         pattern=paste0("^",gsm,"_.*_Red\\.idat\\.gz$"),
                         full.names = TRUE)
  
  new_green <- file.path(here("idat_champ"),
                         paste0(sentrix_id,"_",sentrix_position,"_Grn.idat.gz"))
  new_red <- file.path(here("idat_champ"),
                       paste0(sentrix_id,"_",sentrix_position,"_Red.idat.gz"))
  
  file.copy(green_file, new_green, overwrite = TRUE)
  file.copy(red_file, new_red, overwrite = TRUE)
}


#################################################
#STEP 2 - Data Loading
#################################################

write.csv(csv,file=file.path(here("idat_champ"),"data.csv"))
myLoad <- champ.load(directory=here("idat_champ"), arraytype="EPIC",detPcut=0.05)
CpG.GUI(arraytype="EPIC")



#################################################
#STEP 3 - Quality Control
#################################################

myQC <- champ.QC(beta=myLoad$beta,
                 pheno=myLoad$pd$Condition,
                 resultsDir = here("Images"))
QC.GUI(arraytype = "EPIC")


#################################################
#STEP 4 - Normalization
#################################################

myNorm <- champ.norm(beta=myLoad$beta,
                     arraytype = "EPIC",
                     resultsDir=here("Results"),
                     plotBMIQ=TRUE)
dim(myNorm)
sum(is.na(myNorm))
QC.GUI(beta=myNorm, pheno=myLoad$pd$Condition)


#################################################
#STEP 5 - Singular Value Decomposisiton (SVD)
#################################################

class(myLoad$pd)
class(myLoad$beta)
myLoad$beta <- as.data.frame(myLoad$beta)

myLoad$pd$Condition <- factor(myLoad$pd$Condition,
                                    levels = c("Brain_NAWM","Brain_lesion"))
mySVD <- champ.SVD(beta = as.data.frame(myNorm), pd = myLoad$pd)


####################################################
#STEP 6 - Batch Correction using Combat() function
####################################################
dev.off()
myCombat <- champ.runCombat(beta=myNorm,
                            pd=myLoad$pd,
                            variablename="Condition", 
                            batchname = c("Slide"))
mySVD_combat <- champ.SVD(beta = as.data.frame(myCombat),pd = myLoad$pd)


rownames(myLoad$pd) <- myLoad$pd$Sample_Name
all(rownames(myLoad$pd) == colnames(myNorm))

##################################################################
#STEP 7 - Identifying Differentially Methylated Positions (DMPs)
##################################################################

myLoad$pd$Donor <- factor(sub(",.*", "", myLoad$pd$Diagnosis))
table(myLoad$pd$Slide, myLoad$pd$Condition)

myLoad$pd$Donor <- factor(myLoad$pd$Donor)
myLoad$pd$Sex <- factor(myLoad$pd$Sex)
myLoad$pd$Age <- as.numeric(myLoad$pd$Age)
myLoad$pd$PMI <- as.numeric(myLoad$pd$PMI)

#Creating the design matrix
design <- model.matrix(~ Condition+Age+Sex+PMI, data = myLoad$pd)
design


#Creating a block to consider Paired samples
corfit <- duplicateCorrelation(myCombat,
                               design = design,
                               block = myLoad$pd$Donor)
corfit$consensus.correlation


#Fitting linear regression for all probes
fit <- lmFit(myCombat,
             design,
             block = myLoad$pd$Donor,
             correlation = corfit$consensus)
fit <- eBayes(fit)
colnames(design)
res_comb <- topTable(fit,
                     coef = "ConditionBrain_lesion",
                     number = Inf,
                     adjust.method = "BH")


#DMPs
dmp <- rownames(res_comb)[res_comb$adj.P.Val < 0.05, drop = FALSE]
all <- rownames(res_comb)


#calculating the deltaBeta value
nawm_samples <- rownames(myLoad$pd)[myLoad$pd$Condition == "Brain_NAWM"]
means_nawm <- rowMeans(myCombat[all, nawm_samples, drop=FALSE], )
means_nawm_dmp <- rowMeans(myCombat[dmp, nawm_samples, drop=FALSE], )


lesion_samples <- rownames(myLoad$pd)[myLoad$pd$Condition == "Brain_lesion"]
means_lesion <- rowMeans(myCombat[all, lesion_samples, drop=FALSE], )
means_lesion_dmp <- rowMeans(myCombat[dmp, lesion_samples, drop=FALSE], )


#delta beta for all probes
deltaBeta <- means_lesion - means_nawm

#delta beta for significant probes
deltaBeta_dmp <- means_lesion_dmp - means_nawm_dmp
deltaBeta_dmp
  
sum(deltaBeta_dmp > 0)
sum(deltaBeta_dmp < 0)
mean(deltaBeta_dmp > 0)
mean(deltaBeta_dmp < 0)

length(dmp)
dim(myCombat[dmp, , drop = FALSE])

##################################################################
#STEP 8 - Annotation of DMPs 
##################################################################

data(probe.features)
anno <- probe.features
res_comb$deltaBeta <- deltaBeta[rownames(res_comb)]
colnames(anno)
cols <- c("CHR", "MAPINFO", "Strand", "gene", "feature", "cgi", "feat.cgi", "UCSC_CpG_Islands_Name")

idx <- match(dmp, rownames(anno))
sum(is.na(idx))

dmp_anno <- cbind(dmp,
                  anno[idx, cols, drop = FALSE])
write.csv(dmp_anno,here("Annotation.csv"))

head(dmp_anno)
table(dmp_anno$feature) 
table(dmp_anno$cgi)