#!/usr/bin/env Rscript

#Make sure to set the correct working directory with your GWAS files

#Check whats in the working directory
list.files(path = ".", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#Loading required libraries.
#Loading required libraries.
library(BiocManager)
library(gdsfmt)
library(snpStats)
library(downloader)
library(SNPRelate)                     
library(plyr)
library(GenABEL)
library(doParallel)
library(ggplot2)
library(data.table)
library(qqman)
library(rtracklayer)
library(tidyverse)
library(RCurl)


#accessing seperate phenotype file
clinical1 <- read.csv("NeuroGAP-P_Release5_AllSites.csv")

#subsetting controls
clinical2 <- clinical1[clinical1$is_case == 0,]

#Recategorizing to obtain MDD cases and MDD free controls
clinical2$MDD_is_case <- ifelse(clinical2$kten_total >= 13, 1, 0)

#Retaining only those without missing data on MDD phenotype
clinical2 <- clinical2[!is.na(clinical2$MDD_is_case),]
phenodata <- clinical2

#subsetting phenodata file and renaming variables
phenodata <- subset(phenodata, select=c("subj_id","MDD_is_case","msex", "ethnicity_1","age_at_iview","birth_country","living_arrange","education","marital_status","lec_composite_witnessed_by_me", "lec_composite_happenned_to_me","assist_tobacco_amt","assist_alcohol_amt","assist_cannabis_amt","cidi_q1","cidi_q2","cidi_q3","cidi_q17","bmi", "hiv_positive"))
categoricals <- c("MDD_is_case","msex","birth_country","living_arrange","ethnicity_1","education","marital_status","assist_tobacco_amt","assist_alcohol_amt","assist_cannabis_amt","cidi_q1","cidi_q2","cidi_q3","cidi_q17", "hiv_positive","bmi")
phenodata[,categoricals] <- lapply(phenodata[,categoricals] , factor)
#rename columns
names(phenodata)[names(phenodata) == 'subj_id'] <- "id"
names(phenodata)[names(phenodata) == 'MDD_is_case'] <- "phenotype"
names(phenodata)[names(phenodata) == 'cidi_q1'] <- "arthritis"
names(phenodata)[names(phenodata) == 'cidi_q2'] <- "chronic_back_orneck_pain"
names(phenodata)[names(phenodata) == 'cidi_q3'] <- "frequent_or_severe_headaches"
names(phenodata)[names(phenodata) == 'cidi_q17'] <- "cancer"
names(phenodata)[names(phenodata) == 'assist_tobacco_amt'] <- "tobacco_use"
names(phenodata)[names(phenodata) == 'assist_alcohol_amt'] <- "alcohol_use"
names(phenodata)[names(phenodata) == 'assist_cannabis_amt'] <- "cannabis_use"
names(phenodata)[names(phenodata) == 'hiv_positive'] <- "HIV_status"
names(phenodata)[names(phenodata) == 'lec_composite_witnessed_by_me'] <- "Negative_life_events_witnessed"
names(phenodata)[names(phenodata) == 'lec_composite_happenned_to_me'] <- "Negative_life_events_experienced"
names(phenodata)[names(phenodata) == 'age_at_iview'] <- "age"
names(phenodata)[names(phenodata) == 'msex'] <- "gender"
names(phenodata)[names(phenodata) == 'living_arrange'] <- "living_arrangement"
names(phenodata)[names(phenodata) == 'marital_status'] <- "marital_status"
names(phenodata)[names(phenodata) == 'birth_country'] <- "birth_country"
names(phenodata)[names(phenodata) == 'ethnicity_1'] <- "ethnicity"

str(phenodata)

#GWAS
#read in data
data <- read.plink("NeuroGAP_pilot_clean.bed","NeuroGAP_pilot_clean.bim",
                   "NeuroGAP_pilot_clean.fam",na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)

#check out data
class(data)
length(data)

str(data)

#bim file contains information related to SNPs, allele 1 = minor allele, allele 2 = major allele
bim <- data$map
head(bim)
table(bim$chromosome)

#bed file contains genotype data originally put in binary format
bed <- data$genotypes
bed

#Extract IDs that will be used to determine the cases and controls from the phenotype data
IDs <- row.names(bed)
print(head(IDs))

#fam file containing participant ID
fam <- data$fam
head(fam)

#subsetting phenotype file to only contain the same numbers as in genotype data
clinical <- subset(phenodata, phenodata$id %in% IDs)
rownames(clinical) <- clinical$id
print(head(clinical))
nrow(clinical)

#saving the data for steps to follow
genotype <- data$genotypes
#subsetting genotype
genotype <- genotype[rownames(clinical),]

genoBim <- data$map
colnames(genoBim)<-c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))

#create a dataframe of snp summary statistics
snpsum.col <- col.summary(genotype)
snpsum.row <- row.summary(genotype)
head(snpsum.col)
head(snpsum.row)

#SNP level filtering
#call rate; call rate is defined as the proportion of inidividuals in our study for which the corresponding SNP information is not missing 
#a call rate of 95% for a certain SNP means that 95% of the individuals have data for this SNP. 
#Thus to filter by call rate, we look at the column of Call.rate in our dataframe and choose the SNPs that has percentage of missing observations below our chosen threshold
call <- 0.95
use <- with(snpsum.col, (!is.na(Call.rate) & Call.rate >= call))
use[is.na(use)] <- FALSE              
cat(ncol(genotype)-sum(use),"SNPs will be removed due to low call rate.\n") 
#Now that we know what SNPs to keep, we proceed to subset genotype and SNP summary data 
#so that they contains only rows of SNPs passing the call rate
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]
print(genotype) 

#minor allele frequency, Minor allele frequency (also known as variant allele frequency) is defined as frequency of the less common allele at a variable site.
#Inadequate power to infer a statistically significant relationship between the SNP and the trait under study is the result of a large degree of homogeneity at a given SNP across study participants.
#This occurs when we have a very small minor allele frequency (MAF), which means the majority of the individuals have two copies of the same major alleles
minor <- 0.01
use1 <- with(snpsum.col, (!is.na(MAF) & MAF > minor) )
use1[is.na(use1)] <- FALSE               
cat(ncol(genotype)-sum(use1),"SNPs will be removed due to low MAF .\n"  ) 

#remove snps with a minor allele frequency of less that 1%
genotype <- genotype[,use1]
snpsum.col <- snpsum.col[use1,]
print(genotype)

#Sample level filtering
#The second stage of data pre-processing involves filtering samples, i.e. removing individuals due to missing data, sample contamination, 
#correlation (for population-based investigations), and racial/ethnic or gender ambiguity or discordance

#Basic sample level filtering
#heterozygosity
#We remove any sample with inbreeding coefficient |F| > 0.10, which might indicate interbreeding or bad quality sample.
#We also remove any subjects with NA values.
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)
head(snpsum.row)

hetcutoff <- 0.1   
sampleuse <- with(snpsum.row, abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE  
cat(nrow(genotype)-sum(sampleuse), "subjects will be removed due to inbreeding coefficient.\n")

#we subset the genotype and clinical data set using information obtained above to remove the subjects that do not pass the criteria from our sample.
genotype <- genotype[sampleuse,]
clinical<- clinical[rownames(genotype), ]

#Cryptic relatedness, duplicates and gender identity
#Population-based cohort studies often only concerns unrelated individuals and the generalized linear modeling approach assumes independence across individuals. In regional cohort studies (e.g. hospital-based cohort studies) of complex diseases, however,
#individuals from the same family could be unintentionally recruited. 
#To mitigate this problem, we thus employ Identity-by-descent (IBD) analysis which is a common measurement of relatedness and duplication.
ld.thresh <- 0.2    
kin.thresh <- 0.1   
snpgdsBED2GDS("NeuroGAP_pilot_clean.bed","NeuroGAP_pilot_clean.fam", "NeuroGAP_pilot_clean.bim","NeuroGAP_pilot_clean.gds", cvt.chr="char")
genofile <- snpgdsOpen("NeuroGAP_pilot_clean.gds", readonly = FALSE)
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)
#After creating and removing samples with automatically added suffixes, 
#we have to prune the SNPs based on linkage disequilibrium for IBD analysis.
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, 
                          snp.id = colnames(genotype)) 
snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd),"will be used in IBD analysis\n")

#When all the SNPs are pruned, we need to find the IBD coefficients using Method of Moments procedure, including pairwise kinship
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
ibdcoeff <- snpgdsIBDSelection(ibd)    
head(ibdcoeff)
#Now we have IBD coefficient, we can check to see if there are any candidates for relatedness using kinship threshold.
#We remove any samples with high kinship starting with the sample with the most pairings.
ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]
related.samples <- NULL
while ( nrow(ibdcoeff) > 0 ) {
  sample.counts <- arrange(plyr::count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, 'x']
  cat("Removing sample", as.character(rm.sample), 'too closely related to', 
      sample.counts[1, 'freq'],'other samples.\n')
  
  ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}
#Having removed the related samples, we proceed to filter the genotype and clinical data to include only unrelated samples.
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$id %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", 
    kin.thresh,"\n") 
print(genotype)

#Ancestry
#Principle components (PCs) from multi-dimensional scaling (MDS) 
#self-reported race/ethnicity can differ from clusters of individuals that are based on genetics information
#the presence of an individual not appearing to fall within a racial/ethnic cluster may be suggestive of a sample error.
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.ibd, num.thread=1)
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    
                    PC2 = pca$eigenvect[,2],
                    PC3 = pca$eigenvect[,3],
                    PC4 = pca$eigenvect[,4],
                    PC5 = pca$eigenvect[,5],
                    PC6 = pca$eigenvect[,6],
                    PC7 = pca$eigenvect[,7],
                    PC8 = pca$eigenvect[,8],
                    PC9 = pca$eigenvect[,9],
                    PC10 = pca$eigenvect[,10],
                    stringsAsFactors = FALSE)

#obtaining pedigree for smoother plotting
pedigree <- fam[,-c(3:7)]
pedigree$sample.id <- row.names(pedigree)
common_col_names <- intersect(names(pedigree), names(pctab))
pctab <- merge(pctab, pedigree, by=common_col_names, all.x=TRUE)
colnames(pctab)[colnames(pctab) == 'pedigree'] <- 'study_site'
print(head(pctab))

#plotting the pca plot
png(filename = "PCA_plot.png",width = 5*600, height = 5*600, res = 300, pointsize = 12 ,bg = "white")
pctab %>%
  ggplot + aes(x = PC1, y = PC2, color=study_site) + geom_point() +
  xlab('PC1') + ylab('PC2') + theme(axis.text = element_text(size = 11)) +
  theme(axis.text = element_text(face = "bold")) + theme(legend.text = element_text(face = "bold"))
dev.off()

closefn.gds(genofile)

#SNP level filtering (Part II) - Hardy-Weinberg Equilibrium
#After filtering the sample, we proceed to filter the SNPs with Hardy Weinberg Equilibrium (HWE)
#Violation of HWE can be an indication of the presence of population substructure or occurence of a genotyping error. While it is not always true, it is common to remove any SNPs that violated HWE to avoid genotyping error
#measured with a chi-square goodness-of-fit test between the observed and expected genotypes.

hardy <- 10^-6
x <- clinical[!is.na(clinical$phenotype),]
MDD_controls <- x[x$phenotype==0, 'id']
#Filtering to only retain controls
MDD_controls_genotype <- genotype[row.names(genotype) %in% MDD_controls,]
snpsum.colCont <- col.summary(MDD_controls_genotype)
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)
HWEuse[is.na(HWEuse)] <- FALSE          
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n") 
genotype <- genotype[,HWEuse]
print(genotype)
nrow(genotype)

# Data Generation
# Re-computing PCA 
# preparing for computing principal components(PCs).
# loading data and reading in the gds file for SNPRelate fuctions

genofile <- snpgdsOpen("NeuroGAP_pilot_clean.gds", readonly = TRUE)

# generating data comprised of principal components(PCs), which are intended to capture information on latent population substructure that is typically not available in self-reported race and ethnicity variables.
# Substructure refers to the presence of genetic diversity (e.g. different allele frequencies) within an apparently homogeneous population that is due to population genetic history (e.g. migration, selection, and/or ethnic integration).
# Pruning the SNPs based on their linkage disequilibrium(LD)
ld.thresh <- 0.2  # setting the LD threshold to 0.2. 
#set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh, sample.id = geno.sample.ids, 
                          snp.id = colnames(genotype))

snpset.pca <- unlist(snpSUB, use.names = FALSE)
cat(length(snpset.pca), "SNPs will be used the IBD analysis\n") 

# calculating principal components to be included as covariates in the GWA models(the first 10 principal components in GWA models)
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.pca, num.thread=1) ##calculating PCs using the snpgdsPCA function from SNPRelate and the pruned data. 

pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")
colnames(pcs)[1]<-paste("subj_id")
print(head(pcs))
closefn.gds(genofile)

#PreGWAA steps
#phenodata argument
#First we assign phenoSub to be phenodata, which is created from merging the clinical and pcs data. 
#remove anyone that doesn't have MDD phenotype data available, and subset the genotype data for the remaining individuals. 
genodata <- genotype
clinical3 <- clinical
names(pcs)[names(pcs) == "subj_id"] <- "id"
phenoSub <- merge(clinical3,pcs)  
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
phenodata <- phenoSub
genodata <- genodata[as.character(phenodata$id),]
cat(nrow(genodata), "samples included in analysis.\n")
print(head(phenodata))
nrow(phenodata)
str(phenodata)

#subsetting phenodata file to retain only major variables otherwise GWAA will take too long or R studio will crasH. Also converting categoies to factors
#phenodata <- phenoSub
phenodata <- phenodata[,c("id","phenotype","gender", "age","ethnicity","birth_country","education","Negative_life_events_witnessed", "Negative_life_events_experienced","alcohol_use","arthritis","chronic_back_orneck_pain","frequent_or_severe_headaches","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
str(phenodata)

#ensuring rows are the same length
genodata <- genodata[as.character(phenodata$id),]
nrow(genodata)
nrow(phenodata)

#GWAA
#data preparation
#The GWAA function requires two arguments. The genodata argument should specify the entire genotype data object in SnpMatrix format.
#The phenodata argument should be a data frame with a column of sample IDs, corresponding to the row names of genodata, followed by columns for the outcome variable and covariates to be included in the model
#Print the number of SNPs to be checked
cat(paste(ncol(genodata), " SNPs included in analysis.\n"))
#create text file for GWAA output to be written to
columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
write.table(t(columns), "GWAA.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#GWAA function
#The function itself has four major elements
#(1) determines how many cores to use/which method to employ for parallel processing (depending on specified cores and operating system)
#(2) converts the genotype data into the necessary format for an additive model
#(3) fits a linear model for each SNP and covariates (one group at a time) 
#(4) writes the SNP coefficient details (the estimate, standard error, test statistic and p value) from each model out to the GWAA.txt file. Read through the function to get a general idea of the format in regards to these four steps
GWAA <- function(genodata, phenodata, filename = NULL, append = FALSE, workers = getOption("mc.cores", 
                                                                                           2L), flip = TRUE, select.snps = NULL, hosts = NULL) {
  
  if (is.null(hosts)) {
    cl <- makeCluster(workers)
  } else {
    cl <- makeCluster(hosts, "PSOCK")
  }
  show(cl)
  registerDoParallel(cl)
  
  # Function that will change which allele is counted (major or minor)
  flip.matrix <- function(x) {
    zero2 <- which(x == 0)
    two0 <- which(x == 2)
    x[zero2] <- 2
    x[two0] <- 0
    return(x)
  }
  
  foreach(part = 1:nSplits) %do% {
    genoNum <- as(genodata[, snp.start[part]:snp.stop[part]], "numeric")
    # flip.matrix function employed
    if (isTRUE(flip)) 
      genoNum <- flip.matrix(genoNum)
    rsVec <- colnames(genoNum)
    res <- foreach(snp.name = rsVec, .combine = "rbind") %dopar% {
      a <- summary(glm(phenotype ~ . - id, family = binomial, data = cbind(phenodata, 
                                                                           snp = genoNum[, snp.name])))
      a$coefficients["snp", ]
    }
    
    write.table(cbind(rsVec, res), filename, append = TRUE, quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
    cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 
                100 * part/nSplits))
  }
  
  stopCluster(cl)
  return(print("Done."))
}

#Analysis
nSNPs <- ncol(genodata)
nSplits <- 10
genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

#snps number in each subset
start <- Sys.time()
GWAA(genodata, phenodata, filename="GWAA.txt")

end <- Sys.time()
print(end-start)

# Imputation of SNPs. To fill in the gaps in the chromosomes since the genotype data was not obtained from whole genome sequencing but from genotyping arrays
#Imputation of non-typed snps will be done basing on the 1000genomes phase 3 African data
# preparing for imputing non-typed SNPs

thougeno1 <- read.pedfile("afr_phase3_chr1.ped", snps  = "afr_phase3_chr1.map")
thougeno2 <- read.pedfile("afr_phase3_chr2.ped", snps  = "afr_phase3_chr2.map")
thougeno3 <- read.pedfile("afr_phase3_chr3.ped", snps  = "afr_phase3_chr3.map")
thougeno4 <- read.pedfile("afr_phase3_chr4.ped", snps  = "afr_phase3_chr4.map")
thougeno5 <- read.pedfile("afr_phase3_chr5.ped", snps  = "afr_phase3_chr5.map")
thougeno6 <- read.pedfile("afr_phase3_chr6.ped", snps  = "afr_phase3_chr6.map")
thougeno7 <- read.pedfile("afr_phase3_chr7.ped", snps  = "afr_phase3_chr7.map")
thougeno8 <- read.pedfile("afr_phase3_chr8.ped", snps  = "afr_phase3_chr8.map")
thougeno9 <- read.pedfile("afr_phase3_chr9.ped", snps  = "afr_phase3_chr9.map")
thougeno10 <- read.pedfile("afr_phase3_chr10.ped", snps  = "afr_phase3_chr10.map")
thougeno11 <- read.pedfile("afr_phase3_chr11.ped", snps  = "afr_phase3_chr11.map")
thougeno12 <- read.pedfile("afr_phase3_chr12.ped", snps  = "afr_phase3_chr12.map")
thougeno13<- read.pedfile("afr_phase3_chr13.ped", snps  = "afr_phase3_chr13.map")
thougeno14 <- read.pedfile("afr_phase3_chr14.ped", snps  = "afr_phase3_chr14.map")
thougeno15 <- read.pedfile("afr_phase3_chr15.ped", snps  = "afr_phase3_chr15.map")
thougeno16 <- read.pedfile("afr_phase3_chr16.ped", snps  = "afr_phase3_chr16.map")
thougeno17 <- read.pedfile("afr_phase3_chr17.ped", snps  = "afr_phase3_chr17.map")
thougeno18 <- read.pedfile("afr_phase3_chr18.ped", snps  = "afr_phase3_chr18.map")
thougeno19 <- read.pedfile("afr_phase3_chr19.ped", snps  = "afr_phase3_chr19.map")
thougeno20 <- read.pedfile("afr_phase3_chr20.ped", snps  = "afr_phase3_chr20.map")
thougeno21 <- read.pedfile("afr_phase3_chr21.ped", snps  = "afr_phase3_chr21.map")
thougeno22 <- read.pedfile("afr_phase3_chr22.ped", snps  = "afr_phase3_chr22.map")

# obtaining genotype data and each SNPs position per chromosome
genoMatrix1 <- thougeno1$genotypes
support1 <- thougeno1$map
support1 <- support1[,-3]
support1 <- support1[,c(2,1,3,4,5)]
colnames(support1) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix2 <- thougeno2$genotypes
support2 <- thougeno2$map
support2 <- support2[,-3]
support2 <- support2[,c(2,1,3,4,5)]
colnames(support2) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix3 <- thougeno3$genotypes
support3 <- thougeno3$map
support3 <- support3[,-3]
support3 <- support3[,c(2,1,3,4,5)]
colnames(support3) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix4<- thougeno4$genotypes
support4 <- thougeno4$map
support4 <- support4[,-3]
support4 <- support4[,c(2,1,3,4,5)]
colnames(support4) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix5<- thougeno5$genotypes
support5 <- thougeno5$map
support5 <- support5[,-3]
support5 <- support5[,c(2,1,3,4,5)]
colnames(support5) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix6<- thougeno6$genotypes
support6 <- thougeno6$map
support6 <- support6[,-3]
support6 <- support6[,c(2,1,3,4,5)]
colnames(support6) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix7<- thougeno7$genotypes
support7 <- thougeno7$map
support7 <- support7[,-3]
support7 <- support7[,c(2,1,3,4,5)]
colnames(support7) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix8<- thougeno8$genotypes
support8 <- thougeno8$map
support8 <- support8[,-3]
support8 <- support8[,c(2,1,3,4,5)]
colnames(support8) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix9<- thougeno9$genotypes
support9 <- thougeno9$map
support9 <- support9[,-3]
support9 <- support9[,c(2,1,3,4,5)]
colnames(support9) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix10<- thougeno10$genotypes
support10 <- thougeno10$map
support10 <- support10[,-3]
support10 <- support10[,c(2,1,3,4,5)]
colnames(support10) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix11<- thougeno11$genotypes
support11 <- thougeno11$map
support11 <- support11[,-3]
support11 <- support11[,c(2,1,3,4,5)]
colnames(support11) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix12<- thougeno12$genotypes
support12 <- thougeno12$map
support12 <- support12[,-3]
support12 <- support12[,c(2,1,3,4,5)]
colnames(support12) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix13<- thougeno13$genotypes
support13 <- thougeno13$map
support13 <- support13[,-3]
support13 <- support13[,c(2,1,3,4,5)]
colnames(support13) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix14<- thougeno14$genotypes
support14 <- thougeno14$map
support14 <- support14[,-3]
support14 <- support14[,c(2,1,3,4,5)]
colnames(support14) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix15<- thougeno15$genotypes
support15 <- thougeno15$map
support15 <- support15[,-3]
support15 <- support15[,c(2,1,3,4,5)]
colnames(support15) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix16<- thougeno16$genotypes
support16 <- thougeno16$map
support16 <- support16[,-3]
support16 <- support16[,c(2,1,3,4,5)]
colnames(support16) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix17<- thougeno17$genotypes
support17 <- thougeno17$map
support17 <- support17[,-3]
support17 <- support17[,c(2,1,3,4,5)]
colnames(support17) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix18<- thougeno18$genotypes
support18 <- thougeno18$map
support18 <- support18[,-3]
support18 <- support18[,c(2,1,3,4,5)]
colnames(support18) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix19<- thougeno19$genotypes
support19 <- thougeno19$map
support19 <- support19[,-3]
support19 <- support19[,c(2,1,3,4,5)]
colnames(support19) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix20<- thougeno20$genotypes
support20 <- thougeno20$map
support20 <- support20[,-3]
support20 <- support20[,c(2,1,3,4,5)]
colnames(support20) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix21<- thougeno21$genotypes
support21 <- thougeno21$map
support21 <- support21[,-3]
support21 <- support21[,c(2,1,3,4,5)]
colnames(support21) <- c("SNP", "chr", "position", "A1", "A2")

genoMatrix22<- thougeno22$genotypes
support22 <- thougeno22$map
support22 <- support22[,-3]
support22 <- support22[,c(2,1,3,4,5)]
colnames(support22) <- c("SNP", "chr", "position", "A1", "A2")

# imputing non-typed SNPS
# imputing the non-typed SNPs  by determining the subset of SNPs that are in 1000G but not in our data.
# subsetting the data for and then creating "missing" and "present" snpMatrix objects. 
#Both are needed for imputation rules step of the algorithm. The "present" matrix holds the typed SNPs, while the “missing” matrix holds the SNPs that need to be imputed.

#chromosome by chromosome

#1
presSnps <- colnames(genotype)
presDatChr <- genoBim[genoBim$SNP %in% presSnps, ]
targetSnps <- presDatChr$SNP
is.present <- colnames(genoMatrix1) %in% targetSnps
missing <- genoMatrix1[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix1[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support1$position[is.present]
pos.miss <- support1$position[!is.present]
rules <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules <- rules[imputation.r2(rules) >= r2threshold]

rules <- rules[imputation.maf(rules) >= minor]
cat(length(rules),"imputation rules remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules),"imputation rules remain after uncertain imputations were removed\n") 

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix1)
rm(missing)
rm(present)

#2
is.present <- colnames(genoMatrix2) %in% targetSnps
missing <- genoMatrix2[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix2[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support2$position[is.present]
pos.miss <- support2$position[!is.present]
rules2 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules2” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules2 <- rules2[can.impute(rules2)]
cat("Imputation rules2 for", length(rules2), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules2” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules2 <- rules2[imputation.r2(rules2) >= r2threshold]

rules2 <- rules2[imputation.maf(rules2) >= minor]
cat(length(rules2),"imputation rules2 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules2, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules2),"imputation rules2 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix2)
rm(missing)
rm(present)

#3
is.present <- colnames(genoMatrix3) %in% targetSnps
missing <- genoMatrix3[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix3[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support3$position[is.present]
pos.miss <- support3$position[!is.present]
rules3 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules3” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules3 <- rules3[can.impute(rules3)]
cat("Imputation rules3 for", length(rules3), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules3” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules3 <- rules3[imputation.r2(rules3) >= r2threshold]

rules3 <- rules3[imputation.maf(rules3) >= minor]
cat(length(rules3),"imputation rules3 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules3, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules3),"imputation rules3 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix3)
rm(missing)
rm(present)

#4
is.present <- colnames(genoMatrix4) %in% targetSnps
missing <- genoMatrix4[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix4[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support4$position[is.present]
pos.miss <- support4$position[!is.present]
rules4 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules4” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules4 <- rules4[can.impute(rules4)]
cat("Imputation rules4 for", length(rules4), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules4” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules4 <- rules4[imputation.r2(rules4) >= r2threshold]

rules4 <- rules4[imputation.maf(rules4) >= minor]
cat(length(rules4),"imputation rules4 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules4, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules4),"imputation rules4 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix4)
rm(missing)
rm(present)

#5
is.present <- colnames(genoMatrix5) %in% targetSnps
missing <- genoMatrix5[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix5[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support5$position[is.present]
pos.miss <- support5$position[!is.present]
rules5 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules5” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules5 <- rules5[can.impute(rules5)]
cat("Imputation rules5 for", length(rules5), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules5” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules5 <- rules5[imputation.r2(rules5) >= r2threshold]

rules5 <- rules5[imputation.maf(rules5) >= minor]
cat(length(rules5),"imputation rules5 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules5, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules5),"imputation rules5 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix5)
rm(missing)
rm(present)

#6
is.present <- colnames(genoMatrix6) %in% targetSnps
missing <- genoMatrix6[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix6[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support6$position[is.present]
pos.miss <- support6$position[!is.present]
rules6 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules6” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules6 <- rules6[can.impute(rules6)]
cat("Imputation rules6 for", length(rules6), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules6” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules6 <- rules6[imputation.r2(rules6) >= r2threshold]

rules6 <- rules6[imputation.maf(rules6) >= minor]
cat(length(rules6),"imputation rules6 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules6, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules6),"imputation rules6 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix6)
rm(missing)
rm(present)

#7
is.present <- colnames(genoMatrix7) %in% targetSnps
missing <- genoMatrix7[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix7[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support7$position[is.present]
pos.miss <- support7$position[!is.present]
rules7 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules7” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules7 <- rules7[can.impute(rules7)]
cat("Imputation rules7 for", length(rules7), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules7” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules7 <- rules7[imputation.r2(rules7) >= r2threshold]

rules7 <- rules7[imputation.maf(rules7) >= minor]
cat(length(rules7),"imputation rules7 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules7, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules7),"imputation rules7 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix7)
rm(missing)
rm(present)

#8
is.present <- colnames(genoMatrix8) %in% targetSnps
missing <- genoMatrix8[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix8[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support8$position[is.present]
pos.miss <- support8$position[!is.present]
rules8 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules8” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules8 <- rules8[can.impute(rules8)]
cat("Imputation rules8 for", length(rules8), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules8” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules8 <- rules8[imputation.r2(rules8) >= r2threshold]

rules8 <- rules8[imputation.maf(rules8) >= minor]
cat(length(rules8),"imputation rules8 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules8, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules8),"imputation rules8 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix8)
rm(missing)
rm(present)

#9
is.present <- colnames(genoMatrix9) %in% targetSnps
missing <- genoMatrix9[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix9[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support9$position[is.present]
pos.miss <- support9$position[!is.present]
rules9 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules8” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules9 <- rules9[can.impute(rules9)]
cat("Imputation rules9 for", length(rules9), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules9” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules9 <- rules9[imputation.r2(rules9) >= r2threshold]

rules9 <- rules9[imputation.maf(rules9) >= minor]
cat(length(rules9),"imputation rules9 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules9, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules9),"imputation rules9 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix9)
rm(missing)
rm(present)

#10
is.present <- colnames(genoMatrix10) %in% targetSnps
missing <- genoMatrix10[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix10[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support10$position[is.present]
pos.miss <- support10$position[!is.present]
rules10 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules10” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules10 <- rules10[can.impute(rules10)]
cat("Imputation rules10 for", length(rules10), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules10” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules10 <- rules10[imputation.r2(rules10) >= r2threshold]

rules10 <- rules10[imputation.maf(rules10) >= minor]
cat(length(rules10),"imputation rules10 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules10, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules10),"imputation rules10 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix10)
rm(missing)
rm(present)

#11
is.present <- colnames(genoMatrix11) %in% targetSnps
missing <- genoMatrix11[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix11[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support11$position[is.present]
pos.miss <- support11$position[!is.present]
rules11 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules11” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1100G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules11 <- rules11[can.impute(rules11)]
cat("Imputation rules11 for", length(rules11), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules11” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules11 <- rules11[imputation.r2(rules11) >= r2threshold]

rules11 <- rules11[imputation.maf(rules11) >= minor]
cat(length(rules11),"imputation rules11 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules11, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules11),"imputation rules11 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix11)
rm(missing)
rm(present)

#12
is.present <- colnames(genoMatrix12) %in% targetSnps
missing <- genoMatrix12[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix12[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support12$position[is.present]
pos.miss <- support12$position[!is.present]
rules12 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules12” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1200G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules12 <- rules12[can.impute(rules12)]
cat("Imputation rules12 for", length(rules12), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules12” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules12 <- rules12[imputation.r2(rules12) >= r2threshold]

rules12 <- rules12[imputation.maf(rules12) >= minor]
cat(length(rules12),"imputation rules12 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules12, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules12),"imputation rules12 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix12)
rm(missing)
rm(present)

#13
is.present <- colnames(genoMatrix13) %in% targetSnps
missing <- genoMatrix13[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix13[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support13$position[is.present]
pos.miss <- support13$position[!is.present]
rules13 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules13” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1300G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules13 <- rules13[can.impute(rules13)]
cat("Imputation rules13 for", length(rules13), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules13” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules13 <- rules13[imputation.r2(rules13) >= r2threshold]

rules13 <- rules13[imputation.maf(rules13) >= minor]
cat(length(rules13),"imputation rules13 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules13, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules13),"imputation rules13 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix13)
rm(missing)
rm(present)

#14
is.present <- colnames(genoMatrix14) %in% targetSnps
missing <- genoMatrix14[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix14[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support14$position[is.present]
pos.miss <- support14$position[!is.present]
rules14 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules14” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1400G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules14 <- rules14[can.impute(rules14)]
cat("Imputation rules14 for", length(rules14), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules14” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules14 <- rules14[imputation.r2(rules14) >= r2threshold]

rules14 <- rules14[imputation.maf(rules14) >= minor]
cat(length(rules14),"imputation rules14 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules14, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules14),"imputation rules14 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix14)
rm(missing)
rm(present)

#15
is.present <- colnames(genoMatrix15) %in% targetSnps
missing <- genoMatrix15[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix15[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support15$position[is.present]
pos.miss <- support15$position[!is.present]
rules15 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules15” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1500G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules15 <- rules15[can.impute(rules15)]
cat("Imputation rules15 for", length(rules15), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules15” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules15 <- rules15[imputation.r2(rules15) >= r2threshold]

rules15 <- rules15[imputation.maf(rules15) >= minor]
cat(length(rules15),"imputation rules15 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules15, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules15),"imputation rules15 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix15)
rm(missing)
rm(present)

#16
is.present <- colnames(genoMatrix16) %in% targetSnps
missing <- genoMatrix16[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix16[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support16$position[is.present]
pos.miss <- support16$position[!is.present]
rules16 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules16” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1600G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules16 <- rules16[can.impute(rules16)]
cat("Imputation rules16 for", length(rules16), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules16” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules16 <- rules16[imputation.r2(rules16) >= r2threshold]

rules16 <- rules16[imputation.maf(rules16) >= minor]
cat(length(rules16),"imputation rules16 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules16, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules16),"imputation rules16 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix16)
rm(missing)
rm(present)

#17
is.present <- colnames(genoMatrix17) %in% targetSnps
missing <- genoMatrix17[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix17[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support17$position[is.present]
pos.miss <- support17$position[!is.present]
rules17 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules17” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1700G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules17 <- rules17[can.impute(rules17)]
cat("Imputation rules17 for", length(rules17), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules17” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules17 <- rules17[imputation.r2(rules17) >= r2threshold]

rules17 <- rules17[imputation.maf(rules17) >= minor]
cat(length(rules17),"imputation rules17 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules17, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules17),"imputation rules17 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix17)
rm(missing)
rm(present)

#18
is.present <- colnames(genoMatrix18) %in% targetSnps
missing <- genoMatrix18[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix18[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support18$position[is.present]
pos.miss <- support18$position[!is.present]
rules18 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules18” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1800G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules18 <- rules18[can.impute(rules18)]
cat("Imputation rules18 for", length(rules18), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules18” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules18 <- rules18[imputation.r2(rules18) >= r2threshold]

rules18 <- rules18[imputation.maf(rules18) >= minor]
cat(length(rules18),"imputation rules18 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules18, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules18),"imputation rules18 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix18)
rm(missing)
rm(present)

#19
is.present <- colnames(genoMatrix19) %in% targetSnps
missing <- genoMatrix19[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix19[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support19$position[is.present]
pos.miss <- support19$position[!is.present]
rules19 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules19” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 1900G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules19 <- rules19[can.impute(rules19)]
cat("Imputation rules19 for", length(rules19), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules19” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules19 <- rules19[imputation.r2(rules19) >= r2threshold]

rules19 <- rules19[imputation.maf(rules19) >= minor]
cat(length(rules19),"imputation rules19 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules19, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules19),"imputation rules19 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix19)
rm(missing)
rm(present)

#20
is.present <- colnames(genoMatrix20) %in% targetSnps
missing <- genoMatrix20[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix20[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support20$position[is.present]
pos.miss <- support20$position[!is.present]
rules20 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules20” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 2000G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules20 <- rules20[can.impute(rules20)]
cat("Imputation rules20 for", length(rules20), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules20” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules20 <- rules20[imputation.r2(rules20) >= r2threshold]

rules20 <- rules20[imputation.maf(rules20) >= minor]
cat(length(rules20),"imputation rules20 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules20, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules20),"imputation rules20 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix20)
rm(missing)
rm(present)

#21
is.present <- colnames(genoMatrix21) %in% targetSnps
missing <- genoMatrix21[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix21[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support21$position[is.present]
pos.miss <- support21$position[!is.present]
rules21 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules21” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 2100G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules21 <- rules21[can.impute(rules21)]
cat("Imputation rules21 for", length(rules21), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules21” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules21 <- rules21[imputation.r2(rules21) >= r2threshold]

rules21 <- rules21[imputation.maf(rules21) >= minor]
cat(length(rules21),"imputation rules21 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules21, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules21),"imputation rules21 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix21)
rm(missing)
rm(present)

#22
is.present <- colnames(genoMatrix22) %in% targetSnps
missing <- genoMatrix22[, !is.present] #the missing matrix holds the SNPs that need to be imputed.
print(missing)
cat(length(missing),"SNPS are missing\n")

# the present matrix that holds our typed SNPs
present <- genoMatrix22[, is.present]
print(present)
cat(length(present),"SNPs were typed\n")

# acquiring  typed SNPs positions to be used for imputation
pos.pres <- support22$position[is.present]
pos.miss <- support22$position[!is.present]
rules22 <- snp.imputation(present, missing, pos.pres, pos.miss)
# deriving these “rules22” for the additional SNPs that were not typed in our study using the snp.imputation function

# using the genotypes from the 2200G data to fill in the gaps in our data. Each rule represents a predictive model for genotypes of un-typed SNPs associated with near-by typed SNPs.
rules22 <- rules22[can.impute(rules22)]
cat("Imputation rules22 for", length(rules22), "SNPs were estimated\n")

# removing un-typed SNPs that we fail to derive imputation “rules22” for and filtering out the SNPs that have low estimated minor allele frequency (MAF), and low imputation accuracy. 
#  Quality control for imputation certainty and MAF
r2threshold <- 0.7 #using 0.7 for the R^2 threshold
minor <- 0.01 #using 0.01 for the MAF
rules22 <- rules22[imputation.r2(rules22) >= r2threshold]

rules22 <- rules22[imputation.maf(rules22) >= minor]
cat(length(rules22),"imputation rules22 remain after MAF filtering\n")

# Obtaining posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules22, target, as.numeric=FALSE)
print(imputed)  
cat(length(rules22),"imputation rules22 remain after uncertain imputations were removed\n")

# removing the data we created(genoMatrix,missing, present) to open up memory.
rm(genoMatrix22)
rm(missing)
rm(present)
#Model fitting of imputed SNPs
#Here we use the snp.rhs.tests function from snpStats to fit a generalized linear model similarly
#snp.rhs.tests can handle imputed SNPs with the optional parameter 'rules'. This allows the function to take into account the uncertainty that comes with imputed SNPs

#1
rownames(phenodata) <- as.character(phenodata$id)

imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                      + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules)

results1 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results1 <- na.omit(results1)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut1 <- merge(results1, support1[, c("SNP", "chr", "position")])

#2
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules2)

results2 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results2 <- na.omit(results2)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut2 <- merge(results2, support2[, c("SNP", "chr", "position")])

#3
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules3)

results3 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results3 <- na.omit(results3)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut3 <- merge(results3, support3[, c("SNP", "chr", "position")])

#4
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules4)

results4 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results4 <- na.omit(results4)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut4 <- merge(results4, support4[, c("SNP", "chr", "position")])

#5
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules5)

results5 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results5 <- na.omit(results5)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut5 <- merge(results5, support5[, c("SNP", "chr", "position")])

#6
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules6)

results6 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results6 <- na.omit(results6)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut6 <- merge(results6, support6[, c("SNP", "chr", "position")])

#7
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules7)

results7 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results7 <- na.omit(results7)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut7 <- merge(results7, support7[, c("SNP", "chr", "position")])

#8
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules8)

results8 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results8 <- na.omit(results8)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut8 <- merge(results8, support8[, c("SNP", "chr", "position")])

#9
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules9)

results9 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results9 <- na.omit(results9)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut9 <- merge(results9, support9[, c("SNP", "chr", "position")])

#10
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules10)

results10 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results10 <- na.omit(results10)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut10 <- merge(results10, support10[, c("SNP", "chr", "position")])

#11
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules11)

results11 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results11 <- na.omit(results11)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut11 <- merge(results11, support11[, c("SNP", "chr", "position")])

#12
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules12)

results12 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results12 <- na.omit(results12)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut12 <- merge(results12, support12[, c("SNP", "chr", "position")])

#13
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules13)

results13 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results13 <- na.omit(results13)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut13 <- merge(results13, support13[, c("SNP", "chr", "position")])

#14
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules14)

results14 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results14 <- na.omit(results14)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut14 <- merge(results14, support14[, c("SNP", "chr", "position")])

#15
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules15)

results15 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results15 <- na.omit(results15)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut15 <- merge(results15, support15[, c("SNP", "chr", "position")])

#16
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules16)

results16 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results16 <- na.omit(results16)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut16 <- merge(results16, support16[, c("SNP", "chr", "position")])

#17
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules17)

results17 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results17 <- na.omit(results17)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut17 <- merge(results17, support17[, c("SNP", "chr", "position")])

#18
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules18)

results18 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results18 <- na.omit(results18)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut18 <- merge(results18, support18[, c("SNP", "chr", "position")])

#19
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules19)

results19 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results19 <- na.omit(results19)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut19 <- merge(results19, support19[, c("SNP", "chr", "position")])

#20
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules20)

results20 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results20 <- na.omit(results20)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut20 <- merge(results20, support20[, c("SNP", "chr", "position")])

#21
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules21)

results21 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results21 <- na.omit(results21)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut21 <- merge(results21, support21[, c("SNP", "chr", "position")])

#22
imp <- snp.rhs.tests(phenotype ~ gender + age + ethnicity + birth_country + education + Negative_life_events_witnessed + Negative_life_events_experienced +
                       + alcohol_use + arthritis + chronic_back_orneck_pain + frequent_or_severe_headaches + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
                       pc7 + pc8 + pc9 + pc10, family = "binomial", data = phenodata, snp.data = target, 
                     rules = rules22)

results22 <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results22 <- na.omit(results22)

#The results for SNPs that were successfully fit are then combined with the chromosome and position information by merging the results table with the support information, and adding a column for chromosome.
imputeOut22 <- merge(results22, support22[, c("SNP", "chr", "position")])


#combine the imputed files
imputeOut <- rbind(imputeOut1,imputeOut2,imputeOut3,imputeOut4,imputeOut5,imputeOut6,imputeOut7,imputeOut8,imputeOut9,imputeOut10,imputeOut11,imputeOut12,imputeOut13,imputeOut14,imputeOut15,imputeOut16,imputeOut17,imputeOut18,imputeOut19,imputeOut20,imputeOut21,imputeOut22)

write.csv(imputeOut, "impute.csv", row.names = FALSE)

#Post Analytic visualisation and genomic interrogation

#genomic interrogation
# having generated and fit both typed and imputed genotypes. Let's combine the results, and isolate just those SNPs in our region of interest.
#loading the typed SNP results generated by the GWAA function
# Read in GWAS output that was produced by GWAA function

GWASout <- read.table("GWAA.txt", header = TRUE, colClasses = c("character", 
                                                                rep("numeric", 4)))

GWASout$Neg_logP <- -log10(GWASout$p.value)

GWASout <- merge(GWASout, genoBim[, c("SNP", "chr", "position")])

GWASout <- arrange(GWASout, -Neg_logP)

head(GWASout)

#combining the tables of typed and imputed genotypes into a single table.
GWASout$type <- "typed"

imputeOut$Neg_logP <- -log10(imputeOut$p.value)

imputeOut <- arrange(imputeOut, -Neg_logP)

head(imputeOut)

GWAScomb<-rbind.fill(GWASout, imputeOut) #using rbind.fill() function to fill in missing columns with NAs

head(GWAScomb)
tail(GWAScomb)

GWAScomb$chr <- as.numeric(GWAScomb$chr)
GWAScomb <- GWAScomb[order(GWAScomb$chr),]
GWAScomb$type <- GWAScomb$type %>% replace_na('imputed')

GWAScomb <- GWAScomb[!is.na(GWAScomb$chr),]
nrow(GWAScomb)
str(GWAScomb)

#identify the significant snps
significant_snps <- GWAScomb[GWAScomb$p.value<=5e-08,]
nrow(significant_snps)

print(head(significant_snps))
snps <- list(significant_snps$SNP)

#post analytic visualization
#Manhattan plots- creating a manhattan plot to visualize GWA significant results by chromosome location
png(filename = "manhattan_plot.png",width = 8*600, height = 8*600, res = 350, pointsize = 25 ,bg = "white")
print(manhattan(GWAScomb, chr="chromosome", bp="position", snp="SNP", p="p-value" , main = "Manhattan Plot of snps", col = c("darkblue", "lightblue"), annotatePval = 5e-10))
dev.off()

#q-q plot
png(filename = "qq_plot.png",width = 5*600, height = 5*600, res = 300, pointsize = 12 ,bg = "white")
qq(manhattan_plot$p.value)
dev.off()

#preparing for FUMA
names(manhattan_plot)[names(manhattan_plot) == 'p.value'] <- 'P-value'
names(manhattan_plot)[names(manhattan_plot) == 'chr'] <- 'Chromosome'

#saving the data
write.csv(GWAScomb, "GWAS.csv", row.names = FALSE)

#in the terminal
#sed 's/"//g' GWAS.csv | sed 's/,/ /g' > GWAS.txt

#check the significant snps
significant_snps %>% 
  arrange(desc(p.value)) 
