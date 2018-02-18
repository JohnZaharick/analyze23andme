#23andMe provides all customers with a raw data file containing ~610,000 SNPs.
#A fraction of these have documented phenotypes. They can be found on snpedia.com
#This script uses the SNPediaR package to scrape snpedia.com and identify SNPs that have summaries,
#resulting in ~200-300 SNPs for a user to read about instead of searching all 610,000.
#There is the risk of false negatives in that a SNP may have no summary, but still have a description.
#This script is for exploration of raw 23andMe data and not an exhaustive genetic diagnotistic tool.

install.packages("devtools")
library (devtools)
install_github ("genometra/SNPediaR")
library (SNPediaR)
install.packages("dplyr")
library(dplyr)

#manually open the raw data .txt file from 23andMe in a spreadsheet editor and remove the intial comments
#the final version must contain four columns labelled rsid, chromosome, position, genotype
#save this as a .csv, then load it into R
fileToLoad <- "genome_John_Doe.csv"
#the file our list of SNPs, genotypes, and phenotypes will end up in
saveFile <- "finalSNPlist.csv"

#download most recent list of snps on snpedia using SNPediaR package
snpediaDB <- getCategoryElements (category = "Is_a_snp")
snpediaDB <- tolower(snpediaDB)

#load the personal snpedia file
SNPs <- read.csv(fileToLoad, header=TRUE)

#reduce 23andMe list to only SNPs that are in snpedia
#this will take you from ~610,000 SNPs to ~23,000. A bit more managable.
reducedSNPs <- SNPs %>% filter(SNPs$rsid %in% snpediaDB)

#build a function to extract the summary of a SNP using SNPediaR package
#some snpedia pages will cause getPages() to throw an error due to formatting
#try functions will keep getSummary() from throwing an error
getSummary <- function(name){
  tryCatch(snp <- getPages(titles = name),error=function(cond) {print(name)})
  snp2 <- "NA"
  try(snp2 <- extractTags(snp[[1]], tags = "Summary"), silent=T)
  if (class(snp2) == "try-error") {
    snp2 <- "NA"
  }
  return(snp2[[1]])
}

#build function to loop through SNPs, calling snpedia API
cycleList <- function(arg){
  exportList <- vector("list", length(arg)*2)
  n <- 1
  for (i in arg){
    B<-getSummary(i)
    if(!is.na(B)){
      exportList [[n]] <- i
      n <- n + 1
      exportList [[n]] <- B
      n <- n + 1
    }
  }
  exportList
}

#store only the vector of SNPs for scraping snpedia.com's API
snpsToQuery <- reducedSNPs$rsid

#check every SNP on snpedia.com, saving ones that have summaries
queriedSNPs <- cycleList(snpsToQuery)
queriedSNPs <- unlist(queriedSNPs)
SNPsAndPhenotypes = matrix(queriedSNPs, ncol=2, byrow=TRUE)
colnames(SNPsAndPhenotypes) <- c("rsid","phenotype")
#convert to data frame so we can use $ to select column
SNPsAndPhenotypes <- as.data.frame(SNPsAndPhenotypes)

#filter SNPs once again, this time keeping only ones with summaries
reducedSNPs <- reducedSNPs %>% filter(reducedSNPs$rsid %in% SNPsAndPhenotypes$rsid)
#append summaries
reducedSNPs$phenotype <- SNPsAndPhenotypes$phenotype
#create a new .csv file of our saved results
#manually search the identified SNPs on snpedia.com for the best information
#do not rely on only the downloaded summaries.  Many are vague.
write.csv(reducedSNPs, file = saveFile, row.names=FALSE)