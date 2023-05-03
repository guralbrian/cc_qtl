library(dplyr)
library(statgenGWAS)
library(rJava)
library(ggplot2)
library(dMod)
library(fdrtool)

#### Read in Data ####
s.mda.snp <- read.csv("D:/CC_Full/snp_data/Strain_MDA_SNPs.csv")
s.geno.allele <- read.csv("D:/CC_Full/snp_data/Strain_Genotypes.csv")
sac <- read_excel("phenotype_data/Strain Averages CTRL.xlsx")
sai <- read_excel("phenotype_data/Strain Averages ISO.xlsx")

#### Remove the unnecessary rows
s.mda.snp <- s.mda.snp[which(s.mda.snp$chr != "X" & s.mda.snp$chr != "Y" & s.mda.snp$chr!= "MT"),]
s.mda.snp[s.mda.snp == ""] <- NA

#### Identify major and minor alleles ####
input = t(s.mda.snp)
get_major_minor <- function(cur_col) {
  cur_col[cur_col == ""] <- NA
  cur_allele = cur_col[-c(1:3)]
  rank = as.data.frame(table(cur_allele))
  rank = rank[order(rank$Freq, decreasing = T),]
  allele.a = as.character(rank[1,1])
  allele.b = as.character(rank[2,1])
  return(c(allele.a, allele.b))
}

major_minor_alleles <- apply(input, 2, get_major_minor)
major_minor_alleles <- t(major_minor_alleles)  # Transpose the matrix to match the dimensions of s.mda.snp

output <- cbind(s.mda.snp[,c(1:3)], 
                Allele_A = major_minor_alleles[, 1], 
                Allele_B = major_minor_alleles[, 2], 
                s.mda.snp[,c(4:length(s.mda.snp))])

write.csv(output, "/media/raulab/88CA-7528/Strain_MDA_SNPs_with_alleles.csv", row.names = F)

#### Convert to major/minor alleles ####
alleles <- output
for(i in 1:nrow(alleles)){
  cur_row = alleles[i,6:74]
  cur_maj = as.character(alleles[i,4])
  cur_min = as.character(alleles[i,5])
  cur_row[cur_row == cur_maj] = "A"
  cur_row[cur_row == cur_min] = "B"
  alleles[i,6:74] = cur_row
  if(runif(1)<.005){
    print(paste0(i/nrow(output)*100,"% done!"))}
}
write.csv(alleles, "/media/raulab/88CA-7528/Strain_MDA_SNPs_as_alleles.csv", row.names = F)

#### Conver to binary ####
binary <- alleles[6:74]
binary[binary == "A"] <- 1
binary[binary == "B"] <- 0
binary[binary == "H"] <- 0.5
binary <- cbind(alleles[1:5], binary)
write.csv(binary, "D:/CC_Full/snp_data/Strain_MDA_SNPs_as_binary.csv", row.names = F)

####It wants weird 0, 1, 2 instead of 1,0,0.5 ####
tertiary <- alleles[6:74]
tertiary[tertiary == "A"] <- 0
tertiary[tertiary == "B"] <- 2
tertiary[tertiary == "H"] <- 1
tertiary <- cbind(alleles[1:5], tertiary)

write.csv(tertiary, "D:/CC_Full/snp_data/Strain_MDA_SNPs_as_012.csv", row.names = F)
