library(dplyr)
library(statgenGWAS)
library(rJava)
library(ggplot2)
library(dMod)
library(readxl)

# Read in Data 
s_mda_snp <- read.csv("data/processed/snp_data/strain_mda_snps.csv")
s_geno_allele <- read.csv("data/processed/snp_data/strain_genotypes.csv")
sac <- readxl::read_excel("data/processed/phenotype_data/strain_averages_ctrl.xlsx")
sai <- readxl::read_excel("data/processed/phenotype_data/strain_averages_iso.xlsx")
pheno <- read.csv("data/processed/phenotype_data/phenotype_iso_small.csv")

# Remove the unnecessary rows
s_mda_snp <- s_mda_snp[c(which(s_mda_snp$chr != "X" & 
                               s_mda_snp$chr != "Y" & 
                               s_mda_snp$chr != "MT"&
                    duplicated(s_mda_snp$marker) == FALSE)), ] # Only autosomes, no duplicates 

s_mda_snp[s_mda_snp == ""] <- NA
s_mda_snp <- s_mda_snp[!is.na(s_mda_snp$marker),]

# Identify major and minor alleles
input <- t(s_mda_snp)

# Function to look at each row and find the most and second most frequent alleles. Returns a list of each, which should both be added as columns in a new dataframe with the prior snp matrix
GetMajorMinor <- function(cur_col) {
  cur_col[cur_col == ""] <- NA
  cur_allele = cur_col[-c(1:3)]
  rank = as.data.frame(table(cur_allele))
  rank = rank[order(rank$Freq, decreasing = T),]
  allele.a = as.character(rank[1,1])
  allele.b = as.character(rank[2,1])
  return(c(allele.a, allele.b))
}

major_minor_alleles <- apply(input, 2, GetMajorMinor)
major_minor_alleles <- t(major_minor_alleles)  # Transpose the matrix to match the dimensions of s_mda_snp

# Remake the allele matrix
output <- cbind(s_mda_snp[, c(1:3)], 
                Allele_A = major_minor_alleles[, 1], 
                Allele_B = major_minor_alleles[, 2], 
                s_mda_snp[, c(4:length(s_mda_snp))])


write.csv(output, "data/processed/snp_data/strain_mda_snps_alleles.csv", row.names = FALSE)

# Convert to major/minor alleles

# Function to convert alleles (A,C,T,G) to major and minor for each site (A,B)
AlleleConvert <- function(data, major = "A", minor = "B") {
  data |>
    dplyr::mutate(across(dplyr::starts_with("CC"),
                         ~ dplyr::case_when(
                           . == Allele_A ~ major,
                           . == Allele_B ~ minor,
                           TRUE ~ .
                         )))
}

alleles_ab <- AlleleConvert(output, "A", "B")

write.csv(alleles_ab, "data/processed/snp_data/strain_mda_snps_ab.csv", row.names = FALSE)

# Convert to binary
cc <- grep("^CC", colnames(output), value = TRUE) # Select strain allele columns
not_cc <- colnames(output[which(colnames(output) %in% cc == FALSE)]) # And everything else

binary <- output[cc]
binary[binary == "A"] <- 1
binary[binary == "B"] <- 0
binary[binary == "H"] <- 0.5
binary <- cbind(output[not_cc], binary)

write.csv(binary, "data/processed/snp_data/strain_mda_snps_binary.csv", row.names = FALSE)

# It wants weird 0, 1, 2 instead of 1, 0, 0.5
tertiary <- output[cc]
tertiary[tertiary == "A"] <- 0
tertiary[tertiary == "B"] <- 2
tertiary[tertiary == "H"] <- 1
tertiary <- cbind(output[not_cc], tertiary)

write.csv(tertiary, "data/processed/snp_data/strain_mda_snps_tertiary.csv", row.names = FALSE)

# Also format the phenotype data for later

colnames(pheno)[1] <- "genotype"
pheno$Mouse[which(pheno$Mouse != "NA")] <- "yes" # Makes all values the same case
row.names(pheno) <- pheno$genotype

write.csv(pheno, "data/processed/phenotype_data/phenotype_iso_small_clean.csv", row.names = TRUE)
