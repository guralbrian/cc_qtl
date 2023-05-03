library(dplyr)
library(statgenGWAS)
library(rJava)
library(ggplot2)
library(dMod)
library(fdrtool)
library("readxl")
#### Read in Data ####
s.mda.snp <- read.csv("data/processed/snp_data/strain_mda_snps.csv")
s.geno.allele <- read.csv("data/processed/snp_data/strain_genotypes.csv")
sac <- read_excel("data/processed/phenotype_data/strain_averages_ctrl.xlsx")
sai <- read_excel("data/processed/phenotype_data/strain_averages_iso.xlsx")

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


#write.csv(output, "data/processed/snp_data/strain_mda_snps_alleles.csv", row.names = F)

#### Convert to major/minor alleles ####


# Function for the revised method
alleleConvert <- function(data, major = "A", minor = "B") {
  data %>%
    mutate(across(starts_with("CC"),
                  ~ case_when(
                    . == Allele_A ~ major,
                    . == Allele_B ~ minor,
                    TRUE ~ .
                  )))
}

alleles_ab <- alleleConvert(output, "A", "B")

#write.csv(allele_ab, "data/processed/snp_data/strain_mda_snps_ab.csv", row.names = F)

#### Conver to binary ####
binary <- alleles[6:74]
binary[binary == "A"] <- 1
binary[binary == "B"] <- 0
binary[binary == "H"] <- 0.5
binary <- cbind(alleles[1:5], binary)
write.csv(binary, "data/processed/snp_data/strain_mda_snps_binary.csv", row.names = F)

####It wants weird 0, 1, 2 instead of 1,0,0.5 ####
tertiary <- alleles[6:74]
tertiary[tertiary == "A"] <- 0
tertiary[tertiary == "B"] <- 2
tertiary[tertiary == "H"] <- 1
tertiary <- cbind(alleles[1:5], tertiary)

write.csv(tertiary, "data/processed/snp_data/strain_mda_snps_tertiary.csv", row.names = F)

############### Time for the GWAS ##########################


#Make the genotype/marker files file
dup   <- as.numeric(row.names(tertiary[which(duplicated(tertiary$marker)),]))#find duplicate rs ids
clean <- tertiary[-c(dup),] #remove them
clean <- clean[!is.na(clean$marker),] #remove na ids
row.names(clean)<- clean$marker

marker <- t(clean[,-c(4:5)])
geno <- clean[1:5]


#make the phenotype file
pheno <- read.csv("D:/CC_Full/phenotype_data/Phenotype_Iso_small.csv")
colnames(pheno)[1] <- "genotype"
row.names(pheno) <- pheno$genotype
pheno_list <- split(x = pheno[c("genotype", "Inh.A", "Inh.B", "Inh.C", "Inh.D", "Inh.E", "Inh.F", "Inh.G", "Inh.H", "BW.day.0", "BW.day.28", "TH", "LV", "RV", "LA","RA","Lung","Liver","Adrenal",             
                                "THW.by.BW.0","LVW.by.BW.0","RVW.by.BW.0","LAW.by.BW.0","RAW.by.BW.0","LuW.by.BW.0","LiW.by.BW.0","AdrW.by.BW.0")],
                    f = pheno[["Mouse"]])

#make the GWAS input
gDataDrops1 <- createGData(geno = marker, map = geno)
gDataDrops1 <- createGData(gData = gDataDrops1, pheno = pheno_list)
gDataDropsImputed1 <- codeMarkers(gData = gDataDrops1,
                                  nMissGeno = 0.05, 
                                  nMiss = 0.05, 
                                  MAF = 0.05,
                                  impute = T, 
                                  imputeType = "random", 
                                  verbose = TRUE)


#do the GWAS
GWASDrops1 <- runSingleTraitGwas(gData = gDataDropsImputed1,
                                 trials = c("yes"),
                                 traits = c("Inh.A","Inh.B", "Inh.C", "Inh.D", "Inh.E", "Inh.F", "Inh.G", "Inh.H", "BW.day.0", "BW.day.28", "TH", "LV", "RV", "LA","RA","Lung","Liver","Adrenal",
                                            "THW.by.BW.0","LVW.by.BW.0","RVW.by.BW.0","LAW.by.BW.0","RAW.by.BW.0","LuW.by.BW.0","LiW.by.BW.0","AdrW.by.BW.0"),
                                 GLSMethod = "multi",
                                 kinshipMethod ="IBS",
                                 thrType = "fixed",
                                 LODThr = 0)
#plot the GWAS
plot(GWASDrops1, plotType = "manhattan", trait = "LVW.by.BW.0", yThr = 4.2, title = "Normalized Left Ventricle Weight")+
  theme_bw() + theme(legend.position = "none")
plot(GWASDrops1, plotType = "qq", trait = "LVW.by.BW.0", title = "Normalized Left Ventricle Weight")+
  theme_bw() + theme(legend.position = "none")

plot(GWASDrops1, plotType = "manhattan", trait = "LiW.by.BW.0", yThr = 4.2, title = "Normalized Liver Weight")+
  theme_bw() + theme(legend.position = "none")
plot(GWASDrops1, plotType = "qq", trait = "LiW.by.BW.0", title = "Normalized Liver Weight")+
  theme_bw() + theme(legend.position = "none")


plot(GWASDrops1, plotType = "manhattan", trait = "Inh.A", yThr = 4.2, title = "Percent Inheritance from Founder A")+
  theme_bw() + theme(legend.position = "none")
plot(GWASDrops1, plotType = "qq", trait = "Inh.A", title = "Percent Inheritance from Founder A")+
  theme_bw() + theme(legend.position = "none")

plot(GWASDrops1, plotType = "manhattan", trait = "Inh.B", yThr = 4.2, title = "Percent Inheritance from Founder B")+
  theme_bw() + theme(legend.position = "none")
plot(GWASDrops1, plotType = "qq", trait = "Inh.B", title = "Percent Inheritance from Founder B")+
  theme_bw() + theme(legend.position = "none")


#### fdrtools ####
