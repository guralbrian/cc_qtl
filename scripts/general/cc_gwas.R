############### Time for the GWAS ##########################
library(ggplot2)
library(statgenGWAS)
library(graphics)
library(gridExtra)
library(gridGraphics)

tertiary <- read.csv("data/processed/snp_data/strain_mda_snps_tertiary.csv")

#Make the genotype/marker files file
dup   <- as.numeric(row.names(tertiary[which(duplicated(tertiary$marker)),]))#find duplicate rs ids
clean <- tertiary[-c(dup),] #remove them
clean <- clean[!is.na(clean$marker),] #remove na ids
row.names(clean)<- clean$marker

marker <- t(clean[,-c(4:5)])
geno <- clean[1:5]


#make the phenotype file
pheno <- read.csv("data/processed/phenotype_data/phenotype_iso_small.csv")
colnames(pheno)[1] <- "genotype"
row.names(pheno) <- pheno$genotype
pheno_list <- split(x = pheno[c("genotype", "Inh.A", "Inh.B", "Inh.C", "Inh.D", "Inh.E", "Inh.F", "Inh.G", "Inh.H", "BW.day.0", "BW.day.28", "TH", "LV", "RV", "LA","RA","Lung","Liver","Adrenal",             
                                "THW.by.BW.0","LVW.by.BW.0","RVW.by.BW.0","LAW.by.BW.0","RAW.by.BW.0","LuW.by.BW.0","LiW.by.BW.0","AdrW.by.BW.0")],
                    f = pheno[["Mouse"]])



#### pheno list test ####

#! there was an instance of a "Yes", instead of a "yes"
#! It lead to pheno_list having two main lists for "yes" and "Yes"
#! Is there a reason for the case difference or was it a fluke?

#!! these two parts should go to the data prep script
colnames(pheno)[1] <- "genotype"
pheno$Mouse[which(pheno$Mouse != "NA")] <- "yes"
row.names(pheno) <- pheno$genotype


mouse <- which(colnames(pheno) == "Mouse") # get location of mouse column

pheno_list <- split(x = pheno[-mouse],
                    f = pheno[["Mouse"]]) # rework dataframe into nested list


#make the GWAS input
gDataDrops1 <- createGData(geno = marker, map = geno) |>
                  createGData(pheno = pheno_list)

#! Is there a reason we can't do this? 
#gDataDrops2 <- createGData(geno = marker, map = geno, pheno = pheno_list) 


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

#  brian's plots 

#### fdrtools ####
generate_plots <- function(gwas_data, trait_name, title) {
  # Create Manhattan plot
  manhattan_plot <- plot(gwas_data, plotType = "manhattan", trait = trait_name, yThr = 4.2, title = title, lod = 3) +
    theme_bw() + theme(legend.position = "none")
  
  # Create QQ plot
  qq_plot <- plot(gwas_data, plotType = "qq", trait = trait_name, title = title) +
    theme_bw() + theme(legend.position = "none")
  
  return(list(manhattan_plot, qq_plot))
}


traits <- c("LVW.by.BW.0", "LiW.by.BW.0", "Inh.A", "Inh.B")
titles <- c("Normalized Left Ventricle Weight", "Normalized Liver Weight", "Percent Inheritance from Founder A", "Percent Inheritance from Founder B")



plots_list <- lapply(seq_along(traits), function(i) {
  generate_plots(GWASDrops1, traits[i], titles[i])
})

all_plots <- do.call(patchwork::wrap_plots, c(unlist(plots_list, recursive = FALSE), ncol = 2))
