############### Time for the GWAS ##########################
#Make the genotype/marker files file
dup <- as.numeric(row.names(tertiary[which(duplicated(tertiary$marker)),]))#find duplicate rs ids
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