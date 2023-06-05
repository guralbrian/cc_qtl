library(statgenGWAS)
library(dplyr)
library(patchwork)
library(ggplot2)

# Load data
pheno <- read.csv("data/processed/phenotype_data/phenotype_iso_clean.csv", row.names = 1)
tertiary <- read.csv("data/processed/snp_data/strain_mda_snps_tertiary.csv", row.names = 1)

# Format data
marker <- t(tertiary[,-c(4:5)])
geno <- tertiary[1:5]

# Make the phenotype file
mouse <- which(colnames(pheno) == "Mouse") 
pheno_list <- split(x = pheno[-mouse],
                    f = pheno[["Mouse"]])

# Make the GWAS input
gDataDrops1 <- createGData(geno = marker, map = geno) |>
  createGData(pheno = pheno_list)
gDataDropsImputed1 <- codeMarkers(gData = gDataDrops1,
                                  nMissGeno = 0.05, 
                                  nMiss = 0.05, 
                                  MAF = 0.05,
                                  impute = T, 
                                  imputeType = "random", 
                                  verbose = TRUE)


# Run GWAS
GWASDrops1 <- runSingleTraitGwas(gData = gDataDropsImputed1,
                                 trials = c("yes"),
                                 traits = colnames(pheno_list$yes)[-1],
                                 GLSMethod = "multi",
                                 kinshipMethod ="IBS",
                                 thrType = "fixed",
                                 LODThr = 0)
# Plot the GWAS

# Function to make GWAS plots
GeneratePlots <- function(gwas_data, trait_name, title) {
  # Create Manhattan plot
  manhattan_plot <- plot(gwas_data, plotType = "manhattan", trait = trait_name, yThr = 4.2, title = title, lod = 3) +
    theme_bw() + theme(legend.position = "none")
  
  # Create QQ plot
  qq_plot <- plot(gwas_data, plotType = "qq", trait = trait_name, title = title) +
    theme_bw() + theme(legend.position = "none")
  
  return(list(manhattan_plot, qq_plot))
}

# Arguements for GeneratePlots()
traits <- c("LVW.by.BW.0", "LiW.by.BW.0", "Inh.A", "Inh.B")
titles <- c("Normalized Left Ventricle Weight", "Normalized Liver Weight", 
            "Percent Inheritance from Founder A", "Percent Inheritance from Founder B")


# Make the plots, store in list
plots_list <- lapply(seq_along(traits), function(i) {
  GeneratePlots(GWASDrops1, traits[i], titles[i])
})

# Wrap all of the plots together
all_plots <- do.call(patchwork::wrap_plots, c(unlist(plots_list, recursive = FALSE), ncol = 2))
print(all_plots)