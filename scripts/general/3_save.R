# Save the files, plots, and session info

# Create a folder with the current date in the "data/results" directory
current_date <- format(Sys.Date(), "%Y-%m-%d")
dir_path <- file.path("data", "results", current_date)

# Function to check if the directory already exists
# Useful if making >1 output per day
MakelDirUnique = function(prefix){
  if(!file.exists(dir_path)) {return(dir_path)}
  i <- 1
  repeat {
    f = paste(dir_path, i, sep="_")
    if(!file.exists(f)) {return(f)}
    i <- i + 1
  }
}

dir_path <- MakelDirUnique(dir_path)
dir.create(dir_path, showWarnings = FALSE)

# Save the plot to the newly created folder
plot_file_path <- file.path(dir_path, paste0("plot_", current_date, ".png"))
ggsave(filename = plot_file_path, plot = all_plots, width = 16, height = 20, dpi = "print")

# Save the open scripts
SaveCopyScripts <- function(source = "scripts/general",
                            destination = dir_path){
  rstudioapi::documentSaveAll()
  # Copy the saved general scripts to the destination folder
  scripts <- list.files(source)
  for(i in scripts){
    script_dest <- file.path(dir_path, i)
    script <- file.path(source, i)
    file.copy(script, script_dest)
  }
}

#SaveCopyScripts() 

# Save the session info to the newly created folder
writeLines(capture.output(sessionInfo()), 
           file.path(dir_path, paste0("session_info_", current_date, ".txt")))