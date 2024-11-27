## This script reads the posterior probability (PP) tables of each repetition,
## extract the max PP value per site (row) and calculates the mean posterior
## probability per residue and store the value in a data frame

## Author: Ricardo Muñiz Trejo
## Code developed in conjunction with ChatGPT, a language model developed by 
## OpenAI (https://openai.com/). 

## Date: 7-6-24

# specify the branch length used for simulations when calling the script
# and the directories with complete analyses in 'c(x:y)' format

args <- commandArgs(trailingOnly = TRUE)
simulated_branch_length <- as.numeric(args[1])
directories_str <- args[2]
directories <- eval(parse(text = directories_str))
dataset <- args[3]

# set the path to the parent directory and folder names
parent_dir <- getwd()
folders <- c("paml_best", "paml_poisson", "DMSPhyloAA_RBD")
abbreviations <- c("best", "poisson", "RBD")

# Create an empty data frame to store the results
mean_pp <- data.frame(model = character(),
                      rep = integer(),
                      mean_max = numeric())

# specify which repetitions contain the data needed
working_directories <- directories

# loop over the all the repetitions of the simulations
for (i in 1:length(folders)) {
  model <- folders[i]
  abbreviation <- abbreviations[i]
  
  for (j in 1:length(working_directories)) {
    rep_dir <- paste0("rep", working_directories[j])
    file_path <- file.path(parent_dir, rep_dir, model)
    files <- list.files(file_path, pattern = "_PP.txt", full.names = TRUE, recursive = TRUE)
    
    for (k in 1:length(files)) {
      # as PAML and DMSPhyloAA PP tables have slightly different format, if the
      # file name ends in _sim_ancestor_PP.txt, values are separated with a space,
      # otherwise, use a comma as a separator.
      pp_table <- read.table(files[k], sep = ifelse(grepl("_sim_ancestor_PP.txt", files[k]), " ", ","), header = TRUE)
      max_values <- apply(pp_table, 1, max) # get max pp value per row
      mean_max <- mean(max_values) # calculate the mean per pp table
      mean_pp[nrow(mean_pp) + 1, ] <- c(abbreviation, working_directories[j], mean_max)
    }
  }
}

# add an identifier of the branch length and the protein system used for 
# the simulations
mean_pp$blens <- rep(simulated_branch_length, nrow(mean_pp))
mean_pp$system <- rep(dataset, nrow(mean_pp))

# Save data into a csv
write.csv(mean_pp, file = "mean_posterior_probability.csv", row.names = FALSE)

### Print session info
print("-----------------------------------------------------------------------")
print("Mean MAP state posterior probability data collection completed!")
cat("\n")
sessionInfo()
cat("\n")
print("-----------------------------------------------------------------------")
cat("\n")
cat("\n")