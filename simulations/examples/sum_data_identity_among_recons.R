## This script measures how similar are the reconstruction made by different
## substitution models in a pairwise way. Unlike the identity with the true
## ancestral sequence, this script's purpose is to get a sense of how dissimilar
## are the reconstructions among models of the simulated sequences

## Author: Ricardo Muñiz Trejo
## Code developed in conjunction with ChatGPT, a language model developed by 
## OpenAI (https://openai.com/). 

## Date: 7-6-24

# phylotools contains the function read.phylip() which will be used to read the
# true ancestral sequences from the sim_ancestral.phy files
# data.table is used to reduce compute time
library(phylotools)
library(data.table)

# specify the branch length used for simulations when calling the script
# and the directories with complete analyses in 'c(x:y)' format

args <- commandArgs(trailingOnly = TRUE)
simulated_branch_length <- as.numeric(args[1])
directories_str <- args[2]
directories <- eval(parse(text = directories_str))
dataset <- args[3]

# Define function to calculate identity percentage between two sequences of the
# same length
calc_identity <- function(seq1, seq2) {
  matches <- sum(charToRaw(seq1) == charToRaw(seq2))
  mismatches <- sum(charToRaw(seq1) != charToRaw(seq2))
  identity <- (matches / (matches + mismatches)) * 100
  return(identity)
}

# set the path to the parent directory
parent_dir <- getwd()
folders <- c("paml_best", "paml_poisson", "DMSPhyloAA_RBD")
abbreviations <- c("best", "poisson", "RBD")

# specify which repetitions contain the data needed
working_directories <- directories


# Define the directories and files
dirs <- paste0("rep", working_directories)

# define the directories of the different substitution models, the abbreviations,
# and the file names depending on whether PAML or DMSPhyloAA was used for ancestral 
# sequence reconstruction
model_files <- list(
  list(dir = "paml_best", abbreviation = "best", file1 = "node_5__sim_ancestor.txt", file2 = "node_6___sim_ancestor.txt"),
  list(dir = "paml_poisson", abbreviation = "poisson", file1 = "node_5__sim_ancestor.txt", file2 = "node_6___sim_ancestor.txt"),
  list(dir = "DMSPhyloAA_RBD/ASR", abbreviation = "RBD", file1 = "2.txt", file2 = "5.txt")
)

# Initialize an empty data table to store results
results_dt <- data.table(
  rep = integer(),
  node = integer(),
  model1 = character(),
  model2 = character(),
  identity = numeric()
)

# Loop over directories and files to perform pairwise comparisons
for (dir in 1:length(dirs)) {
  model1 <- 1
  while (model1 < length(model_files)) {
    model2 <- model1 + 1
    while (model2 <= length(model_files)) {
      
      # Read in the sequences from the two files that correspond to ancestor 5
      file1_path <- file.path(parent_dir, dirs[dir], model_files[[model1]]$dir, model_files[[model1]]$file1)
      file2_path <- file.path(parent_dir, dirs[dir], model_files[[model2]]$dir, model_files[[model2]]$file1)
      seq1 <- readLines(file1_path)
      seq2 <- readLines(file2_path)
      # Calculate the identity and store results in the data table
      identity <- calc_identity(seq1, seq2)
      results_dt <- rbind(results_dt, data.table(
        rep = dir,
        node = 5,
        model1 = model_files[[model1]]$abbreviation,
        model2 = model_files[[model2]]$abbreviation,
        identity = identity
      ))
      
      # Read in the sequences from the two files that correspond to ancestor 6
      file1_path <- file.path(parent_dir, dirs[dir], model_files[[model1]]$dir, model_files[[model1]]$file2)
      file2_path <- file.path(parent_dir, dirs[dir], model_files[[model2]]$dir, model_files[[model2]]$file2)
      seq1 <- readLines(file1_path)
      seq2 <- readLines(file2_path)
      # Calculate the identity and store results in the data table
      identity <- calc_identity(seq1, seq2)
      results_dt <- rbind(results_dt, data.table(
        rep = dir,
        node = 6,
        model1 = model_files[[model1]]$abbreviation,
        model2 = model_files[[model2]]$abbreviation,
        identity = identity
      ))
      
      model2 <- model2 + 1
    }
    model1 <- model1 + 1
  }
}

# add an identifier of the branch length and the protein system used for 
# the simulations
results_dt[, blens := simulated_branch_length]
results_dt[, system := dataset]

# save table with summarized values
write.csv(results_dt, file = "identity_among_models.csv", row.names = FALSE)

### Print session info
print("-----------------------------------------------------------------------")
print("Identity among reconstructions data collection completed!")
cat("\n")
sessionInfo()
cat("\n")
print("-----------------------------------------------------------------------")
cat("\n")
cat("\n")
