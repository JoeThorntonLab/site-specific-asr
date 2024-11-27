## This script measures the accuracy of the different substitution models by
## obtaining the identity percentage of the reconstructed sequences by each model 
## with the true ancestral sequence

## Author: Ricardo Muñiz Trejo
## Code developed in conjunction with ChatGPT, a language model developed by 
## OpenAI (https://openai.com/). 

## Date: 7-6-24

# phylotools contains the function read.phylip() which will be used to read the
# true ancestral sequences from the sim_ancestral.phy files
library(phylotools)

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

# Create an empty data frame to store the results
accuracy_table <- data.frame(model = character(),
                             rep = integer(),
                             node = integer(),
                             identity = numeric())

# specify which repetitions contain the data needed
working_directories <- directories


# loop over the all the repetitions of the simulations
for (i in 1:length(working_directories)) {
  
  rep_dir <- paste0("rep", working_directories[i])
  
  # set the path to the sim_ancestral.phy file in the rep directory
  sim_ancestral_file <- file.path(parent_dir,rep_dir, "sim_ancestral.phy")
  
  # read in the reference sequences from the sim_ancestral.phy file
  ref_seqs <- read.phylip(sim_ancestral_file)
  
  # set the reference sequences in their own vector
  
  true_anc_5 <- ref_seqs$seq.text[1]
  true_anc_6 <- ref_seqs$seq.text[2]
  
  # loop over the analysis directories
  for (j in 1:length(folders)) {
    
    abbreviation <- abbreviations[j]
    
    # set paths to node 5 and node 6 files in the analysis directory
    # file naming is different between PAML and DMSPhyloAA
    # PAML files end with *_sim_ancestor.txt
    # DMSPhyloAA files are simply names 2.txt and 5.txt
    # 2.txt is analogous to node_5__sim_ancestor.txt and 
    # 5.txt is analogous to node_6___sim_ancestor.txt
    
    node_5_file <- file.path(parent_dir, rep_dir, folders[j], ifelse(folders[j] %in% c("paml_best", "paml_poisson"), "node_5__sim_ancestor.txt", "ASR/2.txt"))
    node_6_file <- file.path(parent_dir, rep_dir, folders[j], ifelse(folders[j] %in% c("paml_best", "paml_poisson"), "node_6___sim_ancestor.txt", "ASR/5.txt"))
    
    # read in the sequences from the node 5 and node 6 files
    node_5_seq <- readLines(node_5_file, warn = FALSE)
    node_6_seq <- readLines(node_6_file, warn = FALSE)
    
    # Calculate and save identity to node 5
    identity <- calc_identity(true_anc_5, node_5_seq)
    accuracy_table[nrow(accuracy_table) + 1, ] <- c(abbreviation, working_directories[i], 5, identity)
    
    # Calculate and save identity to node 6
    identity <- calc_identity(true_anc_6, node_6_seq)
    accuracy_table[nrow(accuracy_table) + 1, ] <- c(abbreviation, working_directories[i], 6, identity)
      
  }
}

# add an identifier of the branch length and the protein system used for 
# the simulations
accuracy_table$blens <- rep(simulated_branch_length, nrow(accuracy_table))
accuracy_table$system <- rep(dataset, nrow(accuracy_table))

# save data into a csv file
write.csv(accuracy_table, file = "accuracy_table.csv", row.names = FALSE) #Print identity table

### Print session info
print("-----------------------------------------------------------------------")
print("ASR accuracy data collection completed!")
cat("\n")
sessionInfo()
cat("\n")
print("-----------------------------------------------------------------------")
cat("\n")
cat("\n")