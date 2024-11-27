## This purpose of this script is to parse through all repetition of simulations
## for all used models to obtain the posterior probability values of every 
## amino acid state for all sites in each node. Then it checks if the reconstructed 
## state was correctly reconstructed by comparing the values with the true ancestral
## sequence. Finally, it bins the posterior probability values in 5% increments
## and calculates the proportion of times that the state was reconstructed correctly

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

## I obtained the function get_midpoint from: 
## https://stackoverflow.com/questions/22312207/how-to-assign-cut-range-midpoints-in-r
## (Consulted 3/23/2023 7:35 PM)

get_midpoint <- function(cut_label) {
  mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
} ## This function returns the midpoint of a range of values obatined by applyind
## 'cut()' to a numeric vector


# set the path to the parent directory
parent_dir <- getwd()
folders <- c("paml_best", "paml_poisson", "DMSPhyloAA_RBD")
abbreviations <- c("best", "poisson", "RBD")


# create a vector with the parameter names
names <- c("model", "rep", "node", "site", "true_seq", "recons_state", "posterior_prob", "matches")


# create an empty data table to store the parameter values

matches_dt <- data.table(matrix(ncol = length(names), nrow = 0))
setnames(matches_dt, names)

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
  
  # number of sites 
  site <- 1:nchar(true_anc_5)
  
  # loop over the analysis directories
  for (j in 1:length(folders)) {
    
    abbreviation <- abbreviations[j]
    
    # set paths to node 5 and node 6 files in the analysis directory
    node_5_table_path <- file.path(parent_dir, rep_dir, folders[j], ifelse(folders[j] %in% c("paml_best", "paml_poisson"), "node_5__sim_ancestor_PP.txt", "ASR/2_PP.txt"))
    node_6_table_path <- file.path(parent_dir, rep_dir, folders[j], ifelse(folders[j] %in% c("paml_best", "paml_poisson"), "node_6___sim_ancestor_PP.txt", "ASR/5_PP.txt"))
    
    # read in the sequences for the true node 5 sequence and PP table for each repetition and each model
    pp_table <- read.csv(node_5_table_path, sep = ifelse(folders[j] %in% c("paml_best", "paml_poisson"), " ", ","))
    true_anc <- strsplit(true_anc_5, split = "")[[1]][1:nrow(pp_table)]
    
    # reorganize PP table in two columns to compare the posterior probability values
    # of each state with the true ancestral sequence
    pp_table <- tidyr::gather(pp_table, key = "state", value = "probability")
    true_anc <- rep(true_anc, times=20)
    
    matches <- pp_table$state == true_anc
    
    matches_dt <- rbind(matches_dt, data.table(model = abbreviation, 
                                               rep = working_directories[i],
                                               node = 5,
                                               site = site,
                                               true_seq = true_anc,
                                               recons_state = pp_table$state,
                                               posterior_prob = pp_table$probability,
                                               matches = matches))
    
    
    # repeat the previous analysis for node 6
    pp_table <- read.csv(node_6_table_path, sep = ifelse(folders[j] %in% c("paml_best", "paml_poisson"), " ", ","))
    true_anc <- strsplit(true_anc_6, split = "")[[1]][1:nrow(pp_table)]
    
    pp_table <- tidyr::gather(pp_table, key = "state", value = "probability")
    true_anc <- rep(true_anc, times=20)
    
    matches <- pp_table$state == true_anc
    
    matches_dt <- rbind(matches_dt, data.table(model = abbreviation, 
                                               rep = working_directories[i],
                                               node = 6,
                                               site = site,
                                               true_seq = true_anc,
                                               recons_state = pp_table$state,
                                               posterior_prob = pp_table$probability,
                                               matches = matches))
  }
  
  # filter out all cases in which the posterior probability is zero
  matches_dt <- subset(matches_dt, posterior_prob != 0)
  
}



# bin the posterior probability values in 5% increments
matches_dt$bin <- cut(matches_dt$posterior_prob, breaks = seq(-0.05,1.05,0.05))

# Apply get_midpoint() to get the midpoint of each bin
matches_dt$bin_midpoint <- sapply(as.character(matches_dt$bin), get_midpoint)

# summarise data to facilitate future data processing and data visualization
matches_dt <- matches_dt[, count := .N, by = .(bin_midpoint, model)]
matches_summary <- matches_dt[, .(prop_true = sum(matches), total_elements = .N), 
                              by = .(bin_midpoint, model)][, prop_correct := prop_true / total_elements][order(model, bin_midpoint)]

# add an identifier of the branch length and the protein system used for 
# the simulations
matches_summary[, blens := simulated_branch_length]
matches_summary[, system := dataset]

# save table with summarized values
write.csv(matches_summary, file = "pp_bins_summary.csv", row.names = FALSE)

### Print session info
print("-----------------------------------------------------------------------")
print("All states posterior probability data collection completed!")
cat("\n")
sessionInfo()
cat("\n")
print("-----------------------------------------------------------------------")
cat("\n")
cat("\n")