## This script obtains the ancestral reconstruction by maximum parsimony of all
## the simulations to determine the strength of phylogenetic signal. It also 
## compares the ASR accuracy of the maximum likelihood reconstructions with 
## their phylogenetic signal. It produces a total of four output files.

## Author: Ricardo Muñiz Trejo
## Code developed in conjunction with ChatGPT, a language model developed by 
## OpenAI (https://openai.com/). 

## Date: 7-6-24

# phylotools contains the function read.phylip() which will be used to read the
# true ancestral sequences from the sim_ancestral.phy files. phangorn is used
# to obtain the parsimony reconstructions. data.table is used to optimize
# data storage and time
library(phylotools)
library(phangorn)
library(data.table)

# specify the branch length used for simulations when calling the script,
# the directories with complete analyses in 'c(x:y)' format, and the protein
# system used

args <- commandArgs(trailingOnly = TRUE)
blens <- as.numeric(args[1])
directories_str <- args[2]
directories <- eval(parse(text = directories_str))
dataset <- args[3]

# Define function to get the matches between a reconstructed sequence and the
# true sequence
get_mismatch <- function(seq1, seq2) {
  matches <- charToRaw(seq1) == charToRaw(seq2)
  return(matches)
}

# Define a function to check if there's at least one occurrence of 1 in a row
check_parsimony <- function(table) {
  ifelse(1 %in% table, "unambiguous", "ambiguous")
}

# Define a function to process a matrix and generate the mp_sequence vector
get_mp_reconstruction <- function(table) {
  mp_sequence <- c()  # Initialize mp_sequence vector
  
  # Loop over each row of the matrix
  for (i in 1:nrow(table)) {
    # Find the maximum value in the current row
    max_val <- max(table[i,])
    
    # Check if the maximum value is 1
    if (max_val == 1) {
      max_index <- which(table[i,] == max_val)[1] 
      mp_sequence <- c(mp_sequence, colnames(table)[max_index])
    } else {
      mp_sequence <- c(mp_sequence, "X")
    }
  }
  
  return(mp_sequence)  # Return the mp_sequence vector
}

AA <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')


# set the path to the parent directory
parent_dir <- getwd()
folders <- c("paml_best", "DMSPhyloAA_RBD")

# Create an empty data frame to store the results
results_dt <- data.table(rep = integer(),
                             node = integer(),
                             identical = logical(),
                             correct = logical(),
                             correct_ss = logical(),
                             parsimony = character(),
                             mp_correct = logical(),
                             mp_is_map = logical(),
                             mp_is_map_ss = logical())

# specify which repetitions contain the data needed
working_directories <- directories


# loop over the all the repetitions of the simulations
for (i in 1:length(working_directories)) {
  
  rep_dir <- paste0("rep", working_directories[i])
  
  ## Read true ancestral sequences ##
  
  # set the path to the sim_ancestral.phy file in the rep directory
  sim_ancestral_file <- file.path(parent_dir,rep_dir, "sim_ancestral.phy")
  
  # read in the reference sequences from the sim_ancestral.phy file
  ref_seqs <- read.phylip(sim_ancestral_file)
  
  # set the reference sequences in their own vector
  
  true_anc_5 <- ref_seqs$seq.text[1]
  true_anc_6 <- ref_seqs$seq.text[2]
  
  ## Generate most parsimonious reconstruction (MPR) and determine ##
  ## ambiguity in the sequence                                     ##
  
  # read alignment into a phyDat object
  alignment <- read.phyDat(
    file = file.path(parent_dir,rep_dir, "sim_alignment.phy"),
    format="phylip",
    type = "AA"
  )
  # read tree
  tree <- read.tree(file.path(parent_dir,rep_dir, "best_ml_tree.txt"))
  
  # obtain MP reconstructions
  anc.mpr <- ancestral.pars(tree, alignment, "MPR")
  
  # rename marix columns to amino acid states
  colnames(anc.mpr[[5]]) <- AA; colnames(anc.mpr[[6]]) <- AA
  
  # extract indexes of site patterns
  site_pattern_index <- attributes(anc.mpr)$index
  
  # obtain node 5 mpr and ambiguity list
  parsimony_site_pattern <- apply(anc.mpr[[5]], 1, check_parsimony)
  state_site_pattern <- get_mp_reconstruction(anc.mpr[[5]])
  
  node_5_ambiguity <- parsimony_site_pattern[site_pattern_index]
  node_5_mpr <- state_site_pattern[site_pattern_index]
  node_5_mpr <- paste(node_5_mpr, collapse = "")
  
  
  # obtain node 6 mpr and ambiguity list
  parsimony_site_pattern <- apply(anc.mpr[[6]], 1, check_parsimony)
  state_site_pattern <- get_mp_reconstruction(anc.mpr[[6]])
  
  node_6_ambiguity <- parsimony_site_pattern[site_pattern_index]
  node_6_mpr <- state_site_pattern[site_pattern_index]
  node_6_mpr <- paste(node_6_mpr, collapse = "")
  
  
  ## Read reconstructed sequences ##
  SH_node_5_file <- readLines(file.path(parent_dir, rep_dir, folders[1], 
                                        "node_5__sim_ancestor.txt"), warn = F)
  SH_node_6_file <- readLines(file.path(parent_dir, rep_dir, folders[1], 
                                        "node_6___sim_ancestor.txt"), warn = F)
  
  
  SS_node_5_file <- readLines(file.path(parent_dir, rep_dir, folders[2], 
                                        "ASR/2.txt"), warn = F)
  SS_node_6_file <- readLines(file.path(parent_dir, rep_dir, folders[2], 
                                        "ASR/5.txt"), warn = F)
  
  
  ## Compare sequences ##
  
  # are SH and SS reconstructions identical?
  identical_5 <- get_mismatch(SH_node_5_file, SS_node_5_file)
  identical_6 <- get_mismatch(SH_node_6_file, SS_node_6_file)
  
  # are SH reconstructions correct?
  correct_5 <- get_mismatch(SH_node_5_file, true_anc_5)
  correct_6 <- get_mismatch(SH_node_6_file, true_anc_6)
  
  # are SS reconstructions correct?
  correct_5_ss <- get_mismatch(SS_node_5_file, true_anc_5)
  correct_6_ss <- get_mismatch(SS_node_6_file, true_anc_6)
  
  
  # is the MP state the correct reconstruction?
  mp_correct_5 <- get_mismatch(true_anc_5, node_5_mpr)
  mp_correct_6 <- get_mismatch(true_anc_6, node_6_mpr)
  
  # is the MP state assigned as the MAP one in SH?
  mp_map_5 <- get_mismatch(SH_node_5_file, node_5_mpr)
  mp_map_6 <- get_mismatch(SH_node_6_file, node_6_mpr)
  
  # is the MP state assigned as the MAP one in SS?
  mp_map_5_ss <- get_mismatch(SS_node_5_file, node_5_mpr)
  mp_map_6_ss <- get_mismatch(SS_node_6_file, node_6_mpr)
  
  ## Populate data table ##
  
  # Add results from node 5
  results_dt <- rbind(results_dt, data.table(
    rep = working_directories[i],
    node = 5,
    identical = identical_5,
    correct = correct_5,
    correct_ss = correct_5_ss,
    parsimony = node_5_ambiguity,
    mp_correct = mp_correct_5,
    mp_is_map = mp_map_5,
    mp_is_map_ss = mp_map_5_ss
  ))
  
  # Add results from node 6
  results_dt <- rbind(results_dt, data.table(
    rep = working_directories[i],
    node = 6,
    identical = identical_6,
    correct = correct_6,
    correct_ss = correct_6_ss,
    parsimony = node_6_ambiguity,
    mp_correct = mp_correct_6,
    mp_is_map = mp_map_6,
    mp_is_map_ss = mp_map_6_ss
  ))
  
}


#############################################################
#### Comparison of errors between models ####################
#############################################################

all_sites <- nrow(results_dt)

identical_wrong <- nrow(results_dt[identical == TRUE & 
                                     correct == FALSE])/all_sites

disagreement_wrong <- nrow(results_dt[identical == FALSE & 
                                        correct == FALSE & 
                                        correct_ss == FALSE])/all_sites


disagreement_ss_right <- nrow(results_dt[identical == FALSE & 
                                           correct == FALSE & 
                                           correct_ss == TRUE])/all_sites

disagreement_sh_right <- nrow(results_dt[identical == FALSE & 
                                           correct == TRUE & 
                                           correct_ss == FALSE])/all_sites

error_dist <- data.frame(category = c("Same error", "Different error", 
                                      "Only SH correct", "Only SS correct"), 
                         proportion = c(identical_wrong, 
                                        disagreement_wrong, 
                                        disagreement_sh_right, 
                                        disagreement_ss_right),
                         system = dataset,
                         blens = blens
)

write.csv(error_dist, file = "error_comparison.csv", 
          row.names = FALSE)


#############################################################
#### What is the probability of error given the type of #####
#### parsimoniuos reconstruction?                       ##### 
#############################################################

# Probability of error when signal is true
true_signal_sites <- nrow(results_dt[parsimony == "unambiguous" & 
                                       mp_correct == TRUE])

error_true_signal_sh <- nrow(results_dt[parsimony == "unambiguous" & 
                                          mp_correct == TRUE & 
                                          correct == FALSE])/true_signal_sites

error_true_signal_ss <- nrow(results_dt[parsimony == "unambiguous" & 
                                          mp_correct == TRUE & 
                                          correct_ss == FALSE])/true_signal_sites


# Probability of error when signal is misleading

# Misleading signal happens when the unambiguous MP state is wrong

misleading_signal_sites <- nrow(results_dt[parsimony == "unambiguous" & 
                                             mp_correct == FALSE])

error_misleading_signal_sh <- nrow(results_dt[parsimony == "unambiguous" & 
                                                mp_correct == FALSE & 
                                                correct == FALSE])/misleading_signal_sites

error_misleading_signal_ss <- nrow(results_dt[parsimony == "unambiguous" & 
                                                mp_correct == FALSE & 
                                                correct_ss == FALSE])/misleading_signal_sites


# Probability of error when there phylogenetic signal is absent

absent_sites <- nrow(results_dt[parsimony == "ambiguous"])

error_absent_sh <- nrow(results_dt[parsimony == "ambiguous" & 
                                        correct == FALSE])/absent_sites

error_absent_ss <- nrow(results_dt[parsimony == "ambiguous" & 
                                        correct_ss == FALSE])/absent_sites

prob_error_type <- data.frame(category = rep(c("True", "Misleading", "Absent"), times = 2), 
                              proportion = c(error_true_signal_sh, 
                                             error_misleading_signal_sh, 
                                             error_absent_sh,
                                             error_true_signal_ss, 
                                             error_misleading_signal_ss, 
                                             error_absent_ss),
                              reconstruction = c("SH", "SH", "SH", "SS", "SS", "SS"),
                              system = dataset,
                              blens = blens
)


write.csv(prob_error_type, file = "probability_error_by_signal.csv", 
          row.names = FALSE)

##################################################################
#### How often is the maximum a posteriori (MAP) state the    ####
#### same as the most parsimonious (MP) state when the signal ####
#### is misleading?                                           ####
##################################################################


misleading_sites <- nrow(results_dt[parsimony == "unambiguous" & 
                                      mp_correct == FALSE])

misleading_mp_sh <- nrow(results_dt[parsimony == "unambiguous" & 
                                      mp_correct == FALSE & 
                                      mp_is_map == TRUE])/misleading_sites
misleading_mp_ss <- nrow(results_dt[parsimony == "unambiguous" & 
                                      mp_correct == FALSE & 
                                      mp_is_map_ss == TRUE])/misleading_sites


misleading_data <- data.frame(proportion = c(misleading_mp_sh, 
                                             misleading_mp_ss),
                              reconstruction = c("SH", "SS"),
                              system = dataset,
                              blens = blens)

write.csv(misleading_data, file = "parsimony_at_misleading.csv", 
          row.names = FALSE)


###################################################################
#### What is proportion of the different types of parsimonious ####
#### reconstructions per simulation?                           ####
###################################################################

prop_true <- nrow(results_dt[parsimony == "unambiguous" & 
                               mp_correct == TRUE])/all_sites
prop_misleading <- nrow(results_dt[parsimony == "unambiguous" & 
                                     mp_correct == FALSE])/all_sites
prop_absent <- nrow(results_dt[parsimony == "ambiguous"])/all_sites

sites_by_signal <- data.frame(category = c("True", "Misleading", "Absent"),
                              proportion = c(prop_true,
                                             prop_misleading,
                                             prop_absent),
                              system = dataset,
                              blens = blens)

write.csv(sites_by_signal, file = "parsimony_reconstruction.csv", 
          row.names = FALSE)


### Print session info
print("-----------------------------------------------------------------------")
print("Parsimony reconstruction data collection completed!")
cat("\n")
sessionInfo()
cat("\n")
print("-----------------------------------------------------------------------")
cat("\n")
cat("\n")
