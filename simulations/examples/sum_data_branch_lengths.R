## This script collects the total tree length estimated by the substitution 
## models for every repetition and it also gives the mean branch length for the
## four taxa tree that was used for the simulations

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


## The extract_blens_func() takes a tree in newick format, eliminates the puntuation
## marks, extract the numerical values of every branch length and sums them to
## to give the total tree length

extract_blens_func <- function(phylogram){  #To calculate the mean branch lengths
  phylogram <- gsub("t\\d:", "", phylogram) # Remove taxa names
  phylogram <- gsub("\\:|\\,", " ", phylogram) # Remove sep marks and introduce spaces
  phylogram <- gsub("\\(|\\)|\\;", "", phylogram) # Remove parentheses and semi colon
  nums <- as.numeric(scan(text = phylogram, what = "double", sep = " ", quiet = TRUE))
  nums <- signif(nums, digits = 6)
  return(sum(nums))
}

# set the path to the parent directory and folder names
parent_dir <- getwd()
folders <- c("paml_best", "paml_poisson", "DMSPhyloAA_RBD")
abbreviations <- c("best", "poisson", "RBD")

# specify which repetitions contain the data needed
working_directories <- directories

# Create an empty data frame to store the results
recorded_blens <- data.frame(model = character(),
                             rep = integer(),
                             total_blens = numeric(),
                             mean_blens = numeric())


# PAML 4.9j automatically stores the total tree length in the mlc.out file, so
# the info can be extracted directly. DMSPhyloAA only saves the tree in a newick
# file in phylogram.txt, so the extract_blens_func is applied to extract the data
# in the correct format.
files_with_trees = c('mlc.out','phylogram.txt')


# loop over the all the repetitions of the simulations
for (i in 1:length(folders)) {
  model <- folders[i]
  abbreviation <- abbreviations[i]
  
  for (j in 1:length(working_directories)) {
    rep_dir <- paste0("rep", working_directories[j])
    file_path <- file.path(parent_dir, rep_dir, model)
    
    # search for all mlc.out and phylogram.txt files
    files <- list.files(file_path, pattern = paste0(files_with_trees, collapse="|"), 
                        full.names = TRUE, recursive = TRUE)
    
    for (k in 1:length(files)) {
      
      if (grepl("mlc.out", files[k]) == TRUE) {
        blens_file <- readLines(files[k])
        blens_line <- grep("^tree length = ", blens_file)
        
        if (length(blens_line) > 0) { # make sure to grab the correct line
          blens <- as.numeric(gsub("tree length = ", "", blens_file[blens_line]))
          mean_blens <- blens/5 # get the average branch length value
          recorded_blens[nrow(recorded_blens) + 1, ] <- c(abbreviation, working_directories[j], blens, mean_blens)
        }
      } else {
        blens_file <- readLines(files[k])
        
        if (length(blens_line) > 0) {
          blens <- extract_blens_func(blens_file) # apply extract_blens_func to get total tree length
          mean_blens <- blens/5 # get the average branch length value
          
          
          recorded_blens[nrow(recorded_blens) + 1, ] <- c(abbreviation, working_directories[j], blens, mean_blens)
        }
      }
    }
    
  }
}

# add an identifier of the branch length and the protein system used for 
# the simulations
recorded_blens$blens <- rep(simulated_branch_length, nrow(recorded_blens))
recorded_blens$system <- rep(dataset, nrow(recorded_blens))

# Save branch length data into a csv
write.csv(recorded_blens, file = "blens_estimation.csv", row.names = FALSE)

### Print session info
print("-----------------------------------------------------------------------")
print("Branch length estimation data collection completed!")
cat("\n")
sessionInfo()
cat("\n")
print("-----------------------------------------------------------------------")
cat("\n")
cat("\n")