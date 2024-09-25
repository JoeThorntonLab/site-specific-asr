
# Script for simulating an alignment using an experimental phylogenetic model.


# PACKAGES REQUIRED.
library(phangorn)


# FUNCTIONS AND CONSTANTS.

base <- c('A', 'C', 'G', 'T')
AA <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

# Genetic code.
gcode <- data.frame(codon = c('TTT', 'CTT', 'ATT', 'GTT', 'TCT', 'CCT', 'ACT', 'GCT', 'TAT', 'CAT', 'AAT', 'GAT', 'TGT', 'CGT', 'AGT', 'GGT',
                              'TTC', 'CTC', 'ATC', 'GTC', 'TCC', 'CCC', 'ACC', 'GCC', 'TAC', 'CAC', 'AAC', 'GAC', 'TGC', 'CGC', 'AGC', 'GGC',
                              'TTA', 'CTA', 'ATA', 'GTA', 'TCA', 'CCA', 'ACA', 'GCA', 'TAA', 'CAA', 'AAA', 'GAA', 'TGA', 'CGA', 'AGA', 'GGA',
                              'TTG', 'CTG', 'ATG', 'GTG', 'TCG', 'CCG', 'ACG', 'GCG', 'TAG', 'CAG', 'AAG', 'GAG', 'TGG', 'CGG', 'AGG', 'GGG'),
                    AA = c('F', 'L', 'I', 'V', 'S', 'P', 'T', 'A', 'Y', 'H', 'N', 'D', 'C', 'R', 'S', 'G',
                           'F', 'L', 'I', 'V', 'S', 'P', 'T', 'A', 'Y', 'H', 'N', 'D', 'C', 'R', 'S', 'G',
                           'L', 'L', 'I', 'V', 'S', 'P', 'T', 'A', '*', 'Q', 'K', 'E', '*', 'R', 'R', 'G',
                           'L', 'L', 'M', 'V', 'S', 'P', 'T', 'A', '*', 'Q', 'K', 'E', 'W', 'R', 'R', 'G'))

# Derive the exchangeability matrices `R` (list), equilibrium amino acid frequencies `PI` (list),
# relative rate of substitution at each site `site_rate`, number of discrete-gamma categories `ngammacat`, and the shape parameter `alpha`.
# Global variables of these names are created.
derive_model <- function() {
  
  # 1. Calculating the mutation rates.
  
  # Inferred nucleotide exchangeabilities.
  DNAR <- matrix(0, 4L, 4L, dimnames = list(base, base))
  DNAR[cbind(c(2L, 3L, 3L, 4L, 4L, 4L), c(1L, 1L, 2L, 1L, 2L, 3L))] <- 
    as.numeric(strsplit(trimws(gsub('\\].*', '', gsub('.*\\[', '', inferred_param[2L]))), '\\s+')[[1L]])
  DNAR <- DNAR + t(DNAR)
  
  # Inferred nucleotide frequencies.
  DNAPI <- as.numeric(strsplit(trimws(gsub('\\].*', '', gsub('.*\\[', '', inferred_param[3L]))), '\\s+')[[1L]])
  names(DNAPI) <- base
  
  mu <- matrix(0, length(AA), length(AA), dimnames = list(AA, AA))
  
  for(i in seq_along(AA)) {
    for(j in seq_along(AA)) {
      
      if(i == j) next 
      
      codon_i <- gcode$codon[gcode$AA == AA[i]]
      codon_i_freq <- vapply(codon_i, function(codon) prod(DNAPI[strsplit(codon, '')[[1L]]]), 0)
      
      codon_j <- gcode$codon[gcode$AA == AA[j]]
      codon_j_freq <- vapply(codon_j, function(codon) prod(DNAPI[strsplit(codon, '')[[1L]]]), 0)
      
      for(k in seq_along(codon_i)) {
        for(l in seq_along(codon_j)) {
          
          diff_site <- which(strsplit(codon_i[k], '')[[1L]] != strsplit(codon_j[l], '')[[1L]])
          if(length(diff_site) != 1L) next # Only single-mutation changes are allowed
          
          base_k <- substr(codon_i[k], diff_site, diff_site)
          base_l <- substr(codon_j[l], diff_site, diff_site)
          
          mu[i, j] <- mu[i, j] + DNAR[base_k, base_l] * DNAPI[base_l] * codon_i_freq[k] / sum(codon_i_freq)
        }
      }
    }
  }
  
  
  # 2. Calculating the equilibrium amino acid frequencies.
  
  # Equilibrium amino acid frequencies under no selection.
  PI_mutation_only <- vapply(AA, function(k) {
    
    sum(vapply(gcode$codon[gcode$AA == k], function(codon) prod(DNAPI[strsplit(codon, '')[[1L]]]), 1))
    
  }, 1)
  
  # Inferred fitness function (logistic function).
  param <- as.numeric(strsplit(trimws(gsub('\\].*', '', gsub('.*\\[', '', inferred_param[5L]))), '\\s+')[[1L]])
  
  # List of equilibrium amino acid frequencies for each site.
  PI <<- lapply(1L:nsite, function(i) {
    
    phi <- dms[i, ]
    growth_rate <- param[1L] / (1 + exp(-param[2L] * (phi - param[3L])))
    PI_unnorm <- PI_mutation_only * exp(2 * growth_rate)
    PI_unnorm / sum(PI_unnorm)
  })
  
  
  # 3. Calculating the unscaled transition rates.
  
  Q <- lapply(1L:nsite, function(i) {
    
    phi <- dms[i, ]
    growth_rate <- param[1L] / (1 + exp(-param[2L] * (phi - param[3L])))
    s <- t(outer(growth_rate, growth_rate, `-`))
    fix_prob <- s / (1 - exp(-2 * s))
    fix_prob[abs(s) < 1e-10] <- 0.5 # Dealing with numerical problems arising due to very small values of s.
    Qi <- mu * fix_prob
    diag(Qi) <- -apply(Qi, 1L, sum, na.rm = TRUE)
    Qi
  })
  
  
  # 4. Scaling the transition rates.
  
  # Total rate of substitution at each site.
  site_rate <- vapply(1L:nsite, function(i) sum(PI[[i]] * -diag(Q[[i]])), 1)
  Q <- lapply(Q, function(x) x / mean(site_rate))
  site_rate <<- site_rate / mean(site_rate)
  
  
  # 5. Deriving the exchangeability matrix.
  
  R <<- lapply(1L:nsite, function(i) {Ri <- t(Q[[i]]) / PI[[i]]; diag(Ri) <- NA_real_; Ri})
  
  
  # 6. Setting the discrete gamma-distributed among-site rate variation model.
  
  alpha <<- as.numeric(strsplit(inferred_param[4L], split = ' ')[[1L]][3L])
  ngammacat <- strsplit(inferred_param[4L], split = ' ')[[1L]][7L]
  ngammacat <<- as.integer(substr(ngammacat, 1L, nchar(ngammacat) - 1L))
}


# COMMANDLINE INPUTS

# Directories to 1) 'inferred_parameters.txt', 2) DMS data in csv format, and 3) phylogeny, and a branch-length scaler.
input <- commandArgs(trailingOnly = TRUE)

inferred_param <- readLines(input[1L])

dms <- as.matrix(read.csv(input[2L]))
nsite <- nrow(dms)

tree <- unroot(read.tree(input[3L]))
tree$edge.length <- tree$edge.length * as.numeric(input[4L]) # Scaling branch lengths


# RUN.

derive_model()

# Simulation.

alignment <- matrix(NA_character_, length(tree$tip.label), nsite)
ancestral <- matrix(NA_character_, tree$Nnode, nsite)

gamma_rates <- discrete.gamma(alpha, ngammacat)

for(i in 1L:nsite) {
  
  # Simulation for site i.
  
  sim <- simSeq(tree, l = 1L, Q = R[[i]], bf = PI[[i]], type = 'AA',
                rate = sample(gamma_rates, 1L) * site_rate[i], ancestral = TRUE)
  
  alignment[, i] <- AA[unname(unlist(sim))[1L:length(tree$tip.label)]]
  ancestral[, i] <- AA[unname(unlist(sim))[(length(tree$tip.label) + 1L):(length(tree$tip.label) + tree$Nnode)]]
}

row.names(alignment) <- names(sim)[1L:length(tree$tip.label)]
row.names(ancestral) <- names(sim)[(length(tree$tip.label) + 1L):(length(tree$tip.label) + tree$Nnode)]


# Exporting alignments as phylip file.

write(paste0(' ', nrow(alignment), ' ', nsite), 'sim_alignment.phy')
for(i in 1L:nrow(alignment))
  write(paste(rownames(alignment)[i], paste(alignment[i, ], collapse = '')), 'sim_alignment.phy', append = TRUE)

write(paste0(' ', nrow(ancestral), ' ', nsite), 'sim_ancestral.phy')
for(i in 1L:nrow(ancestral))
  write(paste(rownames(ancestral)[i], paste(ancestral[i, ], collapse = '')), 'sim_ancestral.phy', append = TRUE)

