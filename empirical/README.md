# Empirical data

In **data/** you can find all the initial files needed to perform the ancestral
reconstructions. 

- \*_alignment.phy Contains the amino acid alignment for the protein system
- \*_tree.txt Contains the phylogeny used in this work for the protein system
- \*_dms_\*.csv Contains the functional data values 
- *_parsimony_gaps.csv. Contains a table that summarizes the results of parsimony-based ASR of the recodified alignment to identify gaps. A for amino acid state, G for gap. Node numbering was modified to match the node numbering of DMSPhyloAA.py

- *_rates. Contains the exchangeabilities of the site-homogeneous model
- *_freqs. Contains the maximum likelihood estimates of the empirical frequencies of the alignment data of the protein. The estimates were obtained using RAxML 8

In **results/** you can find directories with the ancestral reconstructions and parameters of the 
Directories are named by the protein family (sr, rbd, ha), type of reconstruction (ss, sh), and the name of the model.

- AAPI.txt. Contains the estimated empirical frequencies per site. For SH reconstructions, all sites have the same numbers.
- inferred_parameters.txt All ML estimated of the parameters for the SS model fits are contained here. For the SH models 
- labeled_cladogram.txt
- per_site_llk.txt
- phylogram.txt
- state_labeled_trees.txt

In the subdirectory **ASR/** there are two types of files. 