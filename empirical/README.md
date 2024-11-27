# Empirical data

In this directory you can find all the data required to performed the both the
site-specific and site-homogeneous reconstructions for the tree protein families.

In **data/** you can find all the initial files needed to perform the ancestral
reconstructions. 


- \*_alignment.phy Contains the amino acid alignment for the protein system
- \*_tree.txt Contains the phylogeny used in this work for the protein system
- \*_dms_\*.csv Contains the experimental data on the protein function upon amino acid mutation
- *_rates. Contains the exchangeabilities of the site-homogeneous model
- *_freqs. Contains the maximum likelihood estimates of the empirical frequencies of the alignment data of the protein. The estimates were obtained using RAxML 8.2

In addition, this directory contains some additional data like:

- *_parsimony_gaps.csv. Contains a table that summarizes the results of parsimony-based ASR of the recodified alignment to identify gaps obtained using Mesquite 3.70. A for amino acid state, G for gap. Node numbering was modified to match the node numbering of DMSPhyloAA.py
- *_parsimony.csv. Contains a table with the maximum parsimony reconstructions for all the nodes of each protein family obtained using Mesquite 3.70.
- *_phyml.txt. Contains the output of the estimation of rates of evolution for the alignments calculated by PhyML 3.0


In **results/** you can find the ancestral reconstructions and associated 
parameters of the for the fitting of the site-specific and site-homogeneous
models. Directories are named by the protein (sr, rbd, ha), 
type of reconstruction (ss, sh), and the name of the model.

Each subdirectory contains the following info:

- AAPI.txt. Contains the estimated empirical frequencies per site. 
- inferred_parameters.txt. All ML parameters estimates for the SS models.
- labeled_cladogram.txt. Ultrametric tree with assigned numbers for the nodes
- per_site_llk.txt. Maximum likelihood estimate for each site in the alignment
- phylogram.txt. Unlabelled tree with optimized branch lengths for the fitted model 
- state_labeled_trees.txt. Collection of trees with nodes labelled with the ML state. There is one tree per site. 
- Subdirectory **ASR/** which contains files with the ML reconstructions and a table with posterior probabilities (_PP.txt) of recontruction of all states at all sites

Sources of data:
Sequence alignment, tree, and DMS data for the steroid receptor DNA-binding 
domain (SR) was obtained from Park Y, Metzger BPH, Thornton JW. 2022. Epistatic 
drift causes gradual decay of predictability in protein evolution. Science 
376:823–830.

Sequence alignment, tree, and DMS data for the sarbecovirus recognition binding 
domain (RBD) was obtained from Starr TN, Greaney AJ, Hilton SK, Ellis D, 
Crawford KHD, Dingens AS, Navarro MJ, Bowen JE, Tortorici MA, Walls AC, et al. 
2020. Deep Mutational Scanning of SARS-CoV-2 Receptor Binding Domain Reveals 
Constraints on Folding and ACE2 Binding. Cell 182:1295-1310.e20.

Sequence alignment, tree, and DMS data for hemagglutinin (HA) was obtained from 
Hilton SK, Bloom JD. 2018. Modeling site-specific amino-acid preferences deepens
phylogenetic estimates of viral sequence divergence. Virus Evolution 4:vey033.
Although, we have to point out that the DMS for H1 was first described in 
Doud MB, Bloom JD. 2016. Accurate Measurement of the Effects of All Amino-Acid 
Mutations on Influenza Hemagglutinin. Viruses 8:155. The DMS for H3 was described
first in Lee JM, Huddleston J, Doud MB, Hooper KA, Wu NC, Bedford T, Bloom JD. 
2018. Deep mutational scanning of hemagglutinin helps predict evolutionary fates 
of human H3N2 influenza variants. Proceedings of the National Academy of Sciences 
115:E8276–E8285.

