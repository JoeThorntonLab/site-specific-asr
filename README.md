# site-specific-asr

This repository contains the scripts, data, and results found in the paper:

Muñiz Trejo, R, Y. Park & J. W. Thornton (2024). Robustness of ancestral sequence 
reconstruction to among-site model heterogeneity and epistasis

This main directory contains a copy of the **DMSPhyloAA.py** script which was used
to inform the site-specific models and obtain the ancestral reconstructions. The
script **simulate_alignment.R**, which was used to produce the site-specific
simulations of protein sequence ecolution, is also here. 

**empirical/** contains the all the data and results for the 
analysis of the real alignment data of the three protein families (i. e. steroid
receptor (SR) DNA-binding domain, sarbecoviruses recognition binding domain (RBD), and 
hemagglutinin (HA)) studied in this work. The script for the analysis of these data
is included.

**simulations/** contains the model parameter data used to obtain
the simulations of each protein family, as well as the input four-taxa tree. Tables
that summarize the data from the original simulation and the scripts used to produce
those files are also located there. Additionally, ten examples of simulations using
the models parameters of SR at 0.8 substitutions/site are provided. The script for 
the analysis of these data is included.

simulate_alignment.R

Rscript <path to script>/simulate_alignment.R 
<path to data>/**_infered_parameters.txt
<path to data>/dms_dataset.csv
<path to data>/tree.txt
<insert branch length value>

`Rscript /site-specific-asr/simulate_alignment.R /site-specific-asr/simulations/data/rbd_inferred_parameters.txt /site-specific-asr/simulations/data/rbd_dms_sars.csv /site-specific-asr/simulations/data/tree.txt 0.4`

