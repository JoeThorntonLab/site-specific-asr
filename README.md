# site-specific-asr

This repository contains the scripts, data, and results found in the paper:

Muñiz-Trejo, R, Y. Park & J. W. Thornton (2024). Robustness of ancestral sequence reconstruction to among-site and among-lineage evolutionary heterogeneity

This main directory contains a copy of the **DMSPhyloAA.py** script which was used
to inform the site-specific models and obtain the ancestral reconstructions. The
script **simulate_alignment.R**, which was used to produce the site-specific
simulations of protein sequence ecolution, is also here. 

**empirical/** contains the all the data and results for the 
analysis of the real alignment data of the three protein families (i. e. steroid
receptor (SR) DNA-binding domain, sarbecoviruses recognition binding domain (RBD), and 
hemagglutinin (HA)) studied in this work. 

**simulations/** contains the model parameter data used to obtain
the simulations of each protein family, as well as the input four-taxa tree. Tables
that summarize the data from the original simulation and the scripts used to produce
those files are also located there. Additionally, ten examples of simulations using
the models parameters of SR at 0.8 substitutions/site are provided. 

DMSPhyloAA.py syntax to obtain site-specific recontructions:

`python3 /DMSPhyloAA.py -a /alignment.phy -t tree.txt -p /dms_data.csv -o <path-to-results-directory> --ASR 1`

DMSPhyloAA.py syntax to obtain site-homogeneous recontructions:

`python3 /DMSPhyloAA.py -a /alignment.phy -t tree.txt -p /dms_data.csv -o <path-to-results-directory> --ASR 1 --useSiteHomogeneousModel 1 --AAR /model_rates.txt --AAPI /model_frequencies.txt`

Syntax for simulate_alignment.R

`Rscript /simulate_alignment.R /inferred_parameters.txt /dms_data.csv /tree.txt <branch-length-to-simulate>`

