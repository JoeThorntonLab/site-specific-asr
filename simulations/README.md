# Simulated data

We generated a large amount of simulations for each system and each branch length
condition (1000 replicates for SR, 500 for RBD, and 100 for HA). The files for 
producing the simulations for each dataset is available in **data/**.

The large number of replicates per branch length amounts to thousands of files 
and a lot of storage space. Since uploading that much info is unfeasible, we 
summarized the simulation results in tables with the info needed to replicate the
analyses shown in the manuscript. These tables are stored in the **results/**
subdirectory. Some of the naming conventions between the manuscript and the
summarized table remain, but the underlying data remains the same.

In **examples/** we copied the raw data of ten replicates of the simulations using
SS-SARS at 0.8 substitutions/site. We also included the file for the reconstructions 
with the generative model (DMSPhyloAA-RBD), the best-fitting site-homogeneous model
(paml_best), and the Poisson model (paml_poisson). In addition, the scripts used
to produce the summary tables are located in here and work with the ten examples
using the following syntax:

`Rscript /examples/script_name*.R 0.8 "c(1:10)" RBD`

**DISCLAIMER:** The performance and memory usage of the scripts to summarize the data
was optimized using ChatGPT-4o. In particular, the implementation of the data.table 
package was suggested by the artificial intelligence and it help in reducing 
computing time of the large amounts of data

 

