Welcome to the DIALOGUE!

DIALOGUE is a method for mapping multicellular configurations based on single-cell data. By combining penalized matrix decomposition with multilevel modeling, DIALOGUE robustly recovers spatial information, characterizes the cell environment only based on its transcriptome, and identifies generalizable multicellular programs/configurations, including those that underly different phenotypes of interest in health and disease.

# **Requirements**

* R (tested in R version 3.4.0 (2017-04-21) -- "You Stupid Darkness").
* R libraries: lme4, lmerTest, PMA, plyr, matrixStats, psych, stringi, RColorBrewer, unikn, reshape2, ggplot2

# **Quick start**

To install you can either use ```devtools::install_github("DIALOGUE",your_user_name)``` or ```devtools::install("DIALOGUE")```

To run a toy example
```
rA<-readRDS(system.file("extdata", "toy.example.rds", package = "DIALOGUE"))
Find multicellular programs:
R<-DIALOGUE.run(rA = rA,main = "toy.example",k = 2,results.dir = "~/Desktop/DIALOGUE.results/")
```
# General notes

DIALOGUE will identify multicellular progams, such that each program has a cell-type-specific component. It will generate figures to depict the association between the different cell-type-components of each multicellular program.

# Citation

Jerby-Arnon L and Regev _**Mapping multicellular configurations using single-cell data**_.
