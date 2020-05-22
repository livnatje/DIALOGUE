# **Welcome to the DIALOGUE!**

DIALOGUE is a method for mapping multicellular configurations based on single-cell data. By combining penalized matrix decomposition with multilevel modeling, DIALOGUE robustly recovers spatial information, characterizes the cell environment only based on its transcriptome, and identifies generalizable multicellular programs/configurations, including those that underly different phenotypes of interest in health and disease.

# **Requirements**

* R (tested in R version 3.4.0 (2017-04-21) -- "You Stupid Darkness").
* R libraries: lme4, lmerTest, PMA, plyr, matrixStats, psych, stringi, RColorBrewer, unikn, reshape2, ggplot2

# **Quick start**

To install you can either use ```devtools::install_github("DIALOGUE",your_user_name)``` or ```devtools::install("DIALOGUE")```

The data for testing is provided in the 
[Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP958/dialogue#study-download)
(make sure to download and uncompress the package to the DIALOGUE directory).

To run a toy example, download the toy example data
```
rA<-readRDS(system.file("extdata", "toy.example.rds", package = "DIALOGUE"))
```
Find multicellular programs:
```
R<-DIALOGUE.run(rA = rA,main = "toy.example",k = 2,results.dir = "DIALOGUE.results/")
```

You can also reproduce the colon/IBD multicellular program reported in our paper using the following code 
```
rA<-readRDS(system.file("extdata", "IBD.data.rds", package = "DIALOGUE"))
R<-DIALOGUE.run(rA = rA,main = "IBD",k = 2,results.dir = "DIALOGUE.results/")
```

See ```?DIALOGUE::DIALOGUE.run``` for more information.

# General notes

DIALOGUE will identify multicellular programs, such that each program will have different cell-type-specific components. It will generate figures to depict the association between the different cell-type-components of each multicellular program.

# Citation

Jerby-Arnon and Regev _**Mapping multicellular configurations using single-cell data**_.

