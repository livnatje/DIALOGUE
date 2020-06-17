# **Welcome to the DIALOGUE!**

DIALOGUE is a dimensionality reduction approach that uses cross-cell-type associations to identify multicellular programs and map the cell transcriptome as a function of its environment. Given single-cell data, it combines penalized matrix decomposition with multilevel modeling, to identify generalizable multicellular programs and examines their association with specific phenotypes of interest. By doing so it also robustly recovers spatial information and can characterize the cell environment based only based on its transcriptome.

# **Requirements**

* R (tested in R version 3.4.0).
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
``k`` denotes the number of multicellular programs (MCPs) that will be identified. The different MCPs are not correlated with one another, and the cross-cell-type correlations observed within an MCP usually decreases with k, such that the first few MCPs depict most of the multicellular co-expression. DIALOGUE will always find the same MCPs or a subset of them, no matter which k is used.


You can also reproduce the colon/IBD multicellular program reported in our paper using the following code 
```
rA<-readRDS(system.file("extdata", "IBD.data.rds", package = "DIALOGUE"))
R<-DIALOGUE.run(rA = rA,main = "IBD",k = 2,results.dir = "DIALOGUE.results/")
```

See ```?DIALOGUE::DIALOGUE.run``` for more information.

# General notes

DIALOGUE will identify multicellular programs, such that each program will have different cell-type-specific components. It will generate figures to depict the association between the different cell-type-components of each multicellular program.

# Citation

Jerby-Arnon and Regev _**Mapping multicellular programs from single-cell transcriptomes**_.

