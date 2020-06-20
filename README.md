# **Welcome to the DIALOGUE!**

DIALOGUE is a dimensionality reduction method that uses cross-cell-type associations to identify multicellular programs (MCPs) and map the cell transcriptome as a function of its environment. Given single-cell data, it combines penalized matrix decomposition with multilevel modeling to identify generalizable MCPs and examines their association with specific phenotypes of interest.

<img src="https://github.com/livnatje/DIALOGUE/blob/master/Images/DIALOGUEoverview.png" width=400 />

# **Requirements**

* R (tested in R version 3.4.0).
* R libraries: lme4, lmerTest, PMA, plyr, matrixStats, psych, stringi, RColorBrewer, unikn, reshape2, ggplot2

# **Quick start**

To install DIALOGUE you can either use ```devtools::install_github("DIALOGUE",your_user_name)``` or just download the [package](https://singlecell.broadinstitute.org/single_cell/study/SCP958/dialogue#study-download) and use ```devtools::install("DIALOGUE")```

Here you can find a simple [step-by-step example](https://github.com/livnatje/DIALOGUE/wiki/Step-by-step-example).

All you need for the **input** is the single-cell transcriptomes of different cell types, usually together with some more compact representation (e.g., PCs). The **output** will be multicellular programs (MCPs) of co-regulated genes across the different cell types. Each MCP consists of cell-type-specific gene subsets.

For specific cell-cell "interactions" you can run the pairwise version, using the data of two cell types of interest as input. DIALOGUE can also identify MCPs that span multiple cell types (see **_Jerby-Arnon and Regev BioRxiv 2020_** for examples on 5 and 6 cell types). 

# Citation

Jerby-Arnon and Regev _**Mapping multicellular programs from single-cell transcriptomes**_.

