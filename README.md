# **Welcome to the DIALOGUE!**

DIALOGUE is a dimensionality reduction method that uses cross-cell-type associations to identify multicellular programs (MCPs) and map the cell transcriptome as a function of its environment. Given single-cell data, it combines penalized matrix decomposition with multilevel modeling to identify generalizable MCPs and examines their association with specific phenotypes of interest.

<img src="https://github.com/livnatje/DIALOGUE/blob/master/Images/Livnat_dialogue_rev3-02.png" width=600 />

# **Quick start**

To install DIALOGUE you can either use [```devtools::install_github(repo = "https://github.com/livnatje/DIALOGUE")```](https://www.rdocumentation.org/packages/devtools/versions/1.13.6/topics/install_github) or just download its R package and use ```devtools::install("DIALOGUE")```

The **input** consistes of single-cell transcriptomes of different cell types, usually together with a more compact representation (e.g., PCs). The **output** will be multicellular programs (MCPs) of co-regulated genes across the different cell types, their expression across the cells, and association with specific phenotype(s) of interest. Each MCP consists of multiple cell-type-specific gene subsets.

For specific cell-cell "interactions" you can run the pairwise version, using the data of two cell types of interest as input. DIALOGUE can also identify MCPs that span multiple cell types (as we show in our pre-print **_Jerby-Arnon and Regev bioRxiv 2020_**).

See the [tutorial](https://github.com/livnatje/DIALOGUE/wiki/Tutorial) for more details.

### **Requirements**

* R (tested in R version 3.4.0).
* R libraries: lme4, lmerTest, PMA, plyr, matrixStats, psych, stringi, RColorBrewer, unikn, reshape2, ggplot2, grid, beanplot, parallel

# Citation

[Jerby-Arnon & Regev. DIALOGUE maps multicellular programs in tissue from single-cell or spatial transcriptomics data. _**Nature Biotechnology**_ 2022](https://www.nature.com/articles/s41587-022-01288-0).

# Seminar
**Want to hear more about DIALOGUE and other approaches to study multicellular biology?** Checkout our [seminar at 
BCH Digital Science TV](https://www.youtube.com/watch?v=iBtzD0rKSdM&list=PLZH5lNty_E1pHu2cQY83tDYDCyhJFOd7a&index=2&t=2964s)
