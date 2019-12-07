# CIE
## Causal Inference Engine: Directional Enrichment Analysis on Biological Networks

CIE is a platform for inference of active transcriptional regulators of differential gene expression. CIE provides a user-friendly web-app and an R-package to run inference queries on casual regulatory netowrks. Several regulatory networks are provided in the platform, including a TF-gene intraction network constructed from ChIP-seq and tissue specific gene expression data, a netowrk based on STRING-DB, and curated TF-gene interaction networks such as TRRUST and TRED. The R-package provides functionality for custom networks as well. 

## Citations 
1) Farahman, S., O'Connor, C., Macoska, J., Zarringhalam, K. 'Causal Inference En\gine: A platform for directional gene set enr\ichment analysis and inference of active tran\scriptional regulators', NAR doi:  https://doi.o\rg/10.1093/nar/gkz1046

2) Fakhry CT, Choudhary P, Gutteridge A, Sidders B, Chen P, Ziemek D, Zarringhalam, K. 'Interpreting transcriptional changes using causal graphs: new methods and their practical utility on public networks'. BMC Bioinformatics. 2016;17(1):318

## Installing the R Package

Installation of this package can be accomplished through devtools or remotes. One package dependency, org.Hs.eg.db, must be manually fetched.  Before installing, run the following lines of code to install this package.  To check if the package has already been installed run

library(org.Hs.eg.db)

### Installing org.Hs.eg.db

#### R >= 3.5
If you do not already have BiocManager installed, do so with the following command
```R
install.packages("BiocManager")
```
Then install org.Hs.eg.db with
```R
BiocManager::install("org.Hs.eg.db")
```

#### R < 3.5
```R
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("org.Hs.eg.db")
```

### Downloading and installing CIE
#### With remotes
If you do not already have the package remotes on your system, install it with
```R
install.packages("remotes")
```
CIE can then be installed with
```R
remotes::install_github("umbibio/CIE-R-Package")
```
#### With devtools
If you do not already have the package devtools on your system, install it with
```R
install.packages("devtools")
```
CIE can then be installed with
```R
devtools::install_github("umbibio/CIE-R-Package")
```

## Enabling Graphing
Graphing is not enabled by default, nor is the library for it automatically loaded by the package.  If you wish to access these features, download the stable version of rcytoscapejs here.
[ryctoscapejs](https://github.com/cytoscape/cyjShiny/releases/tag/v0.0.7)
Once the source code is decompressed, you can do a local install with
```R
install.packages("path/to/rcystoscapejs/", repos=NULL, type="source")
```
or
```R
devtools::install("path/to/rcytoscapejs/")
```
Before making calls to our graphing function you must manually load the library as shown below
```R
library(rcytoscapejs)
```

## Data for CIE
Pre-processed databases depicting interactions in human cells are available [here](https://markov.math.umb.edu/app/cie).  Click the "Download files" button at the top of the page.  This will take you to a page where you can download a compressed bundle of the selected files.

## Additional information
For more information, and guidance on how to get information from databases so that CIE can analyze differential gene expression data sourced from species other than human, please visit our [Wiki](https://github.com/umbibio/CIE/wiki)
