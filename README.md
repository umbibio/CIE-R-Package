# CIE
## Causal Inference Using Directional and Unidirectional Enrichment Analysis on Biological Networks

Identification of active regulatory pathways under specific molecular and environmental perturbations is crucial for understanding cellular function. However, inferring causal regulatory mechanisms directly from differential gene expression data is very difficult and remains a central challenge in computational biology.  We present CIE, an early version of a user-friendly web-app, and a tool box which facilitate inference of casual regulatory mechanisms based on existing knowledge bases of regulator-target interactions and user-provided differential gene expression data.  

## Installing the R Package

Installation of this package can be accomplished through devtools or remotes.  However there is one package dependency, org.Hs.eg.db, which cannot be automatically fetched.  Before installing, run the following lines of code to install this package.  To check if the package has already been installed run

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
