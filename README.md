# CIE
## Causal Inference Using Directional and Unidirectional Enrichment Analysis on Biological Networks

Identification of active regulatory pathways under specific molecular and environmental perturbations is crucial for understanding cellular function. However, inferring causal regulatory mechanisms directly from differential gene expression data is very difficult and remains a central challenge in computational biology.  We present CIE, an early version of a user-friendly web-app, and a tool box which facilitate inference of casual regulatory mechanisms based on existing knowledge bases of regulator-target interactions and user-provided differential gene expression data.  

## Downloading the R Package

To download our R package, please install the package devtools if it is not already installed on your system.  
`install.packages("devtools")`

Using devtools, you can install our package directly from GitHub  
`devtools::install_github("umbibio/CIE")`

## Data for CIE
Pre-processed databases depicting interactions in human cells are available [here](http://markov.math.umb.edu/CIE).  Click the "Download files" button at the top of the page.  This will take you to a page where you can download a compressed bundle of the selected files.

## Additional information
For more information, and guidance on how to get information from databases so that CIE can analyze differential gene expression data sourced from species other than human, please visit our [Wiki](https://github.com/umbibio/CIE/wiki)
