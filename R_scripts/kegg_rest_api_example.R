# kegg_rest_api_example.R
# Examples using the KEGGREST bioconductor package
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 08 Mar 2017
# Updated on 08 Mar 2017

require(KEGGREST)
# https://bioconductor.riken.jp/packages/3.0/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.R
listDatabases()
head(keggConv("eco", "ncbi-geneid"))
keggConv("ncbi-geneid", "hsa:10458")
keggConv("genes", "ncbi-geneid:1246517")
