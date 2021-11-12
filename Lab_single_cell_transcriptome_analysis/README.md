# Lab_single_cell_transcriptome_analysis

## Script for single cell transcriptome analysis by using R Seurat package. We used Seurat v3.2 package (https://satijalab.org/seurat/) in the R environment to conduct a follow-up study on the RNA count matrix data of SC samples.

## 1. System requirements
- Hardware recommendation: CPU thread >=4, RAM >=16GB
- Common operating systems: Windows, Linux or Mac OS
- R software, version 3.6 or greater.
- Seurat package version 3.2

## 2. Installation guide
- 1.Install R, version 3.6 or greater. (https://www.r-project.org/). Install R Studio (Recommended) (https://rstudio.com/)
- 2.Installation Instructions for Seurat (v3.2 or greater): remotes::install_github("satijalab/seurat", ref = "release/4.0.0")

## 3.Script
- This script contains how to merge many single cell datas and normalise the whole data by Seurat. In addition, it contains many other analysis about single cell, such as how to cluster the data by map, how to find DEgenes of all clusters, how to plot makergenes expression and so on.

## 4.Single cell datas availability
- Chicken granulosa cells(GSE181756) and chicken hearts(GSE149457)
