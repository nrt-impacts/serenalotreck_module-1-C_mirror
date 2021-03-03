# module-1-C
IMPACTS Frontiers Module 1, Group C: Yunfei Long, Davis Mathieu, Serena Lotreck, and Wei Dong

# Overview
In this project, we are seeking to determine if SNP correlation in certain regions of the gene is correlated to variation in gene expression across inbred maize lines. In order to do this, we will test the significance of the correlation between SNP density in a given quadrant and expression variation. 

# Code Structure 
**Genes.py** defines the Gene class. Each instance of this class represents a single gene. <br>
**geneData.py** defines the geneData class. An instance of this class contains the identifiers of all instances of Genes for a given dataset. <br>
**data_extraction.py** defines a function to get the SNPs that fall within each gene in the genome. <br>
**generateGenes.py** is the script in which the geneData instance for the dataset is created, and in which all calculations and plotting is performed. <br>
