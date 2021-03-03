"""
Script to create a geneData object for a genome and plot correlations.

python generateGenes.py myExpressionData.csv mySNPData.csv
"""
from geneData import geneData
import pandas as pd
import sys
from matplotlib import pyplot as plt
from collections import defaultdict


# read in data
#expressionDataset = pd.read_csv(sys.argv[1])
#SNPDataset = pd.read_csv(sys.argv[2])
path_gene_data = 'data/snp_davis_chr10.csv'
path_snp_data = 'data/snp_imputed_chr10_sample_v2.csv'
expressionDataset = pd.read_csv(path_gene_data,sep='\t')
SNPDataset = pd.read_csv(path_snp_data)

# generate an instance of geneData
print('Starting gene extraction')
maizeGenes = geneData(expressionDataset, SNPDataset)
print('Done with gene extraction, maizeGenes created')

# create plots/regressions
print('Starting plots')
quad_names = ['q1','q2','q3','q4']
x_med = defaultdict(list)
y_med = defaultdict(list)
zero_genes_med = []
for quad in quad_names:
    for gene_id in maizeGenes.dic_gene:
        x_med[quad].append(maizeGenes.dic_gene[gene_id].snp_density[quad])
        y_med[quad].append(maizeGenes.dic_gene[gene_id].median)
        if maizeGenes.dic_gene[gene_id].snp_density[quad] == 0:
            zero_genes_med.append(maizeGenes.dic_gene[gene_id].median)

fig, ax = plt.subplots(2,2, sharex = 'col', sharey = 'row')
plt.xlabel('SNP Density')
plt.ylabel('Median Expression')
ax[0,0].scatter(x_med['q1'],y_med['q1'])
ax[0,0].set_title('Quadrant 1')
ax[0,1].scatter(x_med['q2'],y_med['q2'])
ax[0,1].set_title('Quadrant 2')
ax[1,0].scatter(x_med['q3'],y_med['q3'])
ax[1,0].set_title('Quadrant 3')
ax[1,1].scatter(x_med['q4'],y_med['q4'])
ax[1,1].set_title('Quadrant 4')

plt.savefig('medExp',bbox_inches='tight')
print('Median Expression plots finished')

quad_names = ['q1','q2','q3','q4']
x_var = defaultdict(list)
y_var = defaultdict(list)
x_var_high = defaultdict(list)
y_var_high = defaultdict(list)
zero_genes_var = []
for quad in quad_names:
    for gene_id in maizeGenes.dic_gene:
        x_var[quad].append(maizeGenes.dic_gene[gene_id].snp_density[quad])
        y_var[quad].append(maizeGenes.dic_gene[gene_id].variation)
        if maizeGenes.dic_gene[gene_id].snp_density[quad] == 0:
            zero_genes_var.append(maizeGenes.dic_gene[gene_id].variation)


fig, ax = plt.subplots(2,2, sharex = 'col', sharey = 'row')
plt.xlabel('SNP Density')
plt.ylabel('Variation in Expression')
plt.title('Variation in Gene Expression vs. SNP Density')
ax[0,0].scatter(x_var['q1'],y_var['q1'])
ax[0,0].set_title('Quadrant 1')
ax[0,1].scatter(x_var['q2'],y_var['q2'])
ax[0,1].set_title('Quadrant 2')
ax[1,0].scatter(x_var['q3'],y_var['q3'])
ax[1,0].set_title('Quadrant 3')
ax[1,1].scatter(x_var['q4'],y_var['q4'])
ax[1,1].set_title('Quadrant 4')

plt.savefig('varExp',bbox_inches='tight')

print('Variation in Expression plot finished')


x_var = defaultdict(list)
y_var = defaultdict(list)
for quad in quad_names:
    for gene_id in maizeGenes.dic_gene:
        if maizeGenes.dic_gene[gene_id].snp_density[quad] > 0:
            x_var[quad].append(maizeGenes.dic_gene[gene_id].snp_density[quad])
            y_var[quad].append(maizeGenes.dic_gene[gene_id].variation)

for quad in quad_names:
    plt.figure(quad + 'snp_var')
    plt.scatter(x_var[quad], y_var[quad], s=5)
    plt.yscale('log')
    plt.xlabel('snp density')
    plt.ylabel('variance')
    plt.title('snp density and variance for %s' % quad) 
    plt.xlim((0, 0.15))
    plt.savefig(quad + '_snp_density_variance.png')



'''
plot relation between snp density of full gene and expression variation
Genes without snps are discarded in plotting.
It seems that the minimal variation increases with snp density.
'''
snp_density_per_gene_list = \
[maizeGenes.dic_gene[gene_id].snp_density_per_gene for gene_id in maizeGenes.dic_gene if maizeGenes.dic_gene[gene_id].snp_density_per_gene > 0]
variation_list = \
[maizeGenes.dic_gene[gene_id].variation for gene_id in maizeGenes.dic_gene if maizeGenes.dic_gene[gene_id].snp_density_per_gene > 0]
plt.figure()
plt.scatter(snp_density_per_gene_list, variation_list, c='r', s=4)
plt.xlabel('snp density')
plt.ylabel('expression variation')
plt.title('Snp density and expression variation (genes without snps are ignored)')
plt.yscale('log')
plt.savefig('snp2expression.png')

print('Expression variation plots for non-zero SNP genes finished')

'histogram of snp density for each quadrant'
num_bins = 20
for quad in quad_names:
    plt.figure(quad + 'snp')
    plt.hist(x_var[quad], bins=num_bins)
    plt.yscale('log')
    plt.xlabel('snp density')
    plt.ylabel('number of genes')
    plt.title('histogram of snp density for %s' % quad)
    plt.savefig(quad + '_snp_density_histogram.png')

print('Histogram of snp density for each quadrant finished')

'histogram of snp density for full gene'
snp_density_per_gene_list = \
[maizeGenes.dic_gene[gene_id].snp_density_per_gene for gene_id in maizeGenes.dic_gene]

plt.figure('full_snp')
plt.hist(snp_density_per_gene_list, bins=num_bins)
plt.yscale('log')
plt.xlabel('snp density per gene')
plt.ylabel('number of genes')
plt.title('histogram of snp density for full gene')
plt.savefig('full_snp_density_histogram.png')

print('Histogram of snp density for full gene finished')

'''
histogram of snp numbers at different gene locations (relative locations)\
implementation of Davis's plot by using the gene classes
'''
snp_location_list = []
for gene_id in maizeGenes.dic_gene:
    snp_location_list += maizeGenes.dic_gene[gene_id].normalized_snps

num_bins = 50
plt.figure('snp_location')
plt.hist(snp_location_list, bins=num_bins)
plt.xlabel('relative location in a gene')
plt.ylabel('number of snps')
plt.title('histogram of snp locations')
plt.savefig('snp_location_histogram.png')

print('Davis\' plot finished')

'''
scatterplot of mean vs variation for genes with 0 SNPs
'''
plt.figure('Genes with 0 SNPs')
plt.scatter(zero_genes_med,zero_genes_var)
plt.xlabel('Median Expression (FPKM)')
plt.ylabel('Expression Variation (Standard Deviation)')
plt.title('Median vs. Variation of Expression for Genes with 0 SNPs')
plt.savefig('med_var_exp.png',bbox_inches='tight')

print('Median vs Expression scatterplot finished')
# calculate significance of correlations
