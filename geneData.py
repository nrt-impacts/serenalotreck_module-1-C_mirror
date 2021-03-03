"""Module defining the geneData class"""

from Genes import Gene
import pandas as pd

class geneData():
    """
    An instance represents all genes from a genome.

    Instance attributes:
    num_gene: int, number of genes in set
    dic_gene: dictionary, keys are gene_id's and values are Gene instances
    """
    def __init__(self, df_gene, df_snp):
        """
        Initializer for geneData.

        Extracts information from data to create instances of the Gene class
        and stores them in a dictionary.

        Paramters:
        df_gene: csv, expression dataset
        df_snp: csv
        """
        self.num_gene = len(df_gene)
        self.dic_gene = {}

        gene2snp_dic = self.cal_snp_per_gene(df_gene, df_snp)
        for i in range(self.num_gene):
            chrom = df_gene['chromosome'][i]
            positions = [ df_gene['position_left'][i], df_gene['position_right'][i] ]
            gene_direction = df_gene['Direction'][i]
            gene_id = df_gene['gene'][i]
            variation = df_gene['Variance'][i]
            median = df_gene['Median'][i]
            snps = gene2snp_dic[gene_id]
            self.dic_gene[gene_id] = Gene(chrom, positions, gene_direction, gene_id, variation, median, snps)


    def cal_snp_per_gene(self, df_gene, df_snp):
        '''
        Compute SNP positions in each gene
        input:
            df_gene - gene DataFrame
            df_snp - SNP DataFrame
        ouput:
            dic_gene_snp: a dictionary from gene to all its SNPs; key: gene name, value: a list of its SNP positions
        '''

        n_gene, n_snp = len(df_gene), len(df_snp)

        dic_gene_snp = { df_gene['gene'][i] : [] for i in range(n_gene) }

        for i in range(n_snp):
            chrom, pos = df_snp['chromosome'][i], df_snp['position'][i]
            mask = (pos > df_gene['position_left']) & (pos < df_gene['position_right']) \
                    & (df_gene['chromosome']==chrom)

            if mask.sum()==1:
                gene_name = df_gene['gene'][mask].iloc[0]
                dic_gene_snp[ gene_name ].append(pos)

        return dic_gene_snp


if __name__=='__main__':
    path_gene_data = 'data\\condensed_subset_Davis.csv'
#    path_snp_data = 'data\\condensed_SNP_subset_Davis.csv'
    path_snp_data = 'data\\snp_imputed_chr10_sample_v2.csv'

    df_gene = pd.read_csv(path_gene_data)
    df_snp = pd.read_csv(path_snp_data)
    #df_snp = pd.read_csv(path_snp_data, sep='\t')   # for condensed_SNP_subset_Davis.csv

    maizeGenes = geneData(df_gene, df_snp)

    # show info of gene Zm00001d024939
    gene1 = maizeGenes.dic_gene['Zm00001d024939']
    print(gene1.gene_id)
    print('chromosome:', gene1.chrom)
    print('std:', gene1.variation)
    print('snps:', gene1.snps)
    print('start_pos:', gene1.start_pos)
    print('end_pos:', gene1.end_pos)
    print('quad_starts:', gene1.quad_starts)
    print('snp_density:', gene1.snp_density )
