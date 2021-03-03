"""This module defines a Gene class and related functions"""

from collections import defaultdict

class Gene():
    """
    An instance represents one biological gene.

    Instance attributes:
    chrom: int, chromosome on which the gene is located
    positions: tuple of int, numerical position of (start_pos, end_pos)
    snps: list of int, numerical position of all SNPs in the gene
    snp_density: dict, where values are quadrant names and values are SNP density in quadrant
    quad_starts: list of int, start positions for quadrants 1-4 *MUST BE IN ORDER*
    gene_id: str, gene ID
    variation: float, standard deviation of gene expression across all lines
    median: float, median expression of the gene across all lines
    start_pos: int, start position of the gene, calculated from right and left plot and direction
    end_pos: int, end position of the gene, calculated from right and left plot and direction
    """

    def __init__(self, chrom, positions, gene_direction, gene_id, variation, median, snps=[]):
        """
        Initializer, creates a gene instance.

        Parameters:
        chrom: int, chromosome on which the gene is located
        positions: tuple of int, with right and left positions
        gene_direction: string, either '+' or '-', indicates which strand gene is on
        snps: list of int, numerical position of all SNPs in the gene. Can be empty
        gene_id: str, gene ID
        variation: float, standard deviation of gene expression across all lines
        median: float, median of gene expression across all lines
        """
        self.chrom = chrom
        self.gene_direction = gene_direction
        self.gene_id = gene_id
        self.variation = variation
        self.median = median
        self.snps = snps        
        self.start_pos, self.end_pos = self.determine_start_end_pos(positions[0], positions[1], gene_direction)
        self.length = abs(self.start_pos - self.end_pos + 1)
        self.normalized_snps = self.normalize_snps()    # range from 0 to 1
        self.quad_starts = self.calculate_quadrant_positions( [self.start_pos, self.end_pos] , gene_direction )
        self.snp_density = self.calculate_snp_density()
        self.snp_density_per_gene = self.cal_snp_density_per_gene()

    def set_snps(self, snp_list):
        """
        Sets the instance attribute snps to snp_list
        """
        self.snps = snp_list
        
        
    def normalize_snps(self):
        return [ abs( float(x - self.start_pos) / self.length ) for x in self.snps]
    
    
    def calculate_snp_density(self):
        """
        Calculates the snp density in each quadrant.

        Returns: a dictionary where keys are quadrant names and values are snp density in the quadrant
        """
        quad_starts = self.quad_starts
        quad_ends = [*quad_starts[1:4], self.end_pos]
        snpDensity = {"q1":0,"q2":0,"q3":0,"q4":0}

        quad_len = abs(quad_ends[0] - quad_starts[0] + 1)

        for snp in self.snps:
            for i, key in enumerate(snpDensity):
                if quad_starts[i] < snp < quad_ends[i] or quad_starts[i] > snp > quad_ends[i]:
                    snpDensity[key] += 1

        for key in snpDensity:
            snpDensity[key] = float(snpDensity[key]) / quad_len

        return snpDensity


    def cal_snp_density_per_gene(self):
        return float(len(self.snps)) / self.length


    def calculate_quadrant_positions(self, positions, gene_direction):
        """
        Calculates and returns quadrant start positions.

        Parameters:
        start_pos: int, start position of gene
        end_pos: int, end position of gene

        Returns: list of quadrant start positions in order from 1 to 4
        """
        gene_len = abs(positions[0] - positions[1] + 1)
        quad_len = int(gene_len/4)
        q1_start = positions[0]
        if gene_direction == 1:
            q2_start = q1_start + quad_len
            q3_start = q2_start + quad_len
            q4_start = q3_start + quad_len
        else:
            q2_start = q1_start - quad_len
            q3_start = q2_start - quad_len
            q4_start = q3_start - quad_len

        return [q1_start, q2_start, q3_start, q4_start]


    def determine_start_end_pos(self, left_pos, right_pos, gene_direction):
        """
        Determines which of left and right positions are start and end.

        Parameters:
        left_pos: int, position of left most base pair in gene
        right_pos: int, position of right most base pair in gene

        Returns: a tuple of (start_pos, end_pos)
        """
        if gene_direction == '+':
            return (left_pos, right_pos)
        else:
            return (right_pos, left_pos)
