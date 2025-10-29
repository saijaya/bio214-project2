"""
BIO214 Project 2 - Part 2: Gene Set Enrichment Analysis (GSEA)
Author: [Your Name]
Date: [Date]
Collaborators: [List collaborators or state "None"]

This script implements GSEA algorithm to identify enriched KEGG pathways
in endometriosis patients vs. healthy controls.

Usage: python3 gsea.py expfile sampfile keggfile
"""

import sys
import pandas as pd
import numpy as np


class GSEA:
    """
    Class to perform Gene Set Enrichment Analysis.
    """
    
    def __init__(self):
        """
        Initialize the GSEA class.
        Store any instance variables you need here.
        """
        self.expression_data = None
        self.sample_data = None
        self.gene_sets = None
        self.ranked_genes = None
        # Add any other instance variables you need
    
    
    def load_data(self, expfile, sampfile, genesets):
        """
        Load expression data, sample assignments, and KEGG gene sets.
        
        Args:
            expfile (str): Path to expression data file
            sampfile (str): Path to sample assignment file (sample_id, condition)
            genesets (str): Path to KEGG gene sets file (.gmt format)
        
        Returns:
            None (stores data in instance variables)
        
        Hint: Use pandas to read the files efficiently
              For .gmt file: each row is a gene set, tab-delimited
              Column 1 = pathway name, Column 2 = URL (ignore), rest = genes
        """
        # TODO: Read expression data
        # Likely format: rows = genes, columns = samples
        
        # TODO: Read sample assignments
        # Format: sample_id, condition (healthy or endometriosis)
        
        # TODO: Read KEGG gene sets from .gmt file
        # Store as dictionary: {pathway_name: [list of genes]}
        
        pass
    
    
    def get_gene_rank_order(self):
        """
        Rank genes by log fold-change between endometriosis and healthy samples.
        
        Log fold-change = mean_log_expression(endometriosis) - mean_log_expression(healthy)
        
        Returns:
            list: Gene names (strings) sorted by logFC, highest first
        
        Hint: Since expression is already log-normalized, just subtract means
              Use pandas groupby and vectorized operations for efficiency
        """
        # TODO: Separate samples into healthy and endometriosis groups
        
        # TODO: Calculate mean expression for each group per gene
        
        # TODO: Calculate log fold-change (endometriosis - healthy)
        
        # TODO: Sort genes by logFC (descending order)
        
        ranked_genes = []
        
        return ranked_genes
    
    
    def get_enrichment_score(self, geneset):
        """
        Calculate enrichment score for a given gene set using Brownian bridge method.
        
        Steps:
        1. Get ranked gene list (call get_gene_rank_order if needed)
        2. Identify which ranked genes are in the gene set
        3. Calculate up_score and down_score based on gene set size
        4. Walk through ranked list, adding up_score if gene in set, 
           subtracting down_score if not
        5. Return supremum (maximum) of the running sum as ES
        
        Args:
            geneset (str): Name of gene set (e.g., 'KEGG_CITRATE_CYCLE_TCA_CYCLE')
        
        Returns:
            float: Enrichment score, rounded to 2 decimal places
        
        Note: Only consider genes that are in both the gene set AND expression data
              We're looking for UP-regulated gene sets, so don't take absolute value
        """
        # TODO: Get or use ranked gene list
        
        # TODO: Get genes in this gene set (only those in expression data)
        
        # TODO: Calculate appropriate scores
        # up_score = depends on # genes in set
        # down_score = depends on # genes NOT in set
        
        # TODO: Walk through ranked list, calculate running sum
        # Add up_score if gene in set, subtract down_score if not
        
        # TODO: Find supremum (max value) of running sum
        
        # TODO: Round to 2 decimal places
        
        enrichment_score = 0.0
        
        return round(enrichment_score, 2)
    
    
    def get_sig_sets(self, p):
        """
        Identify significantly enriched gene sets using permutation testing.
        
        Steps:
        1. Calculate actual enrichment scores for all gene sets
        2. Permute sample labels 100 times:
           - For each permutation, shuffle sample-to-condition assignments
           - Recalculate ranked gene list
           - Recalculate enrichment scores for all gene sets
        3. Calculate p-value for each gene set:
           p-value = (# times permuted ES >= actual ES) / 100
        4. Apply Bonferroni correction: corrected_p = p-value * number_of_gene_sets
        5. Return gene sets where corrected_p < threshold p
        
        Args:
            p (float): Significance threshold for corrected p-values
        
        Returns:
            list: Names of significant gene sets (strings)
                  Empty list if no gene sets are significant
        
        IMPORTANT: Must run efficiently! Should complete in < 3 minutes
        Use vectorized operations, avoid loops where possible
        """
        # TODO: Calculate actual enrichment scores for all gene sets
        
        # TODO: Permutation testing (100 iterations)
        # For each iteration:
        #   - Shuffle sample labels (np.random.permutation can help)
        #   - Recalculate gene rankings
        #   - Recalculate ES for all gene sets
        #   - Store results
        
        # TODO: Calculate p-values
        # For each gene set: count how many permuted ES >= actual ES
        # p-value = count / 100
        
        # TODO: Apply Bonferroni correction
        # corrected_p = p-value * number_of_gene_sets
        
        # TODO: Filter for significant gene sets (corrected_p < p)
        
        significant_sets = []
        
        return significant_sets


def main():
    """
    Main function to run GSEA from command line.
    This is for your own testing - autograder calls class methods directly.
    """
    if len(sys.argv) != 4:
        print("Usage: python3 gsea.py expfile sampfile keggfile")
        sys.exit(1)
    
    exp_file = sys.argv[1]
    samp_file = sys.argv[2]
    kegg_file = sys.argv[3]
    
    # Create GSEA object and run analysis
    gsea = GSEA()
    gsea.load_data(exp_file, samp_file, kegg_file)
    
    # Get ranked genes
    ranked_genes = gsea.get_gene_rank_order()
    print(f"Ranked {len(ranked_genes)} genes")
    
    # Example: Get enrichment score for a specific pathway
    example_pathway = "KEGG_CITRATE_CYCLE_TCA_CYCLE"
    # es = gsea.get_enrichment_score(example_pathway)
    # print(f"Enrichment score for {example_pathway}: {es}")
    
    # Get significant gene sets at p < 0.05
    sig_sets = gsea.get_sig_sets(0.05)
    print(f"Found {len(sig_sets)} significant gene sets")
    
    # Optional: Output results for your own analysis


if __name__ == "__main__":
    main()

