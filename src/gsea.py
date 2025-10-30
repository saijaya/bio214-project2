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


def read_geneset_definitions(gmt_file):
    """
    Read KEGG gene sets from .gmt file format.

    Args:
        gmt_file (str): Path to .gmt file

    Returns:
        dict: Dictionary mapping pathway name to list of genes
              {pathway_name: [gene1, gene2, ...]}
    """
    gene_sets = {}

    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pathway_name = parts[0]
                genes = parts[2:]
                gene_sets[pathway_name] = set(genes)

    return gene_sets


def read_human_gene_expression_file(human_gene_expression_file):
    """
    Read expression data file.

    Args:
        human_gene_expression_file (str): Path to expression file

    Returns:
        pd.DataFrame: Expression data with genes as index, samples as columns
    """
    # Read tab-separated file, first column is gene names (index)
    human_gene_expression_df = pd.read_csv(human_gene_expression_file, sep='\t', index_col=0)
    return human_gene_expression_df


def read_human_condition_as_df(samples_txt_file):
    """
    Read sample labels file and return as dictionary.

    Args:
        samples_txt_file (str): Path to sample file

    Returns:
        dict: Dictionary mapping sample_id to condition {sample_id: 0 or 1}
    """
    labels_df = pd.read_csv(samples_txt_file, sep='\t', header=None,
                            names=['sample_id', 'condition'])
    print(labels_df)
    return labels_df


def read_human_condition_as_dict(samples_txt_file):
    """
    Read sample labels file and return as dictionary.

    Args:
        samples_txt_file (str): Path to sample file

    Returns:
        dict: Dictionary mapping sample_id to condition {sample_id: 0 or 1}
    """
    labels_df = read_human_condition_as_df(samples_txt_file)
    labels_dict = dict(zip(labels_df['sample_id'], labels_df['condition']))
    return labels_dict

class GSEA:
    """
    Class to perform Gene Set Enrichment Analysis.
    """
    
    def __init__(self):
        """
        Initialize the GSEA class.
        Store any instance variables you need here.
        """
        self.human_gene_expression_df = None
        self.human_condition_dict = None
        self.human_condition_df = None
        self.geneset_definitions_dict = None
        self.ranked_genes_list = None
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
        self.geneset_definitions_dict = read_geneset_definitions(genesets)
        print(len(self.geneset_definitions_dict))
        
        self.human_condition_dict = read_human_condition_as_dict(sampfile)
        print(self.human_condition_dict)

        self.human_condition_df = read_human_condition_as_df(sampfile)
        print(self.human_condition_df)

        self.human_gene_expression_df = read_human_gene_expression_file(expfile)
        print(self.human_gene_expression_df)

    def get_gene_rank_order(self):
        """
        Rank genes by log fold-change between endometriosis and healthy samples.
        
        Log fold-change = mean_log_expression(endometriosis) - mean_log_expression(healthy)
        
        Returns:
            list: Gene names (strings) sorted by logFC, highest first
        
        Hint: Since expression is already log-normalized, just subtract means
              Use pandas groupby and vectorized operations for efficiency
        """
        print(self.human_gene_expression_df)

        pivot_human_gene_expression_df = self.human_gene_expression_df.reset_index().melt(
            id_vars="SYMBOL",
            var_name="sample_id",
            value_name="gene_expression"
        )
        print(pivot_human_gene_expression_df)

        pivot_human_gene_expression_join_condition_df = pd.merge(pivot_human_gene_expression_df,
                                                                 self.human_condition_df,
                                                                 on="sample_id")
        print(pivot_human_gene_expression_join_condition_df)

        mean_expression_per_gene_condition_df = pivot_human_gene_expression_join_condition_df[["SYMBOL", "condition", "gene_expression"]].groupby(["SYMBOL", "condition"]).mean().reset_index()
        print(mean_expression_per_gene_condition_df)

        # Pivot to make conditions (0, 1) become columns
        pivoted_mean_expression_per_gene_condition_df = mean_expression_per_gene_condition_df.pivot(
            index='SYMBOL',
            columns='condition',
            values='gene_expression'
        )

        pivoted_mean_expression_per_gene_condition_df["log_fc_gene"] = pivoted_mean_expression_per_gene_condition_df[1] - pivoted_mean_expression_per_gene_condition_df[0]

        print(pivoted_mean_expression_per_gene_condition_df)

        ranked_gene_df = pivoted_mean_expression_per_gene_condition_df.sort_values(by="log_fc_gene", ascending=False)
        print(ranked_gene_df)
        self.ranked_genes_list = ranked_gene_df.index.tolist()
        print(self.ranked_genes_list[:10])
        print(self.ranked_genes_list[-10:])

        return self.ranked_genes_list

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
        ranked_genes_list = self.ranked_genes_list
        ranked_genes_set = set(ranked_genes_list)
        set_of_genes_in_geneset = self.geneset_definitions_dict[geneset]
        set_of_genes_in_geneset = set([gene for gene in set_of_genes_in_geneset if gene in ranked_genes_set])

        N_num_genes_universal_L = len(ranked_genes_list)
        G_num_genes_in_geneset_S = len(set_of_genes_in_geneset)
        print(f"num genes universal: {N_num_genes_universal_L}")
        print(f"num genes in geneset \"{geneset}\": {G_num_genes_in_geneset_S}")

        up_score = np.sqrt((N_num_genes_universal_L-G_num_genes_in_geneset_S)/G_num_genes_in_geneset_S)
        down_score = np.sqrt(G_num_genes_in_geneset_S/(N_num_genes_universal_L-G_num_genes_in_geneset_S))

        score = 0.0
        score_list = list()

        for gene in ranked_genes_list:
            if gene in set_of_genes_in_geneset:
                score = score + up_score
            else:
                score = score - down_score

            score_list.append(score)

        geneset_enrichment_score = np.round(max(score_list), 2)
        
        return geneset_enrichment_score
    
    
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
    
    human_gene_expression_file = sys.argv[1]
    human_condition_file = sys.argv[2]
    geneset_definitions_file = sys.argv[3]
    
    # Create GSEA object and run analysis
    gsea = GSEA()
    gsea.load_data(human_gene_expression_file, human_condition_file, geneset_definitions_file)
    
    # Get ranked genes
    ranked_genes = gsea.get_gene_rank_order()
    print(f"Ranked {len(ranked_genes)} genes")
    
    # Example: Get enrichment score for a specific pathway
    example_pathway = "KEGG_CITRATE_CYCLE_TCA_CYCLE"
    gene_enrichment_score = gsea.get_enrichment_score(example_pathway)
    print(f"Enrichment score for {example_pathway}: {gene_enrichment_score}")
    
    # Get significant gene sets at p < 0.05
    # sig_sets = gsea.get_sig_sets(0.05)
    # print(f"Found {len(sig_sets)} significant gene sets")
    
    # Optional: Output results for your own analysis


if __name__ == "__main__":
    main()

