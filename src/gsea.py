"""
disclaimer:
I used cursor to help me with some syntax (such as how to melt dataframes etc.)
the same way I would have used stackoverflow in the past.
However, the choice of data structures, implementation and logic is entirely my own.
"""

import sys
import pandas as pd
import numpy as np


def read_geneset_definitions(gmt_file):
    """
    read KEGG gene sets
    input:
        gmt_file (str): .gmt file
    return:
        dict: dict from pathway name to list of genes {pathway_name: [gene1, gene2, ...]}
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
    read expression data file
    input:
        human_gene_expression_file (str): expression file
    return:
        pandas DataFrame: df with genes as index, samples as columns
    """
    human_gene_expression_df = pd.read_csv(human_gene_expression_file, sep='\t', index_col=0)
    return human_gene_expression_df


def read_human_condition_as_df(samples_txt_file):
    """
    read sample labels and return as dictionary.
    input:
        samples_txt_file (str):sample file
    return:
        dict: dict mapping sample_id to condition
    """
    labels_df = pd.read_csv(samples_txt_file, sep='\t', header=None,
                            names=['sample_id', 'condition'])
    return labels_df


def read_human_condition_as_dict(samples_txt_file):
    """
    read sample labels file and return as dict
    input: samples_txt_file (str): Path to sample file
    return: dict: from sample_id to condition
    """
    labels_df = read_human_condition_as_df(samples_txt_file)
    labels_dict = dict(zip(labels_df['sample_id'], labels_df['condition']))
    return labels_dict


class GSEA:
    """
    class to implement gsea algo
    """
    
    def __init__(self):
        """
        initialize class.
        """
        self.human_gene_expression_df = None
        self.human_condition_dict = None
        self.human_condition_df = None
        self.geneset_definitions_dict = None
        self.ranked_genes_list = None

    def load_data(self, expfile, sampfile, genesets):
        """
        load expression, sample and gene sets
        input:
            expfile (str): expression file
            sampfile (str): sample file
            genesets (str): gene sets file
        return:
            None
        """
        self.geneset_definitions_dict = read_geneset_definitions(genesets)
        self.human_condition_dict = read_human_condition_as_dict(sampfile)
        self.human_condition_df = read_human_condition_as_df(sampfile)
        self.human_gene_expression_df = read_human_gene_expression_file(expfile)

    def get_gene_rank_order(self):
        """
        rank genes by log fold-change
        return: list: gen names (strings) by logFC in desc order
        """
        # EFFICIENCY FIX: pivot the expressions dataframe to turn it into dataframe
        # pivot_human_gene_expression_df has columns: SYMBOL, sample_id, gene_expression
        pivot_human_gene_expression_df = self.human_gene_expression_df.reset_index().melt(
            id_vars="SYMBOL",
            var_name="sample_id",
            value_name="gene_expression"
        )
        # EFFICIENCY FIX: join it with sample data on sampe_id to now include condition column to df
        # pivot_human_gene_expression_join_condition_df has columns: SYMBOL, sample_id, gene_expression, condition
        pivot_human_gene_expression_join_condition_df = pd.merge(pivot_human_gene_expression_df,
                                                                 self.human_condition_df,
                                                                 on="sample_id")
        # EFFICIENCY FIX: get mean expression per gene, per condition
        # mean_expression_per_gene_condition_df has columns: SYMBOL, condition, gene_expression (mean)
        mean_expression_per_gene_condition_df = pivot_human_gene_expression_join_condition_df[["SYMBOL", "condition", "gene_expression"]].groupby(["SYMBOL", "condition"]).mean().reset_index()

        # EFFICIENCY FIX: Pivot again to make conditions (0, 1) become columns
        # pivoted_mean_expression_per_gene_condition_df has columns: SYMBOL, condition, 0, 1 and values are mean gene_expression
        pivoted_mean_expression_per_gene_condition_df = mean_expression_per_gene_condition_df.pivot(
            index="SYMBOL",
            columns="condition",
            values="gene_expression"
        )
        # EFFICIENCY FIX: subtract column 0 from column 1 to get log_fc_gene per gene/ SYMBOL
        pivoted_mean_expression_per_gene_condition_df["log_fc_gene"] = pivoted_mean_expression_per_gene_condition_df[1] - pivoted_mean_expression_per_gene_condition_df[0]
        # Rank it by expression in descending order
        ranked_gene_df = pivoted_mean_expression_per_gene_condition_df.sort_values(by="log_fc_gene", ascending=False)
        self.ranked_genes_list = ranked_gene_df.index.tolist()
        return self.ranked_genes_list

    def get_enrichment_score(self, geneset):
        """
        calculate enrichment score for a given gene set
        input:
            geneset (str):gene set name
        return:
            float: enrichment score
        """

        # get ranked genes
        if self.ranked_genes_list is None:
            self.ranked_genes_list = self.get_gene_rank_order()

        #  EFFICIENCY FIX: use sets instead of lists for quick verification
        ranked_genes_set = set(self.ranked_genes_list)
        set_of_genes_in_geneset = self.geneset_definitions_dict[geneset]
        set_of_genes_in_geneset = set([gene for gene in set_of_genes_in_geneset if gene in ranked_genes_set])

        # num genes in L and S (N and G resp)
        N_num_genes_universal_L = len(self.ranked_genes_list)
        G_num_genes_in_geneset_S = len(set_of_genes_in_geneset)

        # divide by 0 avoid
        if G_num_genes_in_geneset_S == 0 or G_num_genes_in_geneset_S == N_num_genes_universal_L:
            return 0.0

        # up and down scores per lecture
        up_score = np.sqrt((N_num_genes_universal_L-G_num_genes_in_geneset_S)/G_num_genes_in_geneset_S)
        down_score = np.sqrt(G_num_genes_in_geneset_S/(N_num_genes_universal_L-G_num_genes_in_geneset_S))

        # iterate over ranked genes to compute max score
        score = 0.0
        geneset_enrichment_score = None
        for gene in self.ranked_genes_list:
            if gene in set_of_genes_in_geneset:
                score = score + up_score
            else:
                score = score - down_score

            if geneset_enrichment_score is None:
                geneset_enrichment_score = score

            if score > geneset_enrichment_score:
                geneset_enrichment_score = score

        geneset_enrichment_score = np.round(geneset_enrichment_score, 2)

        return geneset_enrichment_score

    def get_sig_sets(self, p):
        """
        get enriched gene sets using permutation test
        input:
            p (float): threshold
        return:
            list: significant genes list
        """
        # calculate actual enrichment scores for all gene sets
        original_condition_df = self.human_condition_df.copy()
        actual_lebel_geneset_enrichment_score_dict = dict()
        for geneset_name in self.geneset_definitions_dict:
            actual_lebel_geneset_enrichment_score_dict[geneset_name] = self.get_enrichment_score(geneset_name)

        # permute conditions 100 times for each geneset to get enrichment scores
        permutated_label_geneset_enrichment_scores_dict = dict()
        for i in range(100):

            # permute condition
            self.human_condition_df['condition'] = np.random.permutation(
                original_condition_df['condition'].values
            )

            # set to none to start afresh
            self.ranked_genes_list = None

            for geneset_name in self.geneset_definitions_dict:
                # initialize geneset in permutated_label_geneset_enrichment_scores_dict
                if geneset_name not in permutated_label_geneset_enrichment_scores_dict:
                    permutated_label_geneset_enrichment_scores_dict[geneset_name] = list()

                permutated_label_geneset_enrichment_scores_dict[geneset_name].append(self.get_enrichment_score(geneset_name))

        # restore original after permutations
        self.human_condition_df = original_condition_df
        self.ranked_genes_list = None

        # get p-values
        # for each gene set count how many permuted score > actual score
        num_genesets = len(self.geneset_definitions_dict)

        significant_sets = list()
        for geneset_name in self.geneset_definitions_dict:
            ES_enrichment_score_actual = actual_lebel_geneset_enrichment_score_dict[geneset_name]
            num_perm_es_greater_actual_es = len([perm_score for perm_score in permutated_label_geneset_enrichment_scores_dict[geneset_name] if perm_score >= ES_enrichment_score_actual])

            p_value = num_perm_es_greater_actual_es/100
        
            # Bonferroni
            corrected_p = p_value * num_genesets

            if corrected_p < p:
                significant_sets.append(geneset_name)

        return significant_sets


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 gsea.py expfile sampfile keggfile")
        sys.exit(1)
    
    human_gene_expression_file = sys.argv[1]
    human_condition_file = sys.argv[2]
    geneset_definitions_file = sys.argv[3]
    
    gsea = GSEA()
    gsea.load_data(human_gene_expression_file, human_condition_file, geneset_definitions_file)
    
    example_pathway = "KEGG_CITRATE_CYCLE_TCA_CYCLE"
    gene_enrichment_score = gsea.get_enrichment_score(example_pathway)
    print(f"Enrichment score for {example_pathway}: {gene_enrichment_score}")
    
    sig_sets = gsea.get_sig_sets(0.05)
    print(f"Found {len(sig_sets)} significant gene sets: {sig_sets}")


if __name__ == "__main__":
    main()

