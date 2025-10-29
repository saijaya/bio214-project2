"""
Helper Functions for BIO214 Project 2

This file contains utility functions that may be useful for implementing
PPI and GSEA algorithms. Use these as INSPIRATION to write your own code.

DO NOT submit this file - it's just for your reference!
"""

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.metrics.pairwise import cosine_distances


# ============================================================================
# HELPER FUNCTIONS FOR PPI (Part 1)
# ============================================================================

def read_gene_list(filepath):
    """
    Read a list of genes from a file (one gene per line).
    
    Args:
        filepath (str): Path to gene list file
    
    Returns:
        list: List of gene names (strings)
    """
    genes = []
    with open(filepath, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene:  # Skip empty lines
                genes.append(gene)
    return genes


def load_interaction_network(filepath):
    """
    Load interaction network into a networkx graph.
    
    Args:
        filepath (str): Path to tab-delimited network file
                       Format: gene1 \t gene2 \t weight
    
    Returns:
        networkx.Graph: Undirected graph of interactions
    """
    # Option 1: Using pandas
    df = pd.read_csv(filepath, sep='\t', header=None, names=['gene1', 'gene2', 'weight'])
    G = nx.from_pandas_edgelist(df, source='gene1', target='gene2')
    
    # Option 2: Using networkx directly
    # G = nx.read_edgelist(filepath, delimiter='\t', data=[('weight', float)])
    
    return G


def get_node_embedding_from_model(model, node_name):
    """
    Extract embedding vector for a specific node from trained Word2Vec model.
    
    Args:
        model: Trained Word2Vec model from Node2Vec
        node_name (str): Name of the node
    
    Returns:
        numpy.ndarray: Embedding vector for the node
    """
    # Word2Vec stores embeddings in model.wv (word vectors)
    if node_name in model.wv:
        return model.wv[node_name]
    else:
        return None


def get_all_embeddings(model, graph):
    """
    Extract embeddings for all nodes in a graph.
    
    Args:
        model: Trained Word2Vec model
        graph: NetworkX graph
    
    Returns:
        tuple: (list of node names, list of embedding vectors)
    """
    nodes = list(graph.nodes())
    embeddings = [model.wv[str(node)] for node in nodes]
    return nodes, embeddings


def find_similar_nodes_by_distance(query_node, all_nodes, all_embeddings, threshold):
    """
    Find nodes similar to a query node based on cosine distance threshold.
    
    Args:
        query_node (str): Name of query node
        all_nodes (list): List of all node names
        all_embeddings (numpy.ndarray): Array of all embeddings (n_nodes x embedding_dim)
        threshold (float): Maximum cosine distance for similarity
    
    Returns:
        set: Set of similar node names
    """
    # Find index of query node
    try:
        query_idx = all_nodes.index(query_node)
    except ValueError:
        return set()
    
    # Get query embedding
    query_embedding = all_embeddings[query_idx].reshape(1, -1)
    
    # Calculate cosine distances to all nodes
    distances = cosine_distances(query_embedding, all_embeddings)[0]
    
    # Find nodes within threshold
    similar_indices = np.where(distances <= threshold)[0]
    similar_nodes = {all_nodes[i] for i in similar_indices}
    
    return similar_nodes


# ============================================================================
# HELPER FUNCTIONS FOR GSEA (Part 2)
# ============================================================================

def read_gmt_file(filepath):
    """
    Read KEGG gene sets from .gmt file format.
    
    Args:
        filepath (str): Path to .gmt file
    
    Returns:
        dict: Dictionary mapping pathway name to list of genes
              {pathway_name: [gene1, gene2, ...]}
    """
    gene_sets = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pathway_name = parts[0]
                # parts[1] is URL/description, skip it
                genes = parts[2:]  # All remaining columns are genes
                gene_sets[pathway_name] = genes
    
    return gene_sets


def calculate_log_fold_change(expression_df, sample_labels, condition1='healthy', condition2='endometriosis'):
    """
    Calculate log fold-change for each gene between two conditions.
    
    Args:
        expression_df (pd.DataFrame): Expression data (genes x samples)
        sample_labels (pd.Series or dict): Maps sample IDs to conditions
        condition1 (str): Name of control condition
        condition2 (str): Name of case condition
    
    Returns:
        pd.Series: Log fold-change for each gene (condition2 - condition1)
    """
    # Get samples for each condition
    condition1_samples = [s for s in expression_df.columns if sample_labels.get(s) == condition1]
    condition2_samples = [s for s in expression_df.columns if sample_labels.get(s) == condition2]
    
    # Calculate mean expression for each condition
    mean_condition1 = expression_df[condition1_samples].mean(axis=1)
    mean_condition2 = expression_df[condition2_samples].mean(axis=1)
    
    # Log fold-change (already log-normalized, so just subtract)
    log_fc = mean_condition2 - mean_condition1
    
    return log_fc


def rank_genes_by_metric(gene_metric_dict, ascending=False):
    """
    Rank genes by a metric (e.g., log fold-change).
    
    Args:
        gene_metric_dict (dict): Dictionary mapping gene names to metric values
        ascending (bool): If False, sort descending (highest first)
    
    Returns:
        list: Ordered list of gene names
    """
    # Convert to pandas Series for easy sorting
    metric_series = pd.Series(gene_metric_dict)
    ranked = metric_series.sort_values(ascending=ascending)
    return list(ranked.index)


def calculate_brownian_bridge_scores(n_in_set, n_not_in_set):
    """
    Calculate the up and down scores for Brownian bridge walk.
    
    Args:
        n_in_set (int): Number of genes in the gene set (that are in expression data)
        n_not_in_set (int): Number of genes NOT in the gene set
    
    Returns:
        tuple: (up_score, down_score)
    """
    # These formulas ensure the walk ends at 0 if genes were randomly distributed
    up_score = 1.0 / n_in_set if n_in_set > 0 else 0
    down_score = 1.0 / n_not_in_set if n_not_in_set > 0 else 0
    
    return up_score, down_score


def brownian_bridge_walk(ranked_genes, gene_set):
    """
    Perform Brownian bridge walk to calculate enrichment score.
    
    Args:
        ranked_genes (list): Ordered list of genes (highest logFC first)
        gene_set (set): Set of genes in the pathway
    
    Returns:
        float: Enrichment score (supremum of the walk)
    """
    # Filter gene set to only genes in ranked list
    gene_set_in_data = set(ranked_genes) & gene_set
    
    if len(gene_set_in_data) == 0:
        return 0.0
    
    # Calculate scores
    n_in_set = len(gene_set_in_data)
    n_not_in_set = len(ranked_genes) - n_in_set
    up_score, down_score = calculate_brownian_bridge_scores(n_in_set, n_not_in_set)
    
    # Walk through ranked list
    running_sum = 0.0
    max_sum = 0.0
    
    for gene in ranked_genes:
        if gene in gene_set_in_data:
            running_sum += up_score
        else:
            running_sum -= down_score
        
        # Track maximum (supremum)
        if running_sum > max_sum:
            max_sum = running_sum
    
    return max_sum


def permute_sample_labels(sample_labels_dict):
    """
    Randomly shuffle condition assignments among samples.
    
    Args:
        sample_labels_dict (dict): Original mapping {sample_id: condition}
    
    Returns:
        dict: Shuffled mapping {sample_id: condition}
    """
    sample_ids = list(sample_labels_dict.keys())
    conditions = list(sample_labels_dict.values())
    
    # Shuffle conditions
    shuffled_conditions = np.random.permutation(conditions)
    
    # Create new mapping
    shuffled_labels = dict(zip(sample_ids, shuffled_conditions))
    
    return shuffled_labels


def bonferroni_correction(p_values, alpha=0.05):
    """
    Apply Bonferroni correction for multiple testing.
    
    Args:
        p_values (dict or pd.Series): Dictionary/Series of p-values
        alpha (float): Significance threshold
    
    Returns:
        dict: Dictionary of {name: is_significant (bool)}
    """
    n_tests = len(p_values)
    corrected_threshold = alpha / n_tests
    
    significant = {name: (pval < corrected_threshold) for name, pval in p_values.items()}
    
    return significant


def calculate_permutation_pvalue(actual_score, permuted_scores):
    """
    Calculate p-value from permutation test.
    
    Args:
        actual_score (float): Actual enrichment score
        permuted_scores (list): List of enrichment scores from permutations
    
    Returns:
        float: p-value
    """
    # Count how many permuted scores are >= actual score
    n_greater_equal = sum(1 for score in permuted_scores if score >= actual_score)
    
    # p-value = proportion of times permuted score >= actual
    p_value = n_greater_equal / len(permuted_scores)
    
    return p_value


# ============================================================================
# TESTING AND DEBUGGING HELPERS
# ============================================================================

def print_graph_info(graph):
    """Print basic statistics about a networkx graph."""
    print(f"Number of nodes: {graph.number_of_nodes()}")
    print(f"Number of edges: {graph.number_of_edges()}")
    print(f"Is connected: {nx.is_connected(graph)}")
    print(f"Number of connected components: {nx.number_connected_components(graph)}")


def print_enrichment_scores_table(scores_dict, top_n=10):
    """
    Print top enrichment scores in a nice table format.
    
    Args:
        scores_dict (dict): Dictionary mapping gene set names to scores
        top_n (int): Number of top scores to display
    """
    # Sort by score
    sorted_scores = sorted(scores_dict.items(), key=lambda x: x[1], reverse=True)
    
    print(f"\nTop {top_n} Enrichment Scores:")
    print("-" * 80)
    print(f"{'Pathway':<50} {'Score':>10}")
    print("-" * 80)
    
    for name, score in sorted_scores[:top_n]:
        print(f"{name:<50} {score:>10.2f}")


def visualize_small_graph(graph, max_nodes=50):
    """
    Visualize a small graph (useful for debugging).
    Requires matplotlib.
    """
    if graph.number_of_nodes() > max_nodes:
        print(f"Graph too large to visualize ({graph.number_of_nodes()} nodes)")
        return
    
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(graph, seed=42)
    nx.draw(graph, pos, with_labels=True, node_color='lightblue', 
            node_size=500, font_size=8, font_weight='bold')
    plt.title("Interaction Network")
    plt.tight_layout()
    plt.show()


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("This file contains helper functions for BIO214 Project 2")
    print("Import these functions in your code or use them as inspiration!")
    print("\nExample function categories:")
    print("  - PPI helpers: Reading data, embeddings, similarity search")
    print("  - GSEA helpers: GMT parsing, log FC, Brownian bridge, permutation testing")
    print("  - Debug helpers: Visualization, printing statistics")

