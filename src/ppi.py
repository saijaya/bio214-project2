"""
BIO214 Project 2 - Part 1: Protein-Protein Interaction Network Analysis
Author: [Your Name]
Date: [Date]
Collaborators: [List collaborators or state "None"]

This script implements Node2Vec embedding for predicting disease-related genes
based on protein-protein interaction networks.

Usage: python3 ppi.py diseaseGeneFile interactionNetworkFile
"""

import sys
import os
import networkx as nx
from node2vec import Node2Vec
from gensim.models import Word2Vec
from sklearn.metrics.pairwise import cosine_distances
import numpy as np
import pandas as pd


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
    df = pd.read_csv(filepath, sep=' ', header=None, names=['gene1', 'gene2', 'weight'])
    print("\nInput gene network file:")
    print(df)
    gene_network = nx.from_pandas_edgelist(df, source='gene1', target='gene2')
    return gene_network


class PPI:
    """
    Class to perform protein-protein interaction network analysis using Node2Vec embeddings.
    """
    
    def __init__(self):
        """
        Initialize the PPI class.
        Store any instance variables you need here.
        """
        self.disease_genes = None
        self.gene_network = None
        self.trained_node2vec_model = None

    def load_data(self, diseaseGeneFile, interactionNetworkFile):
        """
        Load disease gene list and interaction network from files.
        
        Args:
            diseaseGeneFile (str): Path to file containing disease genes (one per line)
            interactionNetworkFile (str): Path to tab-delimited interaction network file
                                         (gene1, gene2, weight per line)
        
        Returns:
            None (stores data in instance variables)
        
        Hint: Use networkx to create a graph from the interaction network.
              Consider using nx.read_edgelist() or nx.parse_edgelist()
        """
        self.disease_genes = read_gene_list(diseaseGeneFile)
        print("\ndisease_genes:")
        for gene in self.disease_genes:
            print(gene)

        self.gene_network = load_interaction_network(interactionNetworkFile)
        nodes = list(self.gene_network.nodes())
        print("\ngene_network:")
        print(f"Total nodes: {len(nodes)}")
        print(f"Sample nodes: {nodes[:20]}")
        edges = list(self.gene_network.edges())
        print(f"Total edges: {len(edges)}")
        print(f"Sample edges: {edges[:10]}")

    def calculate_embedding(self):
        """
        Calculate node embeddings using Node2Vec algorithm.
        
        REQUIRED PARAMETERS:
        - Node2Vec: dimensions=64, walk_length=30, num_walks=100, workers=1, seed=42
        - Training: window=3, min_count=1, batch_words=4
        
        Returns:
            tuple: (gene_nodes, gene_embeddings)
                - gene_nodes: list of node names (strings)
                - gene_embeddings: list of corresponding vector embeddings (numpy arrays)
                Order must correspond between the two lists!
        
        Performance tip: Save and load pre-trained trained_node2vec_model to speed up autograder.
        See project description lines 89-108 for details.
        """
        
        # Performance optimization: check if pre-trained trained_node2vec_model exists
        if os.path.exists("node2vec_pretrained"):
            trained_node2vec_model = Word2Vec.load("node2vec_pretrained")
            print("loaded pretrained trained_node2vec_model successfully")
        else:
            # Create Node2Vec object
            # dimensions = 64, walk_length = 30, num_walks = 100, workers = 1, seed = 42
            node2vec_model = Node2Vec(self.gene_network, dimensions=64, walk_length=30,
                                      num_walks=100, workers=1, seed=42)
            
            # Fit the trained_node2vec_model with EXACT training parameters
            # window=3, min_count=1, batch_words=4
            trained_node2vec_model = node2vec_model.fit(window=3, min_count=1, batch_words=4)
            
            # Save the trained_node2vec_model for future use
            trained_node2vec_model.save("node2vec_pretrained")

        self.trained_node2vec_model = trained_node2vec_model
        
        # Extract embeddings for all nodes in the graph
        gene_nodes = list(self.gene_network.nodes())
        print(f"num genes in network: {len(gene_nodes)}")
        gene_embeddings = [trained_node2vec_model.wv[node] for node in gene_nodes]
        print(f"num embeddings: {len(gene_embeddings)}")
        print(f"verify embedding 0: {self.trained_node2vec_model.wv[gene_nodes[0]]}")
        print(f"verify embedding 0: {gene_embeddings[0]}")

        return gene_nodes, gene_embeddings

    def get_close_genes(self, gene_nodes, gene_embeddings, threshold):
        """
        Find genes similar to disease genes based on cosine distance threshold.
        
        Args:
            gene_nodes (list): List of all gene names (strings)
            gene_embeddings (list): List of all gene embeddings (numpy arrays)
                                   Must be in same order as gene_nodes
            threshold (float): Maximum cosine distance for genes to be considered similar
                             Valid range: [0, 2]
        
        Returns:
            set: Set of gene names (strings) that are similar to disease genes
                 Includes the original disease genes
                 If no genes are predicted, return just the disease genes
        
        Hint: Use sklearn.metrics.pairwise.cosine_distances (much faster than scipy)
              For each disease gene, find all genes within threshold distance
              Take union across all disease genes
        """
        
        # Convert gene_embeddings to numpy array
        network_gene_embeddings_matrix = np.array(gene_embeddings)
        similar_genes = set(self.disease_genes)
        gene_to_index_dict = {gene: index for index, gene in enumerate(gene_nodes)}

        disease_gene_indices = [gene_to_index_dict[gene] for gene in self.disease_genes if gene in gene_to_index_dict]

        if len(disease_gene_indices) == 0:
            return similar_genes

        disease_gene_embeddings_matrix = network_gene_embeddings_matrix[disease_gene_indices]
        disease_to_network_gene_distances = cosine_distances(disease_gene_embeddings_matrix, network_gene_embeddings_matrix)

        threshold_filter = (disease_to_network_gene_distances <= threshold).any(axis=0)
        gene_nodes_matrix = np.array(gene_nodes)
        close_genes = gene_nodes_matrix[threshold_filter]
        similar_genes.update(close_genes)

        print(f"Found {len(similar_genes)} similar genes (including {len(self.disease_genes)} disease genes)")

        return similar_genes


def main():
    """
    Main function to run PPI analysis from command line.
    This is for your own testing - autograder calls class methods directly.
    """
    if len(sys.argv) != 3:
        print("Usage: python3 ppi.py diseaseGeneFile interactionNetworkFile")
        sys.exit(1)
    
    disease_file = sys.argv[1]
    network_file = sys.argv[2]
    
    # Create PPI object and run analysis
    ppi = PPI()
    ppi.load_data(disease_file, network_file)
    gene_nodes, gene_embeddings = ppi.calculate_embedding()

    # Example: Get similar genes with threshold of 0.5
    threshold = 0.2
    similar_genes = ppi.get_close_genes(gene_nodes, gene_embeddings, threshold)
    print(similar_genes)
    
    # You can output to file here if needed for your own analysis


if __name__ == "__main__":
    main()

