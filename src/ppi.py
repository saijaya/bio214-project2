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
        self.interaction_network = None
        self.graph = None
        # Add any other instance variables you need

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

        self.interaction_network = load_interaction_network(interactionNetworkFile)
        nodes = list(self.interaction_network.nodes())
        print("\ngene_network:")
        print(f"Total nodes: {len(nodes)}")
        print(f"Sample nodes: {nodes[:20]}")
        edges = list(self.interaction_network.edges())
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
        
        Performance tip: Save and load pre-trained model to speed up autograder.
        See project description lines 89-108 for details.
        """
        
        # Performance optimization: check if pre-trained model exists
        if os.path.exists("node2vec_pretrained"):
            model = Word2Vec.load("node2vec_pretrained")
            print("loaded pretrained model successfully")
        else:
            # TODO: Create Node2Vec object with EXACT parameters
            # node2vec = Node2Vec(self.graph, dimensions=64, walk_length=30, 
            #                     num_walks=100, workers=1, seed=42)
            
            # TODO: Fit the model with EXACT training parameters
            # model = node2vec.fit(window=3, min_count=1, batch_words=4)
            
            # TODO: Save the model for future use
            # model.save("node2vec_pretrained")
            
            pass
        
        # TODO: Extract embeddings for all nodes in the graph
        # Hint: model.wv contains the word vectors (node embeddings)
        # Make sure gene_nodes and gene_embeddings are in corresponding order!
        
        gene_nodes = []  # List of gene names
        gene_embeddings = []  # List of embedding vectors
        
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
        
        # TODO: Convert gene_embeddings to numpy array if needed
        
        # TODO: For each disease gene:
        #       1. Find its embedding
        #       2. Calculate cosine distances to all genes
        #       3. Find genes within threshold distance
        
        # TODO: Take union of similar genes across all disease genes
        
        # TODO: Include original disease genes in the result set
        
        similar_genes = set()
        
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
    # gene_nodes, gene_embeddings = ppi.calculate_embedding()
    
    # Example: Get similar genes with threshold of 0.5
    # threshold = 0.5
    # similar_genes = ppi.get_close_genes(gene_nodes, gene_embeddings, threshold)
    
    # print(f"Found {len(similar_genes)} similar genes (including disease genes)")
    # You can output to file here if needed for your own analysis


if __name__ == "__main__":
    main()

