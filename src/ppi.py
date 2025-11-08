"""
disclaimer:
I used cursor to help me with some syntax (such as how to melt dataframes etc.)
the same way I would have used stackoverflow in the past.
However, the choice of data structures, implementation and logic is entirely my own.
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
    read list of genes from a file
    input:
        filepath (str): Path to gene list file
    return:
        list: list of gene names
    """
    genes = []
    with open(filepath, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.append(gene)
    return genes


def load_interaction_network(filepath):
    """
    load interaction network into a networkx obj
    input:
        filepath (str): network file path
    return:
        networkx.Graph obj
    """
    df = pd.read_csv(filepath, sep=' ', header=None, names=['gene1', 'gene2', 'weight'])
    gene_network = nx.from_pandas_edgelist(df, source='gene1', target='gene2')
    return gene_network


class PPI:
    """
    class to perform protein-protein interaction network analysis
    """
    
    def __init__(self):
        """
        initalize the PPI class.
        """
        self.disease_genes = None
        self.gene_network = None
        self.trained_node2vec_model = None

    def load_data(self, disease_gene_file, interaction_network_file):
        """
        load disease gene list and interaction network from files
        input:
            disease_gene_file (str): disease genes file path
            interaction_network_file (str): interaction network file path
        return:
            None
        """
        self.disease_genes = read_gene_list(disease_gene_file)
        self.gene_network = load_interaction_network(interaction_network_file)

    def calculate_embedding(self):
        """
        calculate node embeddings
        return:
            tuple: (gene_nodes, gene_embeddings)
                - gene_nodes: list of node names (strings)
                - gene_embeddings: list of corresponding vector embeddings (numpy arrays)
        """
        
        # Performance optimization: check if pre-trained trained_node2vec_model exists
        if os.path.exists("node2vec_pretrained"):
            trained_node2vec_model = Word2Vec.load("node2vec_pretrained")
        else:
            # create node2vec object
            # dimensions = 64, walk_length = 30, num_walks = 100, workers = 1, seed = 42
            node2vec_model = Node2Vec(self.gene_network, dimensions=64, walk_length=30,
                                      num_walks=100, workers=1, seed=42)
            
            # fit the trained_node2vec_model with training parameters
            # window=3, min_count=1, batch_words=4
            trained_node2vec_model = node2vec_model.fit(window=3, min_count=1, batch_words=4)
            
            # save trained_node2vec_model
            trained_node2vec_model.save("node2vec_pretrained")

        self.trained_node2vec_model = trained_node2vec_model
        
        # get embeddings for all nodes in graph
        gene_nodes = list(self.gene_network.nodes())
        gene_embeddings = [trained_node2vec_model.wv[node] for node in gene_nodes]

        return gene_nodes, gene_embeddings

    def get_close_genes(self, gene_nodes, gene_embeddings, threshold):
        """
        fid genes similar to disease genes
        input:
            gene_nodes (list): list of all gene names (strings)
            gene_embeddings (list): List of all gene embeddings (numpy arrays)
            threshold (float): max cosine distance for genes to be considered similar
        return:
            set: set of gene names that are similar to disease genes
        """
        
        # convert gene_embeddings to numpy array
        network_gene_embeddings_matrix = np.array(gene_embeddings)
        similar_genes = set(self.disease_genes)
        gene_to_index_dict = {gene: index for index, gene in enumerate(gene_nodes)}
        disease_gene_indices = [gene_to_index_dict[gene] for gene in self.disease_genes if gene in gene_to_index_dict]

        if len(disease_gene_indices) == 0:
            return similar_genes

        # EFFICIENCY FIX: stacking all disease gene embeddings into a single matrix. This makes it significantly faster:
        # going from loops to using numpy's matrix calculation efficiency.
        disease_gene_embeddings_matrix = network_gene_embeddings_matrix[disease_gene_indices]
        disease_to_network_gene_distances = cosine_distances(disease_gene_embeddings_matrix, network_gene_embeddings_matrix)

        threshold_filter = (disease_to_network_gene_distances <= threshold).any(axis=0)
        gene_nodes_matrix = np.array(gene_nodes)
        close_genes = gene_nodes_matrix[threshold_filter]
        similar_genes.update(close_genes)

        return similar_genes


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 ppi.py diseaseGeneFile interactionNetworkFile")
        sys.exit(1)
    
    disease_file = sys.argv[1]
    network_file = sys.argv[2]
    
    ppi = PPI()
    ppi.load_data(disease_file, network_file)
    gene_nodes, gene_embeddings = ppi.calculate_embedding()

    threshold = 0.2
    similar_genes = ppi.get_close_genes(gene_nodes, gene_embeddings, threshold)
    print(similar_genes)


if __name__ == "__main__":
    main()

