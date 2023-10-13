

def debruijn_graph(string, k_mer):
    # Initialize lists to store edges and a set to store nodes (k-1 mers)
    edges = []
    nodes = set()

    # Iterate through the input string to create k-mers and their overlapping edges
    for i in range(len(string) - k_mer + 1):
        # Extract a k-mer by taking a slice from the input string
        kmer = string[i: i + k_mer]
        
        # Create an edge by taking the first k-1 characters as the source node
        # and the last k-1 characters as the target node
        edge = (kmer[:k_mer - 1], kmer[1:])
        
        # Add the edge to the list of edges
        edges.append(edge)
        
        # Add the source and target nodes to the set of nodes
        nodes.add(kmer[:k_mer - 1])
        nodes.add(kmer[1:])

    # Return the list of edges and the set of nodes
    return edges, nodes

# Example usage: Generate the De Bruijn graph for the input 'ACGCGTCG' with k-mer size 3
nodes, edges = debruijn_graph('ACGCGTCG', 3)

# 'nodes' will contain the set of nodes (k-1 mers)
# 'edges' will contain the list of edges connecting these nodes















    


