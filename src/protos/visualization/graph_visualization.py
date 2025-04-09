import networkx as nx
import matplotlib.pyplot as plt


def visualize_all_edges(all_edges, avg_connectivity=None, name=''):
    """
    Visualize the all-edges dict as a network graph.

    Parameters:
    all_edges (set): Set of tuples representing edges between residues.
    avg_connectivity (dict, optional): Dictionary of average connectivity for each residue.
    """
    # Create a graph
    G = nx.Graph()

    # Add edges to the graph
    G.add_edges_from(all_edges)

    # Set up the plot
    plt.figure(figsize=(12, 8))

    # Compute layout
    pos = nx.spring_layout(G, k=0.5, iterations=50)

    # Draw nodes
    if avg_connectivity:
        node_sizes = [avg_connectivity.get(node, 1) * 100 for node in G.nodes()]
        node_colors = [avg_connectivity.get(node, 0) for node in G.nodes()]
        node_collection = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                                                 cmap=plt.cm.Greys)
    else:
        node_collection = nx.draw_networkx_nodes(G, pos, node_size=100)

    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=0.3)

    # Draw labels with white color for better visibility
    nx.draw_networkx_labels(G, pos, font_size=8, font_color='red')

    # Add colorbar if using avg_connectivity
    if avg_connectivity:
        sm = plt.cm.ScalarMappable(cmap=plt.cm.Greys, norm=plt.Normalize(vmin=min(avg_connectivity.values()),
                                                                         vmax=max(avg_connectivity.values())))
        sm.set_array([])

    plt.title(f"Graph of Retinal Binding Domain in {name}")
    plt.axis('off')
    plt.tight_layout()
    return plt