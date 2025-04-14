#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import json

def parse_args():
    parser = argparse.ArgumentParser(description='Extended Overlap Network Analysis between GRN and GCN')
    parser.add_argument('-g', '--grn', required=True, help='Gene Regulatory Network file (tab-delimited)')
    parser.add_argument('-c', '--gcn', required=True, help='Gene Co-expression Network file (tab-delimited)')
    parser.add_argument('-o', '--output', default='overlap_network', help='Output prefix for files')
    parser.add_argument('-v', '--visualize', action='store_true', help='Generate network visualizations')
    return parser.parse_args()

def read_network_file(file_path, directed=True):
    """Read network file and create appropriate graph"""
    print(f"Reading network from {file_path}...")
    df = pd.read_csv(file_path, sep='\t')
    
    if directed:
        G = nx.DiGraph()
        # Assuming first column is source, second is target
        for _, row in df.iterrows():
            G.add_edge(row[0], row[1])
    else:
        G = nx.Graph()
        # Assuming first column is node1, second is node2
        for _, row in df.iterrows():
            G.add_edge(row[0], row[1])
    
    print(f"Created {'directed' if directed else 'undirected'} graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    return G

def create_extended_overlap_network(G_directed, G_undirected):
    """Create extended overlap network showing all connections of genes in both networks."""
    print("\nCreating extended overlap network...")
    
    # Find overlapping genes
    directed_nodes = set(G_directed.nodes())
    undirected_nodes = set(G_undirected.nodes())
    overlapping_genes = directed_nodes.intersection(undirected_nodes)
    
    print(f"- Total GRN genes: {len(directed_nodes)}")
    print(f"- Total GCN genes: {len(undirected_nodes)}")
    print(f"- Overlapping genes: {len(overlapping_genes)}")
    
    if len(overlapping_genes) == 0:
        print("No overlapping genes found between networks!")
        return None, None, set()
    
    # Create copies of the original networks for the extended overlap
    G_directed_overlap = G_directed.copy()
    G_undirected_overlap = G_undirected.copy()
    
    # Keep only the overlapping genes and their direct connections
    nodes_to_remove_directed = set()
    nodes_to_remove_undirected = set()
    
    # For directed network
    for node in G_directed:
        if node not in overlapping_genes:
            # Check if it's connected to any overlapping gene
            has_connection = False
            for overlap_gene in overlapping_genes:
                if G_directed.has_edge(node, overlap_gene) or G_directed.has_edge(overlap_gene, node):
                    has_connection = True
                    break
            if not has_connection:
                nodes_to_remove_directed.add(node)
    
    # For undirected network
    for node in G_undirected:
        if node not in overlapping_genes:
            # Check if it's connected to any overlapping gene
            has_connection = False
            for overlap_gene in overlapping_genes:
                if G_undirected.has_edge(node, overlap_gene):
                    has_connection = True
                    break
            if not has_connection:
                nodes_to_remove_undirected.add(node)
    
    G_directed_overlap.remove_nodes_from(nodes_to_remove_directed)
    G_undirected_overlap.remove_nodes_from(nodes_to_remove_undirected)
    
    print(f"- Extended GRN overlap: {G_directed_overlap.number_of_nodes()} nodes, {G_directed_overlap.number_of_edges()} edges")
    print(f"- Extended GCN overlap: {G_undirected_overlap.number_of_nodes()} nodes, {G_undirected_overlap.number_of_edges()} edges")
    print(f"- Core overlapping genes: {len(overlapping_genes)}")
    
    return G_directed_overlap, G_undirected_overlap, overlapping_genes

def analyze_overlap_network(G_directed, G_undirected, overlapping_genes, output_prefix):
    """Analyze the overlap network and save results"""
    overlap_df = pd.DataFrame(columns=['Gene', 'GRN_OutDegree', 'GRN_InDegree', 'GCN_Degree'])
    
    for gene in overlapping_genes:
        overlap_df = overlap_df._append({
            'Gene': gene,
            'GRN_OutDegree': G_directed.out_degree(gene),
            'GRN_InDegree': G_directed.in_degree(gene),
            'GCN_Degree': G_undirected.degree(gene)
        }, ignore_index=True)
    
    # Sort by total connectivity (sum of all degrees)
    overlap_df['Total_Connectivity'] = (
        overlap_df['GRN_OutDegree'] + 
        overlap_df['GRN_InDegree'] + 
        overlap_df['GCN_Degree']
    )
    overlap_df = overlap_df.sort_values('Total_Connectivity', ascending=False)
    
    # Save the overlap analysis to file
    output_file = f"{output_prefix}_analysis.tsv"
    overlap_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved overlap network analysis to {output_file}")
    
    # Save network files as edge lists instead of pickle
    grn_edges = f"{output_prefix}_grn_edges.tsv"
    gcn_edges = f"{output_prefix}_gcn_edges.tsv"
    
    # Save GRN (directed)
    with open(grn_edges, 'w') as f:
        f.write("source\ttarget\n")
        for u, v in G_directed.edges():
            f.write(f"{u}\t{v}\n")
            
    # Save GCN (undirected)
    with open(gcn_edges, 'w') as f:
        f.write("node1\tnode2\n")
        for u, v in G_undirected.edges():
            f.write(f"{u}\t{v}\n")
            
    print(f"Saved network edges to {grn_edges} and {gcn_edges}")
    
    # Save overlapping genes list
    with open(f"{output_prefix}_overlap_genes.txt", 'w') as f:
        for gene in sorted(overlapping_genes):
            f.write(f"{gene}\n")
    
    return overlap_df

def visualize_overlap_network(G_directed, G_undirected, overlapping_genes, output_prefix, max_nodes=200):
    """Visualize only the genes that are present in both networks and their direct interactions, 
    limited to a specified number of top nodes by connectivity"""
    print("\nPreparing visualization of overlapping genes only...")
    
    # Create subgraphs containing only the overlapping genes
    overlap_nodes = list(overlapping_genes)
    
    # If we have more nodes than our limit, select the top ones by total degree
    if len(overlap_nodes) > max_nodes:
        print(f"Network has {len(overlap_nodes)} overlapping genes, limiting visualization to top {max_nodes} by connectivity")
        
        # Calculate total degree (connectivity) for each node in both networks
        node_connectivity = {}
        for node in overlap_nodes:
            directed_in_degree = G_directed.in_degree(node) if node in G_directed else 0
            directed_out_degree = G_directed.out_degree(node) if node in G_directed else 0
            undirected_degree = G_undirected.degree(node) if node in G_undirected else 0
            
            node_connectivity[node] = directed_in_degree + directed_out_degree + undirected_degree
        
        # Sort nodes by total connectivity and take the top max_nodes
        overlap_nodes = [node for node, _ in sorted(node_connectivity.items(), 
                                                  key=lambda x: x[1], reverse=True)[:max_nodes]]
        print(f"Selected top {len(overlap_nodes)} genes by connectivity")
    
    # Create subgraphs
    G_directed_overlap = G_directed.subgraph(overlap_nodes).copy()
    G_undirected_overlap = G_undirected.subgraph(overlap_nodes).copy()
    
    print(f"Creating visualization with {len(overlap_nodes)} overlapping genes:")
    print(f"- GRN overlap: {G_directed_overlap.number_of_nodes()} nodes, {G_directed_overlap.number_of_edges()} edges")
    print(f"- GCN overlap: {G_undirected_overlap.number_of_nodes()} nodes, {G_undirected_overlap.number_of_edges()} edges")
    
    # Create combined graph for visualization
    G_combined = nx.DiGraph()
    
    # Add all overlapping nodes
    for node in overlap_nodes:
        G_combined.add_node(node)
    
    # Add edges from both networks with type attribute
    for u, v in G_directed_overlap.edges():
        G_combined.add_edge(u, v, type='regulatory')
    
    for u, v in G_undirected_overlap.edges():
        if not G_combined.has_edge(u, v) and not G_combined.has_edge(v, u):
            G_combined.add_edge(u, v, type='coexpression')
    
    # Create visualization
    print("Generating network layout...")
    plt.figure(figsize=(16, 14))
    
    # Use a layout algorithm appropriate for the network size
    if len(overlap_nodes) > 100:
        print("Using faster layout algorithm...")
        pos = nx.kamada_kawai_layout(G_combined)
    else:
        pos = nx.spring_layout(G_combined, k=0.3, iterations=50, seed=42)
    
    # Draw nodes
    nx.draw_networkx_nodes(G_combined, pos, node_color='lightblue', 
                         node_size=80, alpha=1.0, label='Genes in both networks')
    
    # Draw edges by type
    regulatory_edges = [(u, v) for u, v, d in G_combined.edges(data=True) if d.get('type') == 'regulatory']
    coexp_edges = [(u, v) for u, v, d in G_combined.edges(data=True) if d.get('type') == 'coexpression']
    
    nx.draw_networkx_edges(G_combined, pos, edgelist=regulatory_edges, edge_color='blue', 
                         alpha=0.6, arrows=True, label='Regulatory connection')
    nx.draw_networkx_edges(G_combined, pos, edgelist=coexp_edges, edge_color='green', 
                         alpha=0.4, arrows=False, style='dashed', label='Co-expression')
    
    # Add labels to nodes (if not too many)
    if len(overlap_nodes) <= 50:
        nx.draw_networkx_labels(G_combined, pos, font_size=8)
    else:
        # Identify high-degree nodes for labeling
        node_degrees = dict(G_combined.degree())
        high_degree_nodes = {node: node for node, degree in sorted(node_degrees.items(), 
                                                               key=lambda x: x[1], reverse=True)[:30]}
        nx.draw_networkx_labels(G_combined, pos, labels=high_degree_nodes, font_size=8, font_weight='bold')
    
    plt.title(f"Top {len(overlap_nodes)} Overlapping Genes: Genes Present in Both GRN and GCN", fontsize=16)
    plt.legend(scatterpoints=1, loc='upper right')
    plt.axis('off')
    plt.tight_layout()
    
    output_file = f"{output_prefix}_overlap_only_visualization.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved overlap-only network visualization to {output_file}")
    
    # Also save a list of the overlapping genes that were visualized
    gene_output = f"{output_prefix}_visualized_genes.txt"
    with open(gene_output, 'w') as f:
        for gene in sorted(overlap_nodes):
            f.write(f"{gene}\n")
    print(f"Saved list of {len(overlap_nodes)} visualized genes to {gene_output}")

def output_network_diagnostics(G_directed, G_undirected, overlapping_genes):
    """Output diagnostic information about specific genes in the network"""
    print("\n--- OVERLAP NETWORK DIAGNOSTIC ---")
    
    # Select a few sample genes to diagnose
    sample_genes = list(overlapping_genes)[:5] if len(overlapping_genes) > 5 else list(overlapping_genes)
    
    for gene in sample_genes:
        print(f"\nDiagnosing gene: {gene}")
        print(f"Present in GRN network: {gene in G_directed}")
        print(f"Present in GCN network: {gene in G_undirected}")
        
        if gene in G_directed:
            print(f"GRN in-degree: {G_directed.in_degree(gene)}")
            print(f"GRN out-degree: {G_directed.out_degree(gene)}")
            
            # Show regulators
            regulators = list(G_directed.predecessors(gene))
            if regulators:
                print(f"  Regulated by: {regulators}")
                
            # Show targets
            targets = list(G_directed.successors(gene))
            if targets:
                print(f"  Regulates: {targets}")
        
        if gene in G_undirected:
            print(f"GCN degree: {G_undirected.degree(gene)}")
            
            # Show co-expressed genes
            coexp = list(G_undirected.neighbors(gene))
            if len(coexp) > 5:
                print(f"  Co-expressed with: {coexp[:5]}...")
                print(f"  (and {len(coexp)-5} others)")
            else:
                print(f"  Co-expressed with: {coexp}")
    
    print("\n--- END DIAGNOSTIC ---")

def main():
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)
    
    # Read networks
    G_directed = read_network_file(args.grn, directed=True)
    G_undirected = read_network_file(args.gcn, directed=False)
    
    # Create and analyze extended overlap network
    G_directed_overlap, G_undirected_overlap, overlapping_genes = create_extended_overlap_network(G_directed, G_undirected)
    
    if overlapping_genes:
        # Analyze the overlap network
        overlap_df = analyze_overlap_network(
            G_directed_overlap, G_undirected_overlap, overlapping_genes, args.output
        )
        
        # Output diagnostics
        output_network_diagnostics(G_directed_overlap, G_undirected_overlap, overlapping_genes)
        
        # Visualize if requested
        if args.visualize:
            visualize_overlap_network(G_directed_overlap, G_undirected_overlap, overlapping_genes, args.output)
    # After creating the extended overlap network:
    if args.visualize:
        # Call the new visualization function
        visualize_overlap_network(G_directed, G_undirected, overlapping_genes, args.output, max_nodes=200)
        # Print sample of top genes in overlap
        print("\nTop genes in the overlap network by total connectivity:")
        print(overlap_df[['Gene', 'GRN_OutDegree', 'GRN_InDegree', 'GCN_Degree']].head(10).to_string(index=False))
    
    print("\nExtended overlap network analysis complete!")

if __name__ == "__main__":
    main()