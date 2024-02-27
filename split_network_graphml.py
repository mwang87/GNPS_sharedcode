#!/usr/bin/python

import sys
import os
import networkx as nx
import argparse
import pandas as pd

def main():
    # define args
    parser = argparse.ArgumentParser(description='Creating Component Summary')
    parser.add_argument('input_graphml', help='input_graphml')
    parser.add_argument('output_component_summary', help='output_component_summary')
    parser.add_argument('output_component_graphml_folder', help='output_component_graphml_folder')
    args = parser.parse_args()

    # reading the graphml
    G = nx.read_graphml(args.input_graphml)

    # getting the components
    components = list(nx.connected_components(G))

    component_list = []

    for i, component in enumerate(components):
        subgraph = G.subgraph(component)
        component_id = i + 1

        output_component_filename = os.path.join(args.output_component_graphml_folder, f"component_{component_id}.graphml")
        nx.write_graphml(subgraph, output_component_filename)

        component_dict = {}
        component_dict["component_id"] = component_id
        component_dict["number_of_nodes"] = len(component)
        component_dict["number_of_edges"] = subgraph.number_of_edges()

        component_list.append(component_dict)

    # writing the component summary
    component_df = pd.DataFrame(component_list)
    component_df.to_csv(args.output_component_summary, sep="\t", index=False)







if __name__ == "__main__":
    main()

