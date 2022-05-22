from typing import Tuple, List
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import string


def read_file(path: str) -> np.array:
    matrix = np.loadtxt(path, dtype=int)
    return matrix


def check_filename(filename: str) -> str:
    index = filename.rfind('.')
    if index == -1 or not filename.endswith('.txt'):
        print('Not proper name of file')
        exit(0)

    return filename[:index]


def weighted_adjacency_list(filename: str, graph: nx.Graph) -> None:
    adjacency_dict, output_list = nx.to_dict_of_lists(graph), []
    for k, value in adjacency_dict.items():
        for val in value:
            output_list.append(str(k) + '->' + str(val) + ':' + str(int(graph[k][val]['weight'])))

    return print(f'{filename}\n{sorted(output_list)}')


def save_result(graph: nx.Graph, path: str) -> None:
    position = nx.spring_layout(graph, seed=11)
    labels = nx.get_edge_attributes(graph, 'weight')

    plt.figure()
    nx.draw(
        graph, position, edge_color='black', width=1, linewidths=1,
        node_size=500, node_color='pink', alpha=0.9,
        labels={node: node for node in graph.nodes()}
    )
    nx.draw_networkx_edge_labels(
        graph, position,
        edge_labels=labels,
        font_color='green'
    )

    plt.axis("off")
    plt.savefig(path, format='png')
    plt.close()
