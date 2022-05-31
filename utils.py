from networkx.drawing.nx_pydot import graphviz_layout
from typing import Tuple, List
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import string


def read_file(path: str) -> np.array:
    with open(path, 'r') as file:
        content = [[float(x) for x in line.split()] for line in file]

    matrix = np.zeros((len(content), len(content)), dtype=float)
    for i in range(len(content)):
        for j in range(len(content[i])):
            matrix[i][j], matrix[j][i] = content[i][j], content[i][j]
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


def save_result(graph: nx.Graph, path: str, mode: str = 'graph') -> None:
    position = nx.spring_layout(graph, seed=11) if mode == 'graph' else graphviz_layout(graph, prog="twopi")
    labels = nx.get_edge_attributes(graph, 'weight')

    plt.figure()
    nx.draw(
        graph, position, edge_color='black', width=1, linewidths=1,
        node_size=300, node_color='pink', alpha=0.9,
        labels={node: node for node in graph.nodes()}
    )
    nx.draw_networkx_edge_labels(
        graph, position,
        edge_labels=labels,
        font_size=8,
        font_color='green'
    )

    plt.axis("off")
    plt.savefig(path, format='png')
    plt.close()
