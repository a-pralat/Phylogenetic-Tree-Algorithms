import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def limb_length(matrix: np.array, j: int) -> float:
    limb, matrix_size = np.inf, np.size(matrix, 0)
    mask = [i for i in range(matrix_size) if i != j]
    for k in range(matrix_size):
        if k != j:
            return min(limb, min(matrix[mask, j] + matrix[j, k] - matrix[mask, k])) / 2


def find_attachment_point(matrix: np.array, j: int) -> tuple[int, int]:
    # i_out, k_out = 0, 0
    mask = [i for i in range(np.size(matrix, 0)) if i != j]
    for u in range(len(mask) - 1):
        for v in range(u + 1, len(mask)):
            i, k = mask[u], mask[v]
            if matrix[i][j] + matrix[j][k] == matrix[i][k]:
                return i, k
    # return i_out, k_out


def additive_phylogeny(matrix: np.array, n: int, inner_n: int) -> nx.Graph:
    graph = nx.Graph()

    if n == 2:
        graph.add_edge(0, 1, weight=matrix[0][1])
        return graph

    # 1. Pick an arbitrary leaf j (we choose THE LAST ONE)
    j = n - 1
    # 2. Compute its limb length, limb_length(matrix, j)
    limb = limb_length(matrix, j)

    # 3. Subtract limb_length(matrix, j) from each row and column to
    # produce matrix_bald in which j is a bald limb (length 0)
    matrix[j, :], matrix[:, j], matrix[j, j] = matrix[j, :] - limb, matrix[:, j] - limb, 0

    # ?. (to 6) i, j, k leaves such that matrix_bald[i,j] + matrix_bald[j,k] = matrix_bald[i,k]
    i, k = find_attachment_point(matrix, j)
    matrix_i_j = matrix[i, j]

    # 4. Remove the j-th row and column of the matrix to form the (n-1) x (n-1) matrix_trim
    matrix = np.delete(matrix, j, 0)
    matrix = np.delete(matrix, j, 1)

    # ?. Necessary counter to rebuild the tree
    while inner_n in list(graph.nodes):
        inner_n += 1

    # ?. Recursion
    result = additive_phylogeny(matrix, n - 1, inner_n)

    # 5. Construct tree (matrix_trim)
    # 6. Identify the point in Tree (matrix_trim) where leaf j should be attached
    # 7. Attach j by an edge of length limb_length(matrix, j) in order to form Tree(matrix)
    new_node = None
    path, distance = nx.shortest_path(result, source=i, target=k), 0
    for i in range(1, len(path) - 1):
        distance += result[path[i - 1]][path[i]]['weight']
        new_node = path[i] if distance == matrix_i_j else None

    if new_node is None:
        new_node = inner_n
        while new_node in list(result.nodes):
            new_node += 1
        distance, small_distance, i = 0, 0, 0
        while distance < matrix_i_j:
            i += 1
            small_distance = distance
            distance += result[path[i - 1]][path[i]]['weight']

        result.remove_edge(path[i - 1], path[i])
        result.add_edge(new_node, path[i], weight=distance - matrix_i_j)
        result.add_edge(new_node, path[i - 1], weight=matrix_i_j - small_distance)

    result.add_edge(new_node, n - 1, weight=limb)

    return result


def read_file(path: str) -> np.array:
    return np.loadtxt(path, dtype=int)


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
    nx.draw_networkx_nodes(graph, position, node_size=200)
    nx.draw_networkx_edges(graph, position, width=2)
    labels = nx.get_edge_attributes(graph, 'weight')
    nx.draw_networkx_edge_labels(graph, position, edge_labels=labels)
    nx.draw_networkx_labels(graph, position, font_size=10, font_family="sans-serif")

    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, format='png')
    plt.close()


if __name__ == "__main__":
    input_dir, output_dir = 'examples/input/', 'examples/output/'
    files = ['n4_1.txt', 'n4_2.txt', 'n5_1.txt', 'n8_1.txt', 'n10_1.txt', 'n10_2.txt',
             'n10_3.txt']

    for f in files:
        input_path = input_dir + f
        file = check_filename(f)

        D = read_file(input_path)
        size = np.size(D, 0)

        output = additive_phylogeny(D, size, size)
        weighted_adjacency_list(file, output)
        save_result(output, f'{output_dir}{file}.png')
