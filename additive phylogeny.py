from utils import *


def limb_length(matrix: np.array, j: int) -> float:
    limb, matrix_size = np.inf, np.size(matrix, 0)
    mask = [i for i in range(matrix_size) if i != j]
    for k in range(matrix_size):
        if k != j:
            return min(limb, min(matrix[mask, j] + matrix[j, k] - matrix[mask, k])) / 2


def find_attachment_point(matrix: np.array, j: int) -> Tuple[int, int]:
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


if __name__ == "__main__":
    input_dir, output_dir = 'examples/input/', 'examples/output/'
    files = ['n4_1.txt', 'n4_2.txt', 'n5_1.txt', 'n8_1.txt', 'n9_1_additive.txt']

    for f in files:
        input_path = input_dir + f
        file = check_filename(f)

        D = read_file(input_path)
        size = np.size(D, 0)

        output = additive_phylogeny(D, size, size)
        weighted_adjacency_list(file, output)
        save_result(output, f'{output_dir}additive phylogeny/{file}.png')
