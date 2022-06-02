from utils import *


class NeighborJoining:
    def __init__(self, matrix: np.array):
        self.matrix, self.size = matrix, np.size(matrix, 0)
        self.alphabet = string.ascii_uppercase
        self.labels = self.init_labels()
        self.result = {}

    def init_labels(self) -> List[str]:
        if self.size > len(self.alphabet):
            raise ValueError('[ERROR] Max number of labels (nodes) is 26')

        return list(self.alphabet[:self.size])

    def find_minimum_element(self) -> Tuple[int, int]:
        row_sum = np.zeros(self.size)
        for i in range(self.size):
            row_sum[i] = np.sum(self.matrix[i])

        minimum, indexes = np.inf, (-1, -1)
        for i in range(self.size):
            for j in range(i):
                value = (self.size - 2) * self.matrix[i][j] - row_sum[i] - row_sum[j]

                if minimum > value:
                    minimum, indexes = value, (i, j)

        return indexes

    def update_matrix(self, i: int, j: int) -> List[int]:
        remaining_idx = [k for k in range(self.size) if k not in {i, j}]
        new_matrix = np.zeros((self.size - 1, self.size - 1))
        new_matrix[:-1, :-1] = self.matrix[np.ix_(remaining_idx, remaining_idx)]

        self.matrix, self.size = new_matrix, self.size - 1
        return remaining_idx

    @staticmethod
    def check_if_leaf(value: str) -> bool:
        for letter in value:
            if letter not in string.digits:
                return True
        return False

    @staticmethod
    def add_edges_to_result(graph: nx.Graph, first: str, second: str, first_val: float, second_val: float) -> None:
        graph.add_edge(first + second, first, weight=round(first_val, 2))
        graph.add_edge(first + second, second, weight=round(second_val, 2))

    def build_tree(self, result: dict) -> nx.Graph:
        graph, merged = nx.Graph(), dict()

        for key, val in result.items():
            first = val[0][0] if self.check_if_leaf(val[0][0]) else merged.get(val[0][0])
            second = val[1][0] if self.check_if_leaf(val[1][0]) else merged.get(val[1][0])
            self.add_edges_to_result(graph, first, second, val[0][1], val[1][1])
            merged[key] = first + second

        return graph

    def neighbor_joining(self) -> nx.Graph:
        # Init - increasing numbers used as ids of newly created internal nodes, starting at 1
        number = iter(list(range(1, 1000)))
        nodes = self.labels
        while self.size > 2:
            temp_nodes = []

            # 1. Construct neighbor-join matrix matrix_star from matrix
            # 2. Find a minimum element matrix_star[i][j] in matrix_star
            i, j = self.find_minimum_element()

            # 3. Compute delta_i_j = (TotalDistance_D(i) - TotalDistance_D(j)) / (n-2)
            delta_i_j = (np.sum(self.matrix[i]) - np.sum(self.matrix[j])) / (self.size - 2)

            # 5. Form a matrix by removing i-th and j-th row/column from matrix and adding
            # an m-th row/column such that for any k: D_k,m = (D_i,k + D_j,k - D_i,j) / 2
            matrix_copy = np.copy(self.matrix)

            remaining_idx = self.update_matrix(i, j)
            k = self.size - 1
            for m, m_old in enumerate(remaining_idx):
                distance = (matrix_copy[i][m_old] + matrix_copy[j][m_old] - matrix_copy[i][j]) / 2
                self.matrix[k][m], self.matrix[m][k] = distance, distance
                temp_nodes.append(nodes[m_old])

            # 4. Set LimbLength(i) equal to 1/2(D_i,j + delta_i_j) and LimbLength(j) equal to 1/2(D_i,j - delta_i_j)
            limb_length_i_k = 0.5 * (matrix_copy[i][j] + delta_i_j)
            limb_length_j_k = 0.5 * (matrix_copy[i][j] - delta_i_j)

            # 6. Apply NeighborJoining to D' to obtain Tree(D')
            # get an id for the new node
            node = str(next(number))
            # add the new node to the list of remaining nodes
            temp_nodes.append(node)
            # add k -> i,j edges to tree
            self.result[node] = []
            self.result[node].append((nodes[i], limb_length_i_k))
            self.result[node].append((nodes[j], limb_length_j_k))
            # update list of remaining nodes
            nodes = temp_nodes
        # only two nodes left
        # normally we would add an edge connecting those two remaining nodes (unrooted tree)
        # but since we want a rooted tree, we are adding a root node between them
        # id of the root node
        node = str(next(number))
        # add root node to tree
        self.result[node] = []
        # edges from root node to remaining nodes
        # each edge has length of the distance between them divided by two
        self.result[node].append((nodes[0], self.matrix[0][1] / 2))
        self.result[node].append((nodes[1], self.matrix[0][1] / 2))

        return self.build_tree(self.result)


if __name__ == "__main__":
    # files configuration
    input_dir, output_dir = 'examples/input/', 'examples/output/'
    files = ['n4_1.txt', 'n4_2.txt', 'n5_1.txt', 'n8_1.txt', 'n9_1_additive.txt', 'n9_2_nonadditive.txt']

    for f in files:
        input_path, file = input_dir + f, check_filename(f)

        # execute
        nj = NeighborJoining(read_file(input_path))
        output = nj.neighbor_joining()

        # show output
        weighted_adjacency_list(file, output)

        # choose visualization
        save_result(output, f'{output_dir}neighbor-joining/graph_{file}.png')
        save_result(output, f'{output_dir}neighbor-joining/tree_{file}.png', 'tree')

# Sources:
# https://www.youtube.com/watch?v=Y0QWFFWQzds
