from utils import *


class UPGMA:
    def __init__(self):
        self.all_labels: List = []
        self.node_counter: int = 0

    def init_labels(self, size: int) -> List:
        alphabet = string.ascii_uppercase
        if size > len(alphabet):
            raise ValueError('[ERROR] Max number of labels (nodes) is 26')

        self.all_labels = list(alphabet[:size])
        return list(alphabet[:size])

    @staticmethod
    def closest_clusters(matrix: np.array) -> Tuple[int, int]:
        minimum_distance, size = np.inf, np.size(matrix, 0)
        x, y = np.NINF, np.NINF

        for i in range(size):
            for j in range(i):
                if matrix[i][j] < minimum_distance:
                    minimum_distance = matrix[i][j]
                    x, y = i, j

        return x, y

    def merge_labels(self, labels: List, x: int, y: int) -> str:
        if x > y:
            x, y = y, x

        labels[x] = labels[x] + labels[y]

        self.all_labels.append(labels[x])
        self.node_counter += 1

        labels.pop(y)

        return labels[x]

    @staticmethod
    def update_matrix(matrix: np.array, x: int, y: int) -> np.array:
        size = np.size(matrix, 0)
        new_matrix = np.zeros((size, size))

        if x > y:
            x, y = y, x

        for i in range(size):
            for j in range(i):
                if i not in [x, y] and j not in [x, y]:
                    new_matrix[i][j] = matrix[i][j]
                    new_matrix[j][i] = new_matrix[i][j]  # symmetrical

        for i in range(size):
            if i != x:
                new_matrix[x][i] = (matrix[x][i] + matrix[y][i]) / 2
                new_matrix[i][x] = new_matrix[x][i]  # symmetrical

        new_matrix = np.delete(new_matrix, y, 0)
        new_matrix = np.delete(new_matrix, y, 1)

        return new_matrix

    @staticmethod
    def build_tree(result: dict) -> nx.Graph:
        graph = nx.Graph()

        for key1, value1 in result.items():
            for key2, value2 in result.items():
                if key1 != key2 and key1 + key2 in result.keys():
                    graph.add_edge(key1, key1 + key2, weight=round(result[key1 + key2] - result[key1], 2))
                    graph.add_edge(key2, key1 + key2, weight=round(result[key1 + key2] - result[key2], 2))

        return graph

    def upgma(self, matrix: np.array) -> nx.Graph:
        result, self.node_counter = dict(), np.size(matrix, 0)

        # 1. Form a cluster for each present-dat species, each containing a single leaf
        labels = self.init_labels(np.size(matrix, 0))

        for label in labels:
            result[label] = 0

        # 6. Iterate until a single cluster contains all species
        while len(labels) > 1:
            # 2. Find the two closest clusters C1 and C2 according to the average distance
            x, y = self.closest_clusters(matrix)

            # 3. Merge C1 and C2 into a single cluster C
            new_label = self.merge_labels(labels, x, y)

            # 4. Form a new node for C and connect to C1 and C2 by an edge. Set age of C as D_avg(C1, C2) / 2
            result[new_label] = matrix[x][y] / 2

            # 5. Update the distance matrix by computing the average distance between each pair of clusters
            matrix = self.update_matrix(matrix, x, y)

        return self.build_tree(result)


if __name__ == "__main__":
    # files configuration
    input_dir, output_dir = 'examples/input/', 'examples/output/'
    files = ['n4_1.txt', 'n5_1.txt', 'n5_2.txt', 'n5_3.txt', 'n6_1.txt', 'n8_1.txt', 'n9_1_additive.txt', 'n9_2_nonadditive.txt']

    for f in files:
        input_path, file = input_dir + f, check_filename(f)

        # execute
        upgma = UPGMA()
        output = upgma.upgma(read_file(input_path))

        # show output
        weighted_adjacency_list(file, output)

        # choose visualization
        # save_result(output, f'{output_dir}UPGMA/graph_{file}.png')
        save_result(output, f'{output_dir}UPGMA/tree_{file}.png', 'tree')

# Sources:
# https://www.youtube.com/watch?v=27aOeJ2hSwA
