from utils import *


class AdditivePhylogeny:
    def __init__(self):
        self.new_nodes_start_index = 30

    @staticmethod
    def compute_limb_length(d_ik: float, d_ij: float, d_jk: float) -> float:
        return (d_ik + d_ij - d_jk) / 2

    @staticmethod
    def create_mask(matrix: np.array) -> Tuple[np.array, int]:
        mask, mask_size = [], np.size(matrix, 0)

        # takes half of the matrix
        for i in range(mask_size):
            for j in range(i + 1, mask_size):
                mask.append([i, j, matrix[i][j]])

        return np.array(mask), int(mask_size * (mask_size - 1) / 2)

    def compute_min_limb_length(self, matrix: np.array, j: int) -> float:
        min_limb = np.inf
        mask, mask_size = self.create_mask(matrix)

        for i in range(mask_size):
            x, y, value = int(mask[i, 0]), int(mask[i, 1]), int(mask[i, 2])

            limb = self.compute_limb_length(matrix[x, j], matrix[y, j], value)
            if x != j and y != j and min_limb > limb:
                min_limb = limb

        return min_limb

    def find_attachment_point(self, matrix: np.array, j: int) -> Tuple[int, int]:
        mask, mask_size = self.create_mask(matrix)

        for i in range(mask_size):
            x, y, value = int(mask[i, 0]), int(mask[i, 1]), int(mask[i, 2])

            limb = self.compute_limb_length(matrix[x, j], matrix[y, j], value)
            if x != j and y != j and limb == 0:
                return x, y

    def additive_phylogeny(self, matrix: np.array) -> nx.Graph:
        graph, size = nx.Graph(), np.size(matrix, 0)

        if size == 2:
            graph.add_edge(0, 1, weight=matrix[0][1])
            return graph

        # 1. Pick an arbitrary leaf j (we choose THE LAST ONE)
        j = size - 1
        # 2. Compute its limb length, limb_length(matrix, j)
        limb = self.compute_min_limb_length(matrix, j)

        # 3. Subtract limb_length(matrix, j) from each row and column to
        # produce matrix_bald in which j is a bald limb (length 0)
        matrix[j, :], matrix[:, j], matrix[j, j] = matrix[j, :] - limb, matrix[:, j] - limb, 0

        # ?. (to 6) i, j, k leaves such that matrix_bald[i,j] + matrix_bald[j,k] = matrix_bald[i,k]
        i, k = self.find_attachment_point(matrix, j)
        matrix_i_j = matrix[i, j]

        # 4. Remove the j-th row and column of the matrix to form the (n-1) x (n-1) matrix_trim
        result = self.additive_phylogeny(matrix[:-1, :-1])

        # 5. Construct tree (matrix_trim)
        # 6. Identify the point in Tree (matrix_trim) where leaf j should be attached
        # 7. Attach j by an edge of length limb_length(matrix, j) in order to form Tree(matrix)
        new_node, path = self.new_nodes_start_index, nx.shortest_path(result, source=i, target=k)

        while new_node in result.nodes:
            new_node += 1

        distance, inner_distance, i = 0, 0, 0
        while distance < matrix_i_j:
            inner_distance = distance
            distance += result[path[i]][path[i + 1]]['weight']
            i += 1

        result.remove_edge(path[i - 1], path[i])
        result.add_edge(new_node, path[i], weight=distance - matrix_i_j)
        result.add_edge(new_node, path[i - 1], weight=matrix_i_j - inner_distance)
        result.add_edge(new_node, j, weight=limb)

        return result


if __name__ == "__main__":
    # files configuration
    input_dir, output_dir = 'examples/input/', 'examples/output/'
    files = ['n4_1.txt', 'n4_2.txt', 'n5_1.txt', 'n8_1.txt', 'n9_1_additive.txt']
    for f in files:
        input_path, file = input_dir + f, check_filename(f)

        # execute
        ap = AdditivePhylogeny()
        output = ap.additive_phylogeny(read_file(input_path))

        # show output
        weighted_adjacency_list(file, output)
        save_result(output, f'{output_dir}additive phylogeny/graph_{file}.png')
