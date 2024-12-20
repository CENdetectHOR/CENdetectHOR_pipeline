from collections.abc import Iterable
from matrix_closure import matrix_sparsity, graph_connected_components
import numpy as np
import editdistance
from Bio.SeqRecord import SeqRecord

def build_string_distance_matrix(strings: list[str]) -> np.ndarray:
    seq_triu = np.array([
            (0 if j <= i else editdistance.eval(seq_i, strings[j]))
            for i, seq_i in enumerate(strings)
            for j in range(len(strings))
    ])
    seq_triu.shape = (len(strings), len(strings))
    return seq_triu + seq_triu.T

def build_seqs_distance_matrix(seqs: list[SeqRecord]) -> np.ndarray:
    return build_string_distance_matrix([str(seq.seq) for seq in seqs])

def matrix_prod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    # return np.einsum('ij,jk->ikj', a, b)
    res = a[:,None] * b.T
    return res

def matrix_min(matrix, mask) -> np.ndarray:
    return np.array([np.min(matrix, where=cr, initial=np.max(matrix), axis=1) for cr in mask])
    # return (np.min(np.tile(matrix[:,:,None],len(mask)), where=mask.T, initial=np.max(matrix), axis=1)).T

def min_distance(distance_matrix: np.ndarray) -> any:
    return min(distance_matrix[np.triu_indices(distance_matrix.shape[0],1)])

def merge_clusters(
    distance_matrix: np.ndarray,
    clusters_expansion: np.ndarray = None,
    max_distance: any = None,
    sparsity_threshold: float = 0.97
) -> tuple[np.ndarray, np.ndarray, np.ndarray, any]:
    if max_distance is None:
        max_distance = min_distance(distance_matrix)
    print(f'merge_clusters with distance as {distance_matrix.shape}, clusters as {() if clusters_expansion is None else clusters_expansion.shape}, and max distance {max_distance}')
    adjacency_matrix = distance_matrix <= max_distance
    print(f'adjacency matrix sparsity is {matrix_sparsity(adjacency_matrix)}')

    new_clusters_expansion = graph_connected_components(adjacency_matrix, sparse_matrix=matrix_sparsity(adjacency_matrix) > sparsity_threshold)

    new_rows_distance_matrix = matrix_min(distance_matrix, new_clusters_expansion > 0)
    new_distance_matrix = matrix_min(new_rows_distance_matrix, new_clusters_expansion > 0)

    updated_clusters_expansion = new_clusters_expansion if clusters_expansion is None else new_clusters_expansion @ clusters_expansion
    return new_distance_matrix, updated_clusters_expansion, new_clusters_expansion, max_distance

def distance_values(matrix: np.ndarray):
    return matrix[np.triu_indices(matrix.shape[0], k = 1)]

def get_seq_as_txt(seq: list[int]) -> str:
    return "".join([(chr(65 + symbol) if symbol < 26 else (chr(71 + symbol) if symbol < 52 else '*')) for symbol in seq])

