import numpy as np
from scipy.sparse.csgraph import connected_components, floyd_warshall, dijkstra
import time

def matrix_sparsity(matrix: np.ndarray) -> float:
    return 1.0 - np.count_nonzero(matrix) / matrix.size

def matrix_closure(
    adjacency_matrix: np.ndarray,
    use_floyd_warshall: bool = False,
    use_dijkstra: bool = False
) -> np.ndarray:
    if use_floyd_warshall:
        return 1*np.isfinite(floyd_warshall(csgraph=np.ascontiguousarray(adjacency_matrix), directed=False))
    if use_dijkstra:
        return 1*np.isfinite(dijkstra(csgraph=np.ascontiguousarray(adjacency_matrix), directed=False))
    matrix_closure = adjacency_matrix
    while True:
        new_matrix_closure = matrix_closure | matrix_closure @ matrix_closure
        if np.array_equal(new_matrix_closure, matrix_closure):
            return matrix_closure
        matrix_closure = new_matrix_closure

def triu_closure(
    adjacency_triu: np.ndarray,
    use_floyd_warshall: bool = False,
    use_dijkstra: bool = False
) -> np.ndarray:
    return np.triu(
        matrix_closure(
            adjacency_triu + adjacency_triu.T,
            use_floyd_warshall=use_floyd_warshall, use_dijkstra=use_dijkstra),1)

def graph_connected_components(
    adjacency_matrix: np.ndarray, sparse_matrix: bool=False,
    use_floyd_warshall: bool=False, use_dijkstra: bool=False
) -> np.ndarray:
    if sparse_matrix:
        n_components, labels = connected_components(csgraph=np.ascontiguousarray(adjacency_matrix), directed=False)
        return np.array([labels == component_index for component_index in range(n_components)])
    adjacency_triu = np.triu(adjacency_matrix,1)
    adjacency_triu_closure = triu_closure(adjacency_triu, use_floyd_warshall=use_floyd_warshall, use_dijkstra=use_dijkstra)
    indexes_to_suppress = np.sum(adjacency_triu_closure, axis=0)
    np.fill_diagonal(adjacency_triu_closure, 1)
    return np.delete(adjacency_triu_closure, indexes_to_suppress > 0, axis=0)

def format_time_interval(secs: int) -> str:
    return f'{secs // 3600}:{(secs % 3600) // 60}:{secs % 60}'

def test_graph_connected_components(adjacency_matrix: np.ndarray) -> np.ndarray:
    start_time = time.process_time()
    result_prod = graph_connected_components(adjacency_matrix=adjacency_matrix, sparse_matrix=False, use_floyd_warshall=False)
    print(f'Time for prod: {format_time_interval(time.process_time() - start_time)}')

    start_time = time.process_time()
    result_connected = graph_connected_components(adjacency_matrix=adjacency_matrix, sparse_matrix=True)
    print(f'Time for connected_components: {format_time_interval(time.process_time() - start_time)}')
    print(f'Matrix with connected_components() equal: {np.array_equal(result_prod, result_connected)}')

    start_time = time.process_time()
    result_fw = graph_connected_components(adjacency_matrix=adjacency_matrix, use_floyd_warshall=True)
    print(f'Time for floyd_warshall: {format_time_interval(time.process_time() - start_time)}')
    print(f'Matrix with floyd_warshall() equal: {np.array_equal(result_prod, result_fw)}')

    start_time = time.process_time()
    result_di = graph_connected_components(adjacency_matrix=adjacency_matrix, use_dijkstra=True)
    print(f'Time for dijkstra: {format_time_interval(time.process_time() - start_time)}')
    print(f'Matrix with dijkstra() equal: {np.array_equal(result_prod, result_di)}')

    return result_prod


