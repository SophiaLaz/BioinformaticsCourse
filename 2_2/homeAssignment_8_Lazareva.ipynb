from typing import List, Tuple


def build_graph(reads: List[str], threshold: int) -> List[Tuple[int, int, int]]:
    edges = []
    n, m = len(reads), len(reads[0])
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            max_cover = 0
            for l in range(1, m + 1):
                all_ok = (reads[i][m - l:] == reads[j][:l])
                if all_ok:
                    max_cover = l
            if max_cover >= threshold:
                edges.append((max_cover, i, j))
    return edges


def main(reads: List[str], threshold: int) -> str:
    vertexes = reads.copy()
    edges = build_graph(reads, threshold)
    for _ in range(len(reads) - 1):
        if len(edges) == 0:
            break
        t, i, j = max(edges)
        new_vertex = vertexes[i] + vertexes[j][t:]
        new_vertex_idx = len(vertexes)
        vertexes.append(new_vertex)
        vertexes[i] = ""
        vertexes[j] = ""
        new_edges = []
        for et, ei, ej in edges:
            if i == ei or j == ej:
                continue  # Drop edges to second vertex or from first vertex
            if i == ej and j == ei:
                continue  # Drop reverse edges
            if j == ei:
                ei = new_vertex_idx  # Make edge from new vertex
            if i == ej:
                ej = new_vertex_idx  # Make edge to new vertex
            new_edges.append((et, ei, ej))
        edges = new_edges
    result = "".join(vertexes)
    return result


if __name__ == '__main__':
    reads = ["AAA", "AAT", "ATT", "TTA", "TTT"]
    print(main(reads, 1))
