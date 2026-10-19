def in_out_degree(graph):
    in_degree = {}
    out_degree = {}
    for vertex in graph.keys():
        in_degree[vertex] = 0
    for vertex, neighbours in graph.items():
        out_degree[vertex] = len(neighbours)
        for neighbour in neighbours:
            in_degree[neighbour] += 1
    return in_degree, out_degree


def dfs(graph, start):
    visited = set()
    order = []
    stack = [start]

    while stack:
        vertex = stack.pop()

        if vertex not in visited:
            visited.add(vertex)
            order.append(vertex)

            for neighbour in sorted(graph[vertex], reverse=True):
                if neighbour not in visited:
                    stack.append(neighbour)

    return order


def connected_components(graph):
    components = []
    visited = set()

    for vertex in graph:
        if vertex not in visited:
            component = set(dfs(graph, vertex))
            components.append(component)
            visited.update(component)
    return components


def is_weakly_connected(graph):
    shadow = {vertex: [] for vertex in graph}

    for vertex, neighbours in graph.items():
        for neighbour in neighbours:
            shadow[vertex].append(neighbour)
            shadow[neighbour].append(vertex)

    return len(connected_components(shadow)) == 1


def eulerian_type(graph, directed=False):
    # Directed graph
    if directed:
        in_degree, out_degree = in_out_degree(graph)

        # Check weak connectivity, ignoring isolated vertices
        shadow = {vertex: [] for vertex in graph}

        for vertex, neighbours in graph.items():
            for neighbour in neighbours:
                shadow[vertex].append(neighbour)
                shadow[neighbour].append(vertex)

        active = {
            vertex for vertex in graph
            if in_degree[vertex] + out_degree[vertex] > 0
        }

        for component in connected_components(shadow):
            if component & active and not active.issubset(component):
                return "none"

        # Calculate out-degree minus in-degree
        excess = {
            vertex: out_degree[vertex] - in_degree[vertex]
            for vertex in graph
        }

        if all(value == 0 for value in excess.values()):
            return "circuit"

        if (
            sum(value == 1 for value in excess.values()) == 1
            and sum(value == -1 for value in excess.values()) == 1
            and all(value in (-1, 0, 1) for value in excess.values())
        ):
            return "path"

        return "none"


def kmer_composition(text, k):
    kmers = []
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])
    return kmers


def debruijn_graph(kmers):
    debruijn = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix not in debruijn.keys():
            debruijn[prefix] = []
        if suffix not in debruijn.keys():
            debruijn[suffix] = []

        debruijn[prefix].append(suffix)
    return debruijn


def spell_path(nodes):
    path = nodes[0]
    for i in range(1, len(nodes)):
        path += nodes[i][-1]
    return path


if __name__ == "__main__":
    print(kmer_composition("GATTACA", 4))
    kmers = ["CAT", "ATC", "CAT"]
    graph = debruijn_graph(kmers)

    print(graph)
    print(spell_path(["GAT", "ATT", "TTA"]))
