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


def hierholzer(graph, start, verbose=False):
    unused = {vertex: sorted(neighbours) for vertex, neighbours in graph.items()}

    def make_walk(node):
        walk = [node]

        while unused[node]:
            neighbour = unused[node].pop(0)
            node = neighbour
            walk.append(node)

        return walk

    walk = make_walk(start)

    if verbose:
        print("Walk 1:", walk)

    walk_number = 1

    while any(unused.values()):
        restart_index = None

        for i, node in enumerate(walk):
            if unused[node]:
                restart_index = i
                break

        if restart_index is None:
            raise ValueError("Unused edges remain, but no restart node was found.")

        restart_node = walk[restart_index]
        new_walk = make_walk(restart_node)
        walk_number += 1

        if verbose:
            print(f"Walk {walk_number}, restarted at {restart_node}:", new_walk)

        walk = (
            walk[:restart_index]
            + new_walk
            + walk[restart_index + 1:]
        )

        if verbose:
            print("After splice:", walk)

    return walk


def eulerian_path(graph, verbose=False):
    if eulerian_type(graph, directed=True) == "none":
        raise ValueError("The graph has no Eulerian path.")

    in_degree, out_degree = in_out_degree(graph)

    start = min(graph)

    for vertex in graph:
        if out_degree[vertex] - in_degree[vertex] == 1:
            start = vertex
            break

    return hierholzer(graph, start, verbose)


def contigs(graph):
    in_degree, out_degree = in_out_degree(graph)

    unused = {
        vertex: sorted(neighbours)
        for vertex, neighbours in graph.items()
    }

    contig_paths = []

    starts = sorted(
        vertex for vertex in graph
        if in_degree[vertex] != 1 or out_degree[vertex] != 1
    )

    for start in starts:
        while unused[start]:
            path = [start]
            current = unused[start].pop(0)
            path.append(current)

            while in_degree[current] == 1 and out_degree[current] == 1:
                if not unused[current]:
                    break

                current = unused[current].pop(0)
                path.append(current)

            contig_paths.append(spell_path(path))

    while any(unused.values()):
        start = min(vertex for vertex in unused if unused[vertex])
        path = [start]
        current = start

        while True:
            current = unused[current].pop(0)
            path.append(current)

            if current == start:
                break

        contig_paths.append(spell_path(path))

    return sorted(contig_paths)


if __name__ == "__main__":
    text = "ACGTAC"
    kmers = kmer_composition(text, 3)
    graph = debruijn_graph(kmers)

    print("K-mers:", kmers)
    print("Graph:", graph)
    print("Contigs:", contigs(graph))
