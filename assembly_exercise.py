from Bio import SeqIO
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
    reads = [
        str(record.seq).upper()
        for record in SeqIO.parse("reads.fasta", "fasta")
    ]

    spelled_strings = {}

    for k in [9, 12, 13, 14, 21]:
        kmers = sorted({
            kmer
            for read in reads
            for kmer in kmer_composition(read, k)
        })

        graph = debruijn_graph(kmers)
        edge_count = sum(len(neighbours) for neighbours in graph.values())
        contig_list = contigs(graph)
        contig_count = len(contig_list)

        contig_edges = Counter()
        for contig in contig_list:
            contig_kmers = kmer_composition(contig, k)
            edges = Counter(
                zip(
                    [km[:-1] for km in contig_kmers],
                    [km[1:] for km in contig_kmers]
                )
            )

            for edge, count in edges.items():
                contig_edges[edge] += count

        graph_edges = Counter(
            (u, v)
            for u, neighbours in graph.items()
            for v in neighbours
        )

        print("\nContig edges correct:", contig_edges == graph_edges)

        for contig in contig_list:
            if not set(kmer_composition(contig, k)).issubset(set(kmers)):
                print("Warning: contig contains an unexpected k-mer")

        trimmed_contigs = [
            contig[:-(k - 1)]
            if contig[:k - 1] == contig[-(k - 1):]
            else contig
            for contig in contig_list
        ]

        longest_contig = max(map(len, trimmed_contigs), default=0)

        kind = eulerian_type(graph, directed=True)
        spelled = None
        trimmed_spelled = None

        if kind != "none":
            walk = eulerian_path(graph)
            spelled = spell_path(walk)

            print(
                "Spelled string length correct:",
                len(spelled) == edge_count + k - 1
            )
            print(
                "K-mer composition correct:",
                set(kmer_composition(spelled, k)) == set(kmers)
            )

            if spelled[:k - 1] == spelled[-(k - 1):]:
                trimmed_spelled = spelled[:-(k - 1)]
            else:
                trimmed_spelled = spelled

            spelled_strings[k] = trimmed_spelled
            record = SeqRecord(
                Seq(trimmed_spelled),
                id=f"assembly_k{k}",
                description=f"Assembled genome for k={k}"
            )

            SeqIO.write(record, f"blast_k{k}.fasta", "fasta")

        print(f"\nk = {k}")
        print("|V|:", len(graph))
        print("|E|:", edge_count)
        print("Eulerian type:", kind)
        print("Number of contigs:", contig_count)
        print("Longest contig (trimmed):", longest_contig)
        print(
            "Spelled string length (trimmed):",
            len(trimmed_spelled) if trimmed_spelled is not None else "-"
        )

        if trimmed_spelled is not None:
            print("Trimmed spelled string:", trimmed_spelled)

        print("Trimmed contigs:", trimmed_contigs)

    print("\nDo all spelled strings describe the same circle?")

    for k1, s in spelled_strings.items():
        for k2, t in spelled_strings.items():
            if k1 < k2:
                same_circle = len(s) == len(t) and s in t + t
                print(f"k={k1} vs k={k2}: {same_circle}")
