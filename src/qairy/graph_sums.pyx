from .misc import compute_rs_qairy_index, add_to_dict

from .partitions import (
    cached_increasing_red_partitions_with_min,
    tuple_automorphism_number,
    decreasing_partitions_with_min,
    partitions_with_min,
)

from sage.rings.all import Integer
from sage.graphs.all import Graph
from sage.matrix.all import Matrix
from sage.arith.all import factorial
from sage.misc.all import prod

from itertools import combinations


class QAiryGraph:
    """
    Quantum Airy Graph:
    Used to compute F for rs curves
    """

    def __init__(self, g, edges, tails):
        self.tails = tails
        self.g = g
        self.edges = edges

        # calculated properties
        self.num_V = len(self.g)
        self.num_tails = [len(tails) for tails in self.tails]
        self.val = self._val()
        self.automorphism_number = self._automorphism_number()

        # private variable for keeping track of assignment weights
        self._weight = Integer(0)

    def copy(self):
        new_tails = [t.copy() for t in self.tails]
        new_g = self.g.copy()
        new_edges = self.edges.copy()

        return QAiryGraph(new_g, new_edges, new_tails)

    # number of automorphisms fixing v_0
    def _automorphism_number(self):
        a = Integer(1)
        vert_types = []
        for v in range(1, self.num_V):
            n_t = self.num_tails[v]
            if n_t == 0:
                vert = (self.g[v], self.edges[v])
                vert_types.append(vert)
                a *= vert_types.count(vert)
        return a

    # valence needs to be calculated since tails can be added
    def _val(self):
        val = [self.edges[v] + self.num_tails[v] for v in range(self.num_V)]
        val[0] += sum(self.edges)
        return val

    def weight(self, qairy_coefficient, grading_condition, cache, r, s, b):
        """
        Calculate weight of graph, given tails b, by finding every possible
        edge assignment given:
        1. coefficients for the Airy structure C^(j)[q|a]
        2. grading condition for x's in F_gn[x_1,...x_n]
        3. cache of F_gn values
        """

        def _evaluate_edge_assignments(fat_vertex_assignment, v, b, w):
            # check if any vertices unassigned
            if v < self.num_V:
                g = self.g[v]
                n = self.val[v]
                tail_assignment = tuple(b[tail] for tail in self.tails[v])
                remaining = grading_condition(g, n) - sum(tail_assignment)

                # add every possible edge assignment from the fat vertex to v
                for edge_assignment in cached_increasing_red_partitions_with_min(
                    remaining, self.edges[v], r, 1
                ):
                    sorted_vertex_assignment = tuple(
                        sorted(edge_assignment + tail_assignment, reverse=True)
                    )

                    try:
                        val = cache[(g, n)][sorted_vertex_assignment]
                        _evaluate_edge_assignments(
                            fat_vertex_assignment + edge_assignment,
                            v + 1,
                            b,
                            w * val / tuple_automorphism_number(edge_assignment),
                        )
                    except KeyError:
                        pass
            else:
                w *= qairy_coefficient(b[0], self.g[0], fat_vertex_assignment)
                w *= prod(b[tail] for tail in self.tails[0] if tail > 0)

                self._weight += w

        # check if l + 2*j > i in which case C^(j)[q|a] = -1/r * D^(j)_i[k|a] = 0
        i, _ = compute_rs_qairy_index(r, s, b[0])
        j = self.g[0]
        l = len(self.tails[0]) - 1 + sum(self.edges[v] for v in range(1, self.num_V))

        if l + 2 * j > i:
            return Integer(0)

        self._weight = Integer(0)
        # edge assignments + tail assignments for fat vertex
        fat_vertex_assignment = tuple(-b[tail] for tail in self.tails[0] if tail > 0)
        _evaluate_edge_assignments(fat_vertex_assignment, 1, b, Integer(1))
        return self._weight

    def __str__(self):
        return "{}\n{}\n{}".format(self.g, self.edges, self.tails)

    def __repr__(self):
        return "{}\n{}\n{}".format(self.g, self.edges, self.tails)


# graph with no tails, zero genus, no edges
def zero_qa_graph(k):
    return QAiryGraph(
        [0 for i in range(k)],
        [0 if i == 0 else 0 for i in range(k)],
        [[0] if i == 0 else [] for i in range(k)],
    )


# generate quantum airy graphs for given (g, n)
def generate_qa_graphs(g, n, max_order):
    graph_cache = []
    # generate all graphs on k vertices
    def generate_k(g, n, k):
        # G = zero_qa_graph(k)
        G_data = [
            [0 for i in range(k)],  # genus
            [0 if i == 0 else 1 for i in range(k)],  # edges
            [[0] if i == 0 else [] for i in range(k)],  # tails
        ]
        assign_g(G_data, g, n, k, 0)

    # assign genus (j = G_data.g[0] should be at most floor(r/2) since l + 2*j <= r)
    def assign_g(G_data, g, n, k, i):
        if i < k:

            g_min = G_data[0][max(0, i - 1)]

            if i <= 1:
                g_min = 0

            g_max = g - sum(G_data[0])
            if i == 0:
                g_max = max_order // 2

            for g_v in range(g_min, g_max + 1):
                G_data[0][i] = g_v
                assign_g(G_data, g, n, k, i + 1)
                G_data[0][i] = 0

        else:
            assign_e(G_data, g, n, k, 1)

    # assign edges (l = total edges should be at most r since l + 2*j <= max_order)
    def assign_e(G_data, g, n, k, i):
        if i < k:
            if G_data[0][max(0, i - 1)] == G_data[0][i]:
                e_min = max(
                    1,
                    G_data[1][max(0, i - 1)],
                )
            else:
                e_min = 1

            e_max = min(
                g - (sum(G_data[0]) + sum(G_data[1]) - 1 - (k - 1)),
                max_order - 2 * G_data[0][0],
            )

            for e in range(e_min, e_max + 1):
                G_data[1][i] = e
                assign_e(G_data, g, n, k, i + 1)
                G_data[1][i] = 0

        else:
            if sum(G_data[0]) + sum(G_data[1]) - (k - 1) == g:
                num_tails = [0 for v in range(k)]
                assign_num_tails(G_data, g, n, k, 0, num_tails)

    # assign number of tails
    def assign_num_tails(G_data, g, n, k, i, num_tails):
        if i < k:
            if (
                G_data[1][max(0, i - 1)] == G_data[1][i]
                and G_data[0][max(0, i - 1)] == G_data[0][i]
            ):
                n_min = num_tails[max(0, i - 1)]
            else:
                n_min = 0

            if i == 0:
                n_min = 1
            elif i == 1:
                n_min = 0

            n_max = n - sum(num_tails)

            for n_v in range(n_min, n_max + 1):
                num_tails[i] = n_v
                assign_num_tails(G_data, g, n, k, i + 1, num_tails)
                num_tails[i] = 0

        else:
            if sum(num_tails) == n:

                for v in range(k):
                    valence = G_data[1][v] + num_tails[v]
                    if v == 0:
                        valence = sum(G_data[1]) + num_tails[v]
                    if G_data[0][v] == 0 and valence < 3:
                        return
                    if G_data[0][v] == 1 and valence < 1:
                        return

                # l + 2*j <= i <= r

                if (sum(G_data[1]) + num_tails[0] - 1) + 2 * G_data[0][0] > max_order:
                    return

                assign_tails(
                    G_data,
                    g,
                    n,
                    k,
                    0,
                    num_tails,
                    [i for i in range(1, n)],
                )

    # assign tails
    def assign_tails(G_data, g, n, k, i, num_tails, remaining):
        if i < k:
            n_t = num_tails[i]
            if i == 0:
                n_t -= 1

            for t in list(combinations(remaining, n_t)):
                tails = sorted(t)

                if (
                    i > 1
                    and n_t > 0
                    and num_tails[i - 1] == num_tails[i]
                    and G_data[1][i - 1] == G_data[1][i]
                    and G_data[0][i - 1] == G_data[0][i]
                ):
                    if tails[0] < G_data[2][i - 1][0]:
                        continue

                for tail in tails:
                    remaining.remove(tail)
                    G_data[2][i].append(tail)

                assign_tails(G_data, g, n, k, i + 1, num_tails, remaining)

                for tail in tails:
                    remaining.append(tail)
                    G_data[2][i].remove(tail)
        else:
            graph_cache.append(
                QAiryGraph(
                    G_data[0].copy(),
                    G_data[1].copy(),
                    [t.copy() for t in G_data[2]],
                )
            )

    # k <= r since there can be at most r edges and the
    # graph must be connected
    for k in range(
        1,
        min(2 * g - 2 + n, max_order + 1) + 1,
    ):
        generate_k(g, n, k)

    return graph_cache


class WGraph:
    """
    W Graph:
    Used to compute F for local spectral curves
    """

    __slots__ = [
        "genera",
        "sources",
        "loops",
        "sings",
        "dilatons",
        "adj",
        "num_V",
        "divisions",
        "_graph_data",
        "_val",
    ]
    # initialize source graph with FIXED number of vertices
    def __init__(self, genera, sources, loops, sings, dilatons, adj):
        self.genera = genera  # genus at each vertex
        self.sources = sources  # number of sources at each vertex
        self.loops = loops  # number of loops at each vertex
        self.sings = sings  # singularity associated to each vertex
        self.dilatons = dilatons
        self.adj = adj  # adjacency matrix
        self.num_V = len(genera)  # number of vertices

        self.divisions = [False for i in range(self.num_V)]
        self.divisions[0] = True

        self._graph_data = None
        self._val = None
        # self._automorphism_num = None

    # number of edges
    def num_E(self):
        return sum([sum(row) for row in self.adj]) // 2 + sum(self.loops)

    # total genus
    def g(self):
        return sum(self.genera) + self.num_E() - (self.num_V - 1)

    # total sources
    def num_sources(self):
        return sum(self.sources)

    # total dilatons
    def num_dilatons(self):
        return sum(self.dilatons)

    def val(self, v):
        if self._val != None:
            return self._val[v]
        else:
            self._val = [
                self.sources[v] + sum(self.adj[v]) + 2 * self.loops[v]
                for v in range(self.num_V)
            ]
            return self._val[v]

    # check if stable globally and at each vertex
    def is_stable(self):
        if 2 * self.g() - 2 + self.num_sources() <= 0:
            return False

        for i in range(self.num_V):
            if self.genera[i] == 0:
                n_v = self.sources[i] + 2 * self.loops[i] + sum(self.adj[i])
                if n_v < 3:
                    return False
        return True

    # copy all data
    def copy(self):
        new_genera = self.genera.copy()
        new_sources = self.sources.copy()
        new_loops = self.loops.copy()
        new_sings = self.sings.copy()
        new_dilatons = self.dilatons.copy()
        new_adj = [row.copy() for row in self.adj]
        new_divisions = self.divisions.copy()

        new_G = WGraph(
            new_genera, new_sources, new_loops, new_sings, new_dilatons, new_adj
        )
        new_G.divisions = new_divisions
        # new_G._graph_data = new_G._graph()
        return new_G

    # partition vertices by tuple (g,n,l,s)
    # order this vertex partition
    def get_vertex_partition(self):
        d = {}
        for i in range(self.num_V):
            dat = (
                self.genera[i],
                self.sources[i],
                self.loops[i],
                self.sings[i],
                self.dilatons[i],
            )
            if dat in d:
                d[dat].append(i)
            else:
                d[dat] = [i]
        return [d[dat] for dat in sorted(d)]

    # underlying graph, vertex data, partition of vertices by data
    def _graph(self):
        if self._graph_data != None:
            return self._graph_data
        else:
            graph = Graph(Matrix(self.adj), format="adjacency_matrix")
            vdat = list(
                zip(self.genera, self.sources, self.loops, self.sings, self.dilatons)
            )
            v_part = self.get_vertex_partition()

            self._graph_data = (graph, vdat, v_part)
            return self._graph_data

    # check if source graph is isomorphic to given
    def is_isomorphic(self, H):
        g, gvdat, gpart = self._graph()
        h, hvdat, hpart = H._graph()

        if sorted(gvdat) != sorted(hvdat):
            return False
        else:
            return g.canonical_label(
                partition=gpart, algorithm="sage"
            ) == h.canonical_label(partition=hpart, algorithm="sage")

    # automorphism group of the vertices
    def vertex_automorphism_group(self):
        g, _, gpart = self._graph()
        return g.automorphism_group(partition=gpart, algorithm="sage")

    # order of automorphism group including vertex permutations,
    # loop permutations, multi-edge permutations, source permutations
    def automorphism_number(self, include_tails=True):
        # if self._automorphism_num != None:
        # 	return self._automorphism_num
        # else:
        if include_tails:
            source_aut = prod([factorial(s) for s in self.sources])
            dilaton_aut = prod([factorial(d) for d in self.dilatons])
        else:
            source_aut = Integer(1)
            dilaton_aut = Integer(1)

        v_aut = self.vertex_automorphism_group().order()

        loop_aut = prod([2**l * factorial(l) for l in self.loops])

        edge_aut = Integer(1)
        for i in range(0, self.num_V):
            for j in range(i, self.num_V):
                num_edges = self.adj[i][j]
                if num_edges > 1:
                    edge_aut *= factorial(num_edges)

        # self._automorphism_num = v_aut * loop_aut * edge_aut * source_aut * dilaton_aut

        return v_aut * loop_aut * edge_aut * source_aut * dilaton_aut

    # check if connected
    def is_connected(self):
        visited = [0]
        unvisited = [i for i in range(1, self.num_V)]

        changed = True
        while changed and len(unvisited) > 0:
            changed = False
            for u in unvisited:
                for v in visited:
                    if self.adj[u][v] > 0:
                        if u not in visited:
                            visited.append(u)
                            unvisited.remove(u)
                            changed = True

        return len(visited) == self.num_V

    # string representation
    def __repr__(self):
        return str(
            [self.genera, self.sources, self.loops, self.sings, self.dilatons, self.adj]
        )

    def __str__(self):
        return str(
            [self.genera, self.sources, self.loops, self.sings, self.dilatons, self.adj]
        )


def zero_W_graph(k):
    return WGraph(
        [0 for i in range(k)],
        [0 for i in range(k)],
        [0 for i in range(k)],
        [0 for i in range(k)],
        [0 for i in range(k)],
        [[0 for i in range(k)] for i in range(k)],
    )


def generate_W_graphs(C, g, n, m):
    graph_cache = []

    for k in range(1, (2 * g - 2 + n) + 1):
        graph_cache.extend(generate_W_k_graphs(C, g, n, m, k))

    return graph_cache


def generate_W_k_graphs(C, g, n, m, k):
    graph_cache = []
    ram_graph_cache = []

    def assign_g(G, g, n, m, i):
        if i < G.num_V:
            g_min = G.genera[max(0, i - 1)]  # sequence is increasing
            g_max = g - sum(G.genera)  # sum does not exceed g

            for g_v in range(g_min, g_max + 1):
                G.genera[i] = g_v
                if g_v > g_min:
                    G.divisions[i] = True

                assign_g(G, g, n, m, i + 1)

                G.genera[i] = 0
                G.divisions[i] = False
        else:
            assign_n(G, g, n, m, 0)

    # then we assign n_v (sources at v) to each vertex
    def assign_n(G, g, n, m, i):
        if i < G.num_V:

            if G.divisions[i]:
                n_min = 0
            else:
                n_min = G.sources[i - 1]

            n_max = n - sum(G.sources)  # sum does not exceed n

            if i == G.num_V - 1:
                if n_max < n_min:
                    return
                else:
                    n_min = n_max

            for n_v in range(n_min, n_max + 1):
                G.sources[i] = n_v
                prev_div = G.divisions[i]
                if n_v > n_min:
                    G.divisions[i] = True

                assign_n(G, g, n, m, i + 1)

                G.sources[i] = 0
                G.divisions[i] = prev_div
        else:
            if G.num_V > 1 and G.sources[G.num_V - 1] > G.sources[G.num_V - 2]:
                G.divisions[G.num_V - 1] = True

            assign_l(G, g, n, m, 0)

    # then we assign l_v (loops at v) to each vertex
    def assign_l(G, g, n, m, i):
        if i < G.num_V:
            if G.divisions[i]:
                l_min = 0
            else:
                l_min = G.loops[i - 1]

            l_max = g - sum(G.genera) - sum(G.loops)  # sum does not exceed g

            for l_v in range(l_min, l_max + 1):
                G.loops[i] = l_v
                prev_div = G.divisions[i]
                if l_v > l_min:
                    G.divisions[i] = True
                assign_l(G, g, n, m, i + 1)
                G.loops[i] = 0
                G.divisions[i] = prev_div
        else:
            e_max = (G.num_V - 1) + g - sum(G.genera) - sum(G.loops)

            assign_a(G, g, n, m, 0, 1, e_max, [])

    # finally we create the adjacency matrix
    def assign_a(G, g, n, m, i, j, e_max, local_cache):
        if j >= G.num_V:
            # make sure v_i was connected
            if sum(G.adj[i]) == 0 and G.num_V > 1:
                return
            i += 1
            j = i + 1

        if i < G.num_V - 1:
            if G.divisions[i]:
                ai_min = 0
            else:
                ai_min = G.adj[i - 1][j]

            if j > i + 1 and not G.divisions[j]:
                aj_min = G.adj[i][j - 1]
            else:
                aj_min = 0

            a_min = max(aj_min, ai_min)
            a_max = e_max  # - number of disconnected vertices

            for a in range(a_min, a_max + 1):
                G.adj[i][j] = a
                G.adj[j][i] = a

                prev_div_i = G.divisions[i]
                prev_div_j = G.divisions[j]

                if a > ai_min:
                    G.divisions[i] = True

                if j > i + 1 and a > aj_min:
                    G.divisions[j] = True

                assign_a(G, g, n, m, i, j + 1, e_max - a, local_cache)

                G.adj[i][j] = 0
                G.adj[j][i] = 0

                G.divisions[i] = prev_div_i
                G.divisions[j] = prev_div_j
        else:
            # check if G is connected and stable
            if G.is_connected() and G.is_stable() and G.g() == g:
                # check if G is new
                new_graph = True
                G._graph_data = None
                for other in local_cache:
                    if G.is_isomorphic(other):
                        new_graph = False
                        break

                if new_graph:
                    local_cache.append(G.copy())
                    graph_cache.append(G.copy())

    # assign singularities; called after finding all
    # non-isomorphic graphs without singularity assignments
    def assign_m(G, g, n, m, i, local_cache):
        if i < G.num_V:
            if G.divisions[i]:
                m_min = 0
            else:
                m_min = G.sings[i - 1]

            # m_min = 0
            m_max = m - 1

            for m_v in range(m_min, m_max + 1):
                G.sings[i] = m_v
                prev_div = G.divisions[i]
                if m_v > m_min:
                    G.divisions[i] = True

                assign_m(G, g, n, m, i + 1, local_cache)

                G.sings[i] = 0
                G.divisions[i] = prev_div
        else:

            # check if G is new
            new_graph = True
            G._graph_data = None
            for other in local_cache:
                if G.is_isomorphic(other):
                    new_graph = False
                    break

            if new_graph:
                local_cache.append(G.copy())
                assign_d(G, g, n, m, 0, [])

    def assign_d(G, g, n, m, i, local_cache):
        if i < G.num_V:
            if G.divisions[i]:
                d_min = 0
            else:
                d_min = G.dilatons[i - 1]

            (r, s) = C.ramification_profile[G.sings[i]]
            g_v = G.genera[i]
            n_v = G.val(i)
            d_max = s * (2 * g_v - 2 + n_v) - n_v

            if r == 2:
                d_max //= 2

            for d_v in range(d_min, d_max + 1):
                G.dilatons[i] = d_v
                prev_div = G.divisions[i]
                if d_v > d_min:
                    G.divisions[i] = True

                assign_d(G, g, n, m, i + 1, local_cache)

                G.dilatons[i] = 0
                G.divisions[i] = prev_div
        else:
            # check if G is new
            new_graph = True
            G._graph_data = None
            for other in local_cache:
                if G.is_isomorphic(other):
                    new_graph = False
                    break

            if new_graph:
                local_cache.append(G.copy())
                ram_graph_cache.append(G.copy())

    # generate graphs
    G = zero_W_graph(k)
    assign_g(G, g, n, m, 0)

    for G in graph_cache:
        G._graph_data = None
        assign_m(G, g, n, m, 0, [])

    return ram_graph_cache


def get_poly_dict(C, G_list, rs_cache):
    poly_dict = {}

    def generate_assignments(G):
        assignments = {}

        # for each vertex add all possible assignments
        for v in range(0, G.num_V):
            assignments[v] = []

            g = G.genera[v]

            # half edges
            n = G.sources[v]
            d = G.dilatons[v]
            l = 2 * G.loops[v]
            e = sum(G.adj[v])

            val = n + d + l + e

            (r, s) = C.ramification_profile[G.sings[v]]

            min_d = s + 1
            if min_d % r == 0:
                min_d += 1

            # grading condition
            total = s * (2 * g - 2 + val)

            # assign 1 to n, l, e and min_d to d
            extra = total - n - min_d * d - l - e

            n_max = extra if n > 0 else 0

            # assign source labels
            for n_total in range(n_max + 1):
                n_assignments = (
                    decreasing_partitions_with_min(n_total, n) if n > 0 else [tuple()]
                )

                for n_assignment in n_assignments:
                    d_max = extra - n_total if d > 0 else 0

                    # assign dilaton labels
                    for d_total in range(d_max + 1):
                        d_assignments = (
                            decreasing_partitions_with_min(d_total, d)
                            if d > 0
                            else [tuple()]
                        )

                        for d_assignment in d_assignments:
                            l_max = extra - n_total - d_total if l > 0 else 0

                            # assign loop labels
                            for l_total in range(l_max + 1):
                                l_assignments = (
                                    partitions_with_min(l_total, l)
                                    if l > 0
                                    else [tuple()]
                                )

                                for l_assignment in l_assignments:
                                    e_total = (
                                        extra - n_total - d_total - l_total
                                        if e > 0
                                        else 0
                                    )

                                    e_assignments = (
                                        partitions_with_min(e_total, e)
                                        if e > 0
                                        else [tuple()]
                                    )

                                    # assign edge labels
                                    for e_assignment in e_assignments:
                                        extra_assignment = (
                                            n_assignment
                                            + d_assignment
                                            + l_assignment
                                            + e_assignment
                                        )

                                        if sum(extra_assignment) == extra:
                                            assignment = tuple(
                                                extra_assignment[i] + min_d
                                                if i >= n and i < n + d
                                                else extra_assignment[i] + 1
                                                for i in range(val)
                                            )
                                            if all(
                                                label % r != 0 for label in assignment
                                            ):
                                                assignments[v].append(assignment)

            # if there are no possibilities return
            if len(assignments[v]) == 0:
                return

        choose_assignment(G, assignments, {v: [] for v in range(G.num_V)}, 0)

    def choose_assignment(G, assignments, assignment, v):
        if v < G.num_V:
            for v_assignment in assignments[v]:
                prev_v_assignment = assignment[v].copy()

                # set assignment for v
                assignment[v] = list(v_assignment)

                # proceed to next vertex
                choose_assignment(G, assignments, assignment, v + 1)

                # backtrack
                assignment[v] = prev_v_assignment
        else:
            evaluate_assignment(G, assignment)

    # half edges are assumed to be in the following order:
    # [sources, dilaton leaves, loops, edges]
    # adds weight of term to poly_dict
    def evaluate_assignment(G, A):
        w = Integer(1)

        e_offset = {
            v: G.sources[v] + G.dilatons[v] + 2 * G.loops[v] for v in range(G.num_V)
        }

        term_dict = {var_name: [] for var_name in ["F02", "F01", "x"]}

        for v in range(G.num_V):
            g = G.genera[v]

            # half edges
            n = G.sources[v]
            d = G.dilatons[v]
            l = 2 * G.loops[v]
            e = sum(G.adj[v])

            val = n + d + l + e
            chi = 2 * g - 2 + val

            sing = G.sings[v]
            (r, s) = C.ramification_profile[sing]

            assignment = A[v]
            sorted_b = tuple(sorted(assignment, reverse=True))

            # get F[g, val] from cache
            if (r, s) in rs_cache:
                if (g, val) in rs_cache[(r, s)]:
                    if sorted_b in rs_cache[(r, s)][(g, val)]:
                        w *= rs_cache[(r, s)][(g, val)][sorted_b]
                    else:
                        # if weight for any vertex is zero return zero
                        return Integer(0)
            else:
                F_rs_str = "F_rs[{},{}][{},{}]".format(r, s, g, val)
                if F_rs_str not in term_dict:
                    term_dict[F_rs_str] = []

                term_dict[F_rs_str].append((sorted_b, 1))
                # raise Exception("F[{},{}] not catch for ({},{}).".format(g, val, r, s))

            n_assignment = assignment[:n]
            d_assignment = assignment[n : n + d]
            l_assignment = assignment[n + d : n + d + l]

            # sources
            for n_label in n_assignment:
                term_dict["x"].append(((sing + 1, n_label), 1))

            w /= tuple_automorphism_number(n_assignment)

            # dilaton leaves weighted by F01[>s]
            for d_label in d_assignment:
                term_dict["F01"].append(((sing + 1, d_label), 1))
                w /= d_label

            w /= tuple_automorphism_number(d_assignment)

            # loops weighted by Phi
            for i in range(l // 2):
                (l1, l2) = l_assignment[2 * i : 2 * i + 2]

                term_dict["F02"].append(
                    (
                        (sing + 1, sing + 1, l1, l2),
                        1,
                    )
                )

                w /= l1 * l2

            # edges weighted by Phi
            for other in range(v + 1, G.num_V):
                for e in range(G.adj[v][other]):
                    e1 = assignment[e_offset[v]]
                    e2 = A[other][e_offset[other]]

                    term_dict["F02"].append(
                        (
                            (
                                sing + 1,
                                G.sings[other] + 1,
                                e1,
                                e2,
                            ),
                            1,
                        )
                    )

                    e_offset[v] += 1
                    e_offset[other] += 1

                    w /= e1 * e2

            # dilaton leaves weighted by F01[s]
            term_dict["F01"].append(((sing + 1, s), -chi))

            w /= (-1) ** (chi)

        # automorphism factor
        w /= G.automorphism_number(include_tails=False)

        term = flatten_term_dict(term_dict)

        add_to_dict(poly_dict, term, w)

    def flatten_term_dict(term_dict):
        # collect powers
        powers = {}
        for var_name in term_dict:
            powers[var_name] = {}

            if var_name == "F02":
                for (a1, a2, l1, l2), power in term_dict[var_name]:

                    if l1 > l2 or (l1 == l2 and a1 > a2):
                        indices = (a2, a1, l2, l1)
                    else:
                        indices = (a1, a2, l1, l2)

                    add_to_dict(powers[var_name], indices, power)

            else:
                for (indices, power) in term_dict[var_name]:
                    add_to_dict(powers[var_name], indices, power)

        term = []
        for var_name in powers:
            term.extend(
                [
                    (var_name, indices, powers[var_name][indices])
                    for indices in powers[var_name]
                ]
            )

        return tuple(sorted(term))

    for G in G_list:
        generate_assignments(G)

    return poly_dict


def get_poly_dict_list(C, G_list, rs_cache):
    poly_dict_list = []

    def generate_assignments(G):
        assignments = {}

        # for each vertex add all possible assignments
        for v in range(0, G.num_V):
            assignments[v] = []

            g = G.genera[v]

            # half edges
            n = G.sources[v]
            d = G.dilatons[v]
            l = 2 * G.loops[v]
            e = sum(G.adj[v])

            val = n + d + l + e

            (r, s) = C.ramification_profile[G.sings[v]]

            min_d = s + 1
            if min_d % r == 0:
                min_d += 1

            # grading condition
            total = s * (2 * g - 2 + val)

            # assign 1 to n, l, e and min_d to d
            extra = total - n - min_d * d - l - e

            n_max = extra if n > 0 else 0

            # assign source labels
            for n_total in range(n_max + 1):
                n_assignments = (
                    decreasing_partitions_with_min(n_total, n) if n > 0 else [tuple()]
                )

                for n_assignment in n_assignments:
                    d_max = extra - n_total if d > 0 else 0

                    # assign dilaton labels
                    for d_total in range(d_max + 1):
                        d_assignments = (
                            decreasing_partitions_with_min(d_total, d)
                            if d > 0
                            else [tuple()]
                        )

                        for d_assignment in d_assignments:
                            l_max = extra - n_total - d_total if l > 0 else 0

                            # assign loop labels
                            for l_total in range(l_max + 1):
                                l_assignments = (
                                    partitions_with_min(l_total, l)
                                    if l > 0
                                    else [tuple()]
                                )

                                for l_assignment in l_assignments:
                                    e_total = (
                                        extra - n_total - d_total - l_total
                                        if e > 0
                                        else 0
                                    )

                                    e_assignments = (
                                        partitions_with_min(e_total, e)
                                        if e > 0
                                        else [tuple()]
                                    )

                                    # assign edge labels
                                    for e_assignment in e_assignments:
                                        extra_assignment = (
                                            n_assignment
                                            + d_assignment
                                            + l_assignment
                                            + e_assignment
                                        )

                                        if sum(extra_assignment) == extra:
                                            assignment = tuple(
                                                extra_assignment[i] + min_d
                                                if i >= n and i < n + d
                                                else extra_assignment[i] + 1
                                                for i in range(val)
                                            )
                                            if all(
                                                label % r != 0 for label in assignment
                                            ):
                                                assignments[v].append(assignment)

            # if there are no possibilities return
            if len(assignments[v]) == 0:
                return

        choose_assignment(G, assignments, {v: [] for v in range(G.num_V)}, 0)

    def choose_assignment(G, assignments, assignment, v):
        if v < G.num_V:
            for v_assignment in assignments[v]:
                prev_v_assignment = assignment[v].copy()

                # set assignment for v
                assignment[v] = list(v_assignment)

                # proceed to next vertex
                choose_assignment(G, assignments, assignment, v + 1)

                # backtrack
                assignment[v] = prev_v_assignment
        else:
            evaluate_assignment(G, assignment)

    # half edges are assumed to be in the following order:
    # [sources, dilaton leaves, loops, edges]
    # adds weight of term to poly_dict
    def evaluate_assignment(G, A):
        w = Integer(1)
        e_offset = {
            v: G.sources[v] + G.dilatons[v] + 2 * G.loops[v] for v in range(G.num_V)
        }

        term_dict = {var_name: [] for var_name in ["F02", "F01", "x"]}

        for v in range(G.num_V):
            g = G.genera[v]

            # half edges
            n = G.sources[v]
            d = G.dilatons[v]
            l = 2 * G.loops[v]
            e = sum(G.adj[v])

            val = n + d + l + e
            chi = 2 * g - 2 + val

            sing = G.sings[v]
            (r, s) = C.ramification_profile[sing]

            assignment = A[v]
            sorted_b = tuple(sorted(assignment, reverse=True))

            # get F[g, val] from cache
            if (r, s) in rs_cache:
                if (g, val) in rs_cache[(r, s)]:
                    if sorted_b in rs_cache[(r, s)][(g, val)]:
                        w *= rs_cache[(r, s)][(g, val)][sorted_b]
                    else:
                        # if weight for any vertex is zero return zero
                        return Integer(0)
            else:
                F_rs_str = "F_rs[{},{}][{},{}]".format(r, s, g, val)
                if F_rs_str not in term_dict:
                    term_dict[F_rs_str] = []

                term_dict[F_rs_str].append((sorted_b, 1))
                # raise Exception("F[{},{}] not catch for ({},{}).".format(g, val, r, s))

            n_assignment = assignment[:n]
            d_assignment = assignment[n : n + d]
            l_assignment = assignment[n + d : n + d + l]

            # sources
            for n_label in n_assignment:
                term_dict["x"].append(((sing + 1, n_label), 1))

            w /= tuple_automorphism_number(n_assignment)

            # dilaton leaves weighted by F01[>s]
            for d_label in d_assignment:
                term_dict["F01"].append(((sing + 1, d_label), 1))
                w /= d_label

            w /= tuple_automorphism_number(d_assignment)

            # loops weighted by Phi
            for i in range(l // 2):
                (l1, l2) = l_assignment[2 * i : 2 * i + 2]

                term_dict["F02"].append(
                    (
                        (sing + 1, sing + 1, l1, l2),
                        1,
                    )
                )

                w /= l1 * l2

            # edges weighted by Phi
            for other in range(v + 1, G.num_V):
                for e in range(G.adj[v][other]):
                    e1 = assignment[e_offset[v]]
                    e2 = A[other][e_offset[other]]

                    term_dict["F02"].append(
                        (
                            (
                                sing + 1,
                                G.sings[other] + 1,
                                e1,
                                e2,
                            ),
                            1,
                        )
                    )

                    e_offset[v] += 1
                    e_offset[other] += 1

                    w /= e1 * e2

            # dilaton leaves weighted by F01[s]
            term_dict["F01"].append(((sing + 1, s), -chi))

            w /= (-1) ** (chi)

        # automorphism factor
        w /= G.automorphism_number(include_tails=False)

        term = flatten_term_dict(term_dict)

        poly_dict_list.append((term, w))

    def flatten_term_dict(term_dict):
        # collect powers
        powers = {}
        for var_name in term_dict:
            powers[var_name] = {}

            if var_name == "F02":
                for (a1, a2, l1, l2), power in term_dict[var_name]:

                    if l1 > l2 or (l1 == l2 and a1 > a2):
                        indices = (a2, a1, l2, l1)
                    else:
                        indices = (a1, a2, l1, l2)

                    add_to_dict(powers[var_name], indices, power)

            else:
                for (indices, power) in term_dict[var_name]:
                    add_to_dict(powers[var_name], indices, power)

        term = []
        for var_name in powers:
            term.extend(
                [
                    (var_name, indices, powers[var_name][indices])
                    for indices in powers[var_name]
                ]
            )

        return tuple(sorted(term))

    for G in G_list:
        generate_assignments(G)

    return poly_dict_list


def add_contribution_to_poly_dict(poly_dict, C, G, rs_cache):
    def generate_assignments(G):
        assignments = {}

        # for each vertex add all possible assignments
        for v in range(0, G.num_V):
            assignments[v] = []

            g = G.genera[v]

            # half edges
            n = G.sources[v]
            d = G.dilatons[v]
            l = 2 * G.loops[v]
            e = sum(G.adj[v])

            val = n + d + l + e

            (r, s) = C.ramification_profile[G.sings[v]]

            min_d = s + 1
            if min_d % r == 0:
                min_d += 1

            # grading condition
            total = s * (2 * g - 2 + val)

            # assign 1 to n, l, e and min_d to d
            extra = total - n - min_d * d - l - e

            n_max = extra if n > 0 else 0

            # assign source labels
            for n_total in range(n_max + 1):
                n_assignments = (
                    decreasing_partitions_with_min(n_total, n) if n > 0 else [tuple()]
                )

                for n_assignment in n_assignments:
                    d_max = extra - n_total if d > 0 else 0

                    # assign dilaton labels
                    for d_total in range(d_max + 1):
                        d_assignments = (
                            decreasing_partitions_with_min(d_total, d)
                            if d > 0
                            else [tuple()]
                        )

                        for d_assignment in d_assignments:
                            l_max = extra - n_total - d_total if l > 0 else 0

                            # assign loop labels
                            for l_total in range(l_max + 1):
                                l_assignments = (
                                    partitions_with_min(l_total, l)
                                    if l > 0
                                    else [tuple()]
                                )

                                for l_assignment in l_assignments:
                                    e_total = (
                                        extra - n_total - d_total - l_total
                                        if e > 0
                                        else 0
                                    )

                                    e_assignments = (
                                        partitions_with_min(e_total, e)
                                        if e > 0
                                        else [tuple()]
                                    )

                                    # assign edge labels
                                    for e_assignment in e_assignments:
                                        extra_assignment = (
                                            n_assignment
                                            + d_assignment
                                            + l_assignment
                                            + e_assignment
                                        )

                                        if sum(extra_assignment) == extra:
                                            assignment = tuple(
                                                extra_assignment[i] + min_d
                                                if i >= n and i < n + d
                                                else extra_assignment[i] + 1
                                                for i in range(val)
                                            )
                                            if all(
                                                label % r != 0 for label in assignment
                                            ):
                                                assignments[v].append(assignment)

            # if there are no possibilities return
            if len(assignments[v]) == 0:
                return

        choose_assignment(G, assignments, {v: [] for v in range(G.num_V)}, 0)

    def choose_assignment(G, assignments, assignment, v):
        if v < G.num_V:
            for v_assignment in assignments[v]:
                prev_v_assignment = assignment[v].copy()

                # set assignment for v
                assignment[v] = list(v_assignment)

                # proceed to next vertex
                choose_assignment(G, assignments, assignment, v + 1)

                # backtrack
                assignment[v] = prev_v_assignment
        else:
            evaluate_assignment(G, assignment)

    # half edges are assumed to be in the following order:
    # [sources, dilaton leaves, loops, edges]
    # adds weight of term to poly_dict
    def evaluate_assignment(G, A):
        w = Integer(1)
        e_offset = {
            v: G.sources[v] + G.dilatons[v] + 2 * G.loops[v] for v in range(G.num_V)
        }

        term_dict = {var_name: [] for var_name in ["F02", "F01", "x"]}

        for v in range(G.num_V):
            g = G.genera[v]

            # half edges
            n = G.sources[v]
            d = G.dilatons[v]
            l = 2 * G.loops[v]
            e = sum(G.adj[v])

            val = n + d + l + e
            chi = 2 * g - 2 + val

            sing = G.sings[v]
            (r, s) = C.ramification_profile[sing]

            assignment = A[v]
            sorted_b = tuple(sorted(assignment, reverse=True))

            # get F[g, val] from cache
            if (r, s) in rs_cache:
                if (g, val) in rs_cache[(r, s)]:
                    if sorted_b in rs_cache[(r, s)][(g, val)]:
                        w *= rs_cache[(r, s)][(g, val)][sorted_b]
                    else:
                        # if weight for any vertex is zero return zero
                        return Integer(0)
            else:
                F_rs_str = "F_rs[{},{}][{},{}]".format(r, s, g, val)
                if F_rs_str not in term_dict:
                    term_dict[F_rs_str] = []

                term_dict[F_rs_str].append((sorted_b, 1))
                # raise Exception("F[{},{}] not catch for ({},{}).".format(g, val, r, s))

            n_assignment = assignment[:n]
            d_assignment = assignment[n : n + d]
            l_assignment = assignment[n + d : n + d + l]

            # sources
            for n_label in n_assignment:
                term_dict["x"].append(((sing + 1, n_label), 1))

            w /= tuple_automorphism_number(n_assignment)

            # dilaton leaves weighted by F01[>s]
            for d_label in d_assignment:
                term_dict["F01"].append(((sing + 1, d_label), 1))
                w /= d_label

            w /= tuple_automorphism_number(d_assignment)

            # loops weighted by Phi
            for i in range(l // 2):
                (l1, l2) = l_assignment[2 * i : 2 * i + 2]

                term_dict["F02"].append(
                    (
                        (sing + 1, sing + 1, l1, l2),
                        1,
                    )
                )

                w /= l1 * l2

            # edges weighted by Phi
            for other in range(v + 1, G.num_V):
                for e in range(G.adj[v][other]):
                    e1 = assignment[e_offset[v]]
                    e2 = A[other][e_offset[other]]

                    term_dict["F02"].append(
                        (
                            (
                                sing + 1,
                                G.sings[other] + 1,
                                e1,
                                e2,
                            ),
                            1,
                        )
                    )

                    e_offset[v] += 1
                    e_offset[other] += 1

                    w /= e1 * e2

            # dilaton leaves weighted by F01[s]
            term_dict["F01"].append(((sing + 1, s), -chi))

            w /= (-1) ** (chi)

        # automorphism factor
        w /= G.automorphism_number(include_tails=False)

        term = flatten_term_dict(term_dict)

        add_to_dict(poly_dict, term, w)

    def flatten_term_dict(term_dict):
        # collect powers
        powers = {}
        for var_name in term_dict:
            powers[var_name] = {}

            if var_name == "F02":
                for (a1, a2, l1, l2), power in term_dict[var_name]:

                    if l1 > l2 or (l1 == l2 and a1 > a2):
                        indices = (a2, a1, l2, l1)
                    else:
                        indices = (a1, a2, l1, l2)

                    add_to_dict(powers[var_name], indices, power)

            else:
                for (indices, power) in term_dict[var_name]:
                    add_to_dict(powers[var_name], indices, power)

        term = []
        for var_name in powers:
            term.extend(
                [
                    (var_name, indices, powers[var_name][indices])
                    for indices in powers[var_name]
                ]
            )

        return tuple(sorted(term))

    generate_assignments(G)


def add_contribution_of_W_k_graphs(C, g, n, m, k):
    graph_cache = []

    def assign_g(G, g, n, m, i):
        if i < G.num_V:
            g_min = G.genera[max(0, i - 1)]  # sequence is increasing
            g_max = g - sum(G.genera)  # sum does not exceed g

            for g_v in range(g_min, g_max + 1):
                G.genera[i] = g_v
                if g_v > g_min:
                    G.divisions[i] = True

                assign_g(G, g, n, m, i + 1)

                G.genera[i] = 0
                G.divisions[i] = False
        else:
            assign_n(G, g, n, m, 0)

    # then we assign n_v (sources at v) to each vertex
    def assign_n(G, g, n, m, i):
        if i < G.num_V:

            if G.divisions[i]:
                n_min = 0
            else:
                n_min = G.sources[i - 1]

            n_max = n - sum(G.sources)  # sum does not exceed n

            if i == G.num_V - 1:
                if n_max < n_min:
                    return
                else:
                    n_min = n_max

            for n_v in range(n_min, n_max + 1):
                G.sources[i] = n_v
                prev_div = G.divisions[i]
                if n_v > n_min:
                    G.divisions[i] = True

                assign_n(G, g, n, m, i + 1)

                G.sources[i] = 0
                G.divisions[i] = prev_div
        else:
            if G.num_V > 1 and G.sources[G.num_V - 1] > G.sources[G.num_V - 2]:
                G.divisions[G.num_V - 1] = True

            assign_l(G, g, n, m, 0)

    # then we assign l_v (loops at v) to each vertex
    def assign_l(G, g, n, m, i):
        if i < G.num_V:
            if G.divisions[i]:
                l_min = 0
            else:
                l_min = G.loops[i - 1]

            l_max = g - sum(G.genera) - sum(G.loops)  # sum does not exceed g

            for l_v in range(l_min, l_max + 1):
                G.loops[i] = l_v
                prev_div = G.divisions[i]
                if l_v > l_min:
                    G.divisions[i] = True
                assign_l(G, g, n, m, i + 1)
                G.loops[i] = 0
                G.divisions[i] = prev_div
        else:
            e_max = (G.num_V - 1) + g - sum(G.genera) - sum(G.loops)

            assign_a(G, g, n, m, 0, 1, e_max, [])

    # finally we create the adjacency matrix
    def assign_a(G, g, n, m, i, j, e_max, local_cache):
        if j >= G.num_V:
            # make sure v_i was connected
            if sum(G.adj[i]) == 0 and G.num_V > 1:
                return
            i += 1
            j = i + 1

        if i < G.num_V - 1:
            if G.divisions[i]:
                ai_min = 0
            else:
                ai_min = G.adj[i - 1][j]

            if j > i + 1 and not G.divisions[j]:
                aj_min = G.adj[i][j - 1]
            else:
                aj_min = 0

            a_min = max(aj_min, ai_min)
            a_max = e_max  # - number of disconnected vertices

            for a in range(a_min, a_max + 1):
                G.adj[i][j] = a
                G.adj[j][i] = a

                prev_div_i = G.divisions[i]
                prev_div_j = G.divisions[j]

                if a > ai_min:
                    G.divisions[i] = True

                if j > i + 1 and a > aj_min:
                    G.divisions[j] = True

                assign_a(G, g, n, m, i, j + 1, e_max - a, local_cache)

                G.adj[i][j] = 0
                G.adj[j][i] = 0

                G.divisions[i] = prev_div_i
                G.divisions[j] = prev_div_j
        else:
            # check if G is connected and stable
            if G.is_connected() and G.is_stable() and G.g() == g:
                # check if G is new
                new_graph = True
                G._graph_data = None
                for other in local_cache:
                    if G.is_isomorphic(other):
                        new_graph = False
                        break

                if new_graph:
                    local_cache.append(G.copy())
                    graph_cache.append(G.copy())

    # assign singularities; called after finding all
    # non-isomorphic graphs without singularity assignments
    def assign_m(G, g, n, m, i, local_cache):
        if i < G.num_V:
            if G.divisions[i]:
                m_min = 0
            else:
                m_min = G.sings[i - 1]

            # m_min = 0
            m_max = m - 1

            for m_v in range(m_min, m_max + 1):
                G.sings[i] = m_v
                prev_div = G.divisions[i]
                if m_v > m_min:
                    G.divisions[i] = True

                assign_m(G, g, n, m, i + 1, local_cache)

                G.sings[i] = 0
                G.divisions[i] = prev_div
        else:

            # check if G is new
            new_graph = True
            G._graph_data = None
            for other in local_cache:
                if G.is_isomorphic(other):
                    new_graph = False
                    break

            if new_graph:
                local_cache.append(G.copy())
                assign_d(G, g, n, m, 0, [])

    def assign_d(G, g, n, m, i, local_cache):
        if i < G.num_V:
            if G.divisions[i]:
                d_min = 0
            else:
                d_min = G.dilatons[i - 1]

            (r, s) = C.ramification_profile[G.sings[i]]
            g_v = G.genera[i]
            n_v = G.val(i)
            d_max = s * (2 * g_v - 2 + n_v) - n_v

            if r == 2:
                d_max //= 2

            for d_v in range(d_min, d_max + 1):
                G.dilatons[i] = d_v
                prev_div = G.divisions[i]
                if d_v > d_min:
                    G.divisions[i] = True

                assign_d(G, g, n, m, i + 1, local_cache)

                G.dilatons[i] = 0
                G.divisions[i] = prev_div
        else:
            # check if G is new
            new_graph = True
            G._graph_data = None
            for other in local_cache:
                if G.is_isomorphic(other):
                    new_graph = False
                    break

            if new_graph:
                local_cache.append(G.copy())

                add_contribution_to_poly_dict(C.F_dict[(g, n)], C, G, C._rs_cache)

    # generate graphs
    G = zero_W_graph(k)
    assign_g(G, g, n, m, 0)

    for G in graph_cache:
        G._graph_data = None
        assign_m(G, g, n, m, 0, [])
