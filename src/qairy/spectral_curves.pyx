from .graph_sums import (
    generate_W_k_graphs,
    generate_qa_graphs,
    add_contribution_of_W_k_graphs,
    get_poly_dict,
    add_contribution_to_poly_dict,
    WGraph,
)

from .generating_function import GeneratingFunction, load_F_for_rs_curve

from .partitions import (
    cached_decreasing_red_partitions_with_min,
    decreasing_red_partitions_with_min,
    tuple_automorphism_number,
    partitions_with_min,
    get_pairs,
)

from .misc import add_to_dict, compute_rs_qairy_index

from sage.structure.sage_object import SageObject
from sage.arith.all import factorial, gcd
from sage.misc.all import prod, cached_function, cached_method
from sage.rings.all import Integer, E, QQ

from csv import reader as csv_reader
from csv import writer as csv_writer

from itertools import permutations


class LocalSpectralCurve(SageObject):
    """
    Local spectral curve with ramification profile [(r_1,s_1), ... , (r_d,s_d)]
    """

    def __init__(self, ramification_profile, name=None):
        assert is_admissible_ramification_profile(ramification_profile)

        self.ramification_profile = ramification_profile

        if name == None:
            name = "rams=" + str(tuple(self.ramification_profile)).replace(
                " ", ""
            ).replace(",)", ")")

        self.name = name
        self._dict = {}

    def invariants(self, chi, use_rs_values=True):
        return _convert_F_to_w(self.compute_F(chi, use_rs_values))

    def generating_function(self, chi, use_rs_values=True):
        return self.compute_F(chi, use_rs_values)

    def compute_F(self, chi, use_rs_values=True):
        if use_rs_values:
            self._load_caches(chi)

        for x in range(1, chi + 1):
            for g in range((x + 3) // 2):
                n = x + 2 - 2 * g
                self.compute_F_gn(g, n)

        return GeneratingFunction(self._dict, "F_" + self.name)

    def compute_F_gn(self, g, n, use_rs_values=True):
        if use_rs_values:
            self._load_caches(2 * g - 2 + n)

        self._dict[(g, n)] = {}

        for k in range(1, 2 * g - 2 + n + 1):
            G_list = generate_W_k_graphs(self, g, n, len(self.ramification_profile), k)

            for term, w in get_poly_dict(self, G_list, self._rs_cache).items():
                add_to_dict(self._dict[(g, n)], term, w)

        return GeneratingFunction(self._dict, "F_" + self.name)

    def _compute_F_gn_without_saving_graphs(self, g, n, use_rs_values=True):
        if use_rs_values:
            self._load_caches(2 * g - 2 + n)

        self._dict[(g, n)] = {}

        for k in range(1, 2 * g - 2 + n + 1):
            add_contribution_of_W_k_graphs(
                self, g, n, len(self.ramification_profile), k
            )

        return GeneratingFunction(self._dict, "F_" + self.name)

    def _compute_F_gn_low_mem(self, g, n, dirname, use_rs_values=True):
        if use_rs_values:
            self._load_caches(2 * g - 2 + n)

        self._dict[(g, n)] = {}

        graph_file_name = f"{dirname}/{self.name}/G_g={g},n={n}.tsv"

        with open(graph_file_name, "w") as graphs_data:
            for k in range(1, 2 * g - 2 + n + 1):
                graphs_data_writer = csv_writer(graphs_data, delimiter="\t")
                G_list = generate_W_k_graphs(
                    self, g, n, len(self.ramification_profile), k
                )
                for G in G_list:
                    graphs_data_writer.writerow(
                        [G.genera, G.sources, G.loops, G.sings, G.dilatons, G.adj]
                    )

        with open(graph_file_name, "r") as graphs_data:
            graphs_data_reader = csv_reader(graphs_data, delimiter="\t")

            for graph_data in graphs_data_reader:
                genera, sources, loops, sings, dilatons, adj = graph_data

                G = WGraph(
                    eval(genera),
                    eval(sources),
                    eval(loops),
                    eval(sings),
                    eval(dilatons),
                    eval(adj),
                )

                add_contribution_to_poly_dict(
                    self._dict[(g, n)], self, G, self._rs_cache
                )

        return GeneratingFunction(self._dict, "F_" + self.name)

    def _compute_F_gn_graph_list(self, g, n, use_rs_values=True):
        if use_rs_values:
            self._load_caches(2 * g - 2 + n)

        self._dict[(g, n)] = {}
        F_gn_list = []

        for k in range(1, 2 * g - 2 + n + 1):
            G_list = generate_W_k_graphs(self, g, n, len(self.ramification_profile), k)

            for G in G_list:
                G_term = {}
                add_contribution_to_poly_dict(G_term, self, G, self._rs_cache)
                F_gn_list.append(G_term.copy())

        return F_gn_list

    def _load_caches(self, chi):
        if not hasattr(self, "_rs_cache"):
            self._rs_cache = {}

        for (r, s) in self.ramification_profile:
            if (r, s) not in self._rs_cache:
                rs_cache_bound = 0
            else:
                rs_cache_bound = max([2 * g - 2 + n for (g, n) in self._rs_cache[(r, s)]])

            max_d = s * chi - 1 - (chi % 2 == 0)

            if (s + 1) % r == 0:
                max_d //= 2

            if rs_cache_bound < chi + max_d:
                self._rs_cache[(r, s)] = load_F_for_rs_curve(r, s, chi + max_d, False)


class RSCurve(SageObject):
    """
    Spectral curve given by (x = z^r / r, y = z^(s-r)) and the standard bidifferential.
    """

    def __init__(self, r, s, name=None, assume_symmetric=True):
        if name == None:
            name = "r={},s={}".format(r, s)

        self.name = name

        self.r = r
        self.s = s

        self._assume_symmetric = assume_symmetric
        self._dict = {}

    def invariants(self, chi):
        return _convert_F_to_w(self.compute_F(chi))

    def generating_function(self, chi):
        return self.compute_F(chi)

    def compute_F(self, chi):
        self._load_caches(chi)

        for x in range(1, chi + 1):
            for g in range((x + 3) // 2):
                n = x + 2 - 2 * g
                self.compute_F_gn(g, n)

        return GeneratingFunction(self._dict, "F_" + self.name)

    def compute_F_gn(self, g, n):
        self._dict[(g, n)] = {}

        self.graph_cache[(g, n)] = generate_qa_graphs(g, n, self.r)
        self.F_gn_value_cache[(g, n)] = {}

        if self._assume_symmetric:
            possible_x_indices = cached_decreasing_red_partitions_with_min(
                self.grading_condition(g, n), n, self.r, min_val=1
            )
        else:
            possible_x_indices = partitions_with_min(
                self.grading_condition(g, n), n, min_val=1
            )

        for x_indices in possible_x_indices:
            val = self.compute_F_gn_value(g, n, x_indices)

            if val != 0:
                self.F_gn_value_cache[(g, n)][x_indices] = val

                powers = {}
                for l in x_indices:
                    add_to_dict(powers, l, 1)
                term = tuple(("x", (l,), power) for l, power in powers.items())

                self._dict[(g, n)][term] = val / tuple_automorphism_number(x_indices)

        return GeneratingFunction(self._dict, "F_" + self.name)

    def compute_F_gn_value(self, g, n, x_indices):
        val = Integer(0)

        # grading condition
        if sum(x_indices) != self.s * (2 * g - 2 + n) or n < 1:
            return Integer(0)

        x_indices = tuple(sorted(x_indices, reverse=True))

        if (g, n) in self.F_gn_value_cache:
            if x_indices in self.F_gn_value_cache[(g, n)]:
                return self.F_gn_value_cache[(g, n)][x_indices]

        if 2 * g - 2 + n > 1:
            # if s=r+1 then apply string equation
            if self.s == (self.r + 1):
                if x_indices[n - 1] == 1:
                    for j in range(n - 1):
                        if x_indices[j] > self.r:
                            val += x_indices[j] * self.compute_F_gn_value(
                                g,
                                n - 1,
                                tuple(
                                    x_indices[i] - self.r if i == j else x_indices[i]
                                    for i in range(n - 1)
                                ),
                            )
                    return val

            # apply dilaton equation
            for l in range(n):
                if x_indices[l] == self.s:
                    spliced_x_indices = x_indices[:l] + x_indices[l + 1 :]
                    return self.grading_condition(g, n - 1) * self.compute_F_gn_value(
                        g, n - 1, spliced_x_indices
                    )

        for G in self.graph_cache[(g, n)]:
            weight = G.weight(
                self.qairy_coefficient,
                self.grading_condition,
                self.F_gn_value_cache,
                self.r,
                self.s,
                x_indices,
            )

            if weight != 0:
                val += weight / G.automorphism_number

        return val

    # coefficients of quantum Airy structure associated to rs curve
    # C^(j)[q|a] = -1/r * D^(j)_i[k|a]
    @cached_method
    def qairy_coefficient(self, q, j, a):
        @cached_function
        def Psi(r, j, a):
            i = int(len(a) + 2 * j)

            if i > r:
                return Integer(0)

            theta = E(r)

            psi = Integer(0)
            for m in list(permutations([*range(r)], i)):
                term = Integer(1)

                if len(a) > 0:
                    term *= prod(
                        [theta ** (-m[l] * a[l - 2 * j]) for l in range(2 * j, i)]
                    )

                if j > 0:
                    term *= prod(
                        [
                            theta ** (m[2 * l] + m[2 * l + 1])
                            / (theta ** (m[2 * l]) - theta ** (m[2 * l + 1])) ** 2
                            for l in range(j)
                        ]
                    )

                psi += term

            psi /= factorial(i)
            return QQ(psi)

        i, k = compute_rs_qairy_index(self.r, self.s, q)
        l = len(a)

        if l + 2 * j <= i and sum(a) == self.r * (k - i + 1) + self.s * (i - 2 * j - l):
            c = Integer(-1) ** (i - 2 * j - l)
            c *= factorial(i) / (factorial(i - 2 * j - l) * factorial(j) * (2**j))
            c *= Psi(
                self.r,
                j,
                a + tuple(-self.s for t in range(i - 2 * j - l)),
            )

            return c / (-self.r)
        else:
            return Integer(0)

    @cached_method
    def grading_condition(self, g, n):
        return self.s * (2 * g - 2 + n)

    def _load_caches(self, chi):
        # graph_cache loaded as needed
        if not hasattr(self, "graph_cache"):
            self.graph_cache = {}

        if not hasattr(self, "F_gn_value_cache"):
            self.F_gn_value_cache = {}


class AiryCurve(RSCurve):
    """
    Special case in which (r,s) = (2,3) computed using the
    Virasoro constraints.
    """

    def __init__(self):
        super().__init__(2, 3, name="Airy", assume_symmetric=True)

    def compute_F_gn(self, g, n):
        self._load_caches(2 * g - 2 + n)
        self._dict[(g, n)] = {}

        if (g, n) not in [
            (0, 3),
            (1, 1),
        ]:
            self.F_gn_value_cache[(g, n)] = {}

        if self._assume_symmetric:
            possible_x_indices = decreasing_red_partitions_with_min(
                self.grading_condition(g, n), n, 2, min_val=1
            )
        else:
            possible_x_indices = partitions_with_min(
                self.grading_condition(g, n), n, min_val=1
            )

        for x_indices in possible_x_indices:
            self.F_gn_value_cache[(g, n)][x_indices] = self.compute_F_gn_value(
                g, n, x_indices
            )

            powers = {}
            for l in x_indices:
                add_to_dict(powers, l, 1)
            term = tuple(("x", (l,), power) for l, power in powers.items())

            self._dict[(g, n)][term] = self.F_gn_value_cache[(g, n)][
                x_indices
            ] / tuple_automorphism_number(x_indices)

        return GeneratingFunction(self._dict, "F_" + self.name)

    def compute_F_gn_value(self, g, n, d):
        # Step 0: Check if integrating top form
        if sum(d) != self.grading_condition(g, n) or n < 1:
            return Integer(0)

        d = tuple(sorted(d, reverse=True))
        # Step 1: Check if cached (includes initial conditions)
        if (g, n) in self.F_gn_value_cache:
            if d in self.F_gn_value_cache[(g, n)]:
                return self.F_gn_value_cache[(g, n)][d]

        # Step 2: String equation
        if d[n - 1] == 1:
            total = Integer(0)

            for j in range(0, n - 1):
                if d[j] >= 3:
                    d0 = tuple(d[i] - 2 if i == j else d[i] for i in range(0, n - 1))
                    total += d[j] * self.compute_F_gn_value(g, n - 1, d0)
            return total

        # Step 3: Dilaton equation
        if d[n - 1] == 3:
            return self.grading_condition(g, n - 1) * self.compute_F_gn_value(
                g, n - 1, d[: n - 1]
            )

        # Step 4: Virasoro contraints
        if d[n - 1] > 3:
            m = d[n - 1] - 3
            total = Integer(0)

            for j in range(0, n - 1):
                d0 = tuple(d[i] + m if i == j else d[i] for i in range(0, n - 1))
                total += d[j] * self.compute_F_gn_value(g, n - 1, d0)

            for a in range(1, m, 2):
                b = m - a

                d0 = tuple(d[: n - 1]) + (a, b)

                total += self.compute_F_gn_value(g - 1, n + 1, d0) / Integer(2)

                for pair in self.pair_cache[n - 1]:
                    d1 = tuple(d[i] for i in pair[0]) + (a,)
                    d2 = tuple(d[i] for i in pair[1]) + (b,)

                    for g1 in range(0, g + 1):
                        g2 = g - g1
                        total += (
                            self.compute_F_gn_value(g1, len(d1), d1)
                            * self.compute_F_gn_value(g2, len(d2), d2)
                            / Integer(2)
                        )

            return total

    def _load_caches(self, chi):
        if not hasattr(self, "pair_cache"):
            self.pair_cache = {}

        pair_cache_bound = max((i for i in self.pair_cache), default=0)
        for i in range(pair_cache_bound, chi - 2):
            self.pair_cache[i] = get_pairs(tuple(range(i)))

        if not hasattr(self, "F_gn_value_cache"):
            self.F_gn_value_cache = {
                (0, 3): {(1, 1, 1): Integer(1)},
                (1, 1): {(3,): Integer(1) / Integer(8)},
            }


def is_admissible_ramification_profile(ramification_profile):
    return all(
        r >= 1
        and s in range(1, r + 2)
        and ((r + 1) % s == 0 or (r - 1) % s == 0)
        and gcd(r, s) == 1
        for (r, s) in ramification_profile
    )


def _convert_F_to_w(F):
    def convert_x_to_dXi(poly):
        new_poly = {}
        for term in poly:
            x_indices = []
            sym_terms = []
            non_x_part = []
            for var_name, indices, power in term:
                if var_name == "x":
                    x_indices += power * [indices]
                else:
                    non_x_part.append((var_name, indices, power))

            sym_x_indices = list(set(permutations(x_indices)))

            for new_indices in sym_x_indices:
                dXi = [("dXi", (ind, (i + 1,)), 1) for i, ind in enumerate(new_indices)]
                sym_terms.append(tuple(non_x_part + dXi))

            for sym_term in sym_terms:
                new_poly[tuple(sym_term)] = (
                    factorial(len(x_indices)) * poly[term] / len(sym_x_indices)
                )

        return new_poly

    w = F.copy()
    w.name = "w"

    # for (g,n) in w.dict:
    #     for term in w.dict[(g,n)]:
    #         w.dict[(g,n)][term]*=factorial(n)

    for (g, n) in F.dict:
        w.dict[(g, n)] = convert_x_to_dXi(w.dict[(g, n)])

    return w
