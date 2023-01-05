from .misc import (
    add_to_dict,
    get_nested_element,
)

from .partitions import tuple_automorphism_number

from sage.arith.all import factorial
from sage.rings.all import Integer
from sage.all import var, sage_eval

from copy import deepcopy
from sympy.utilities.iterables import multiset_permutations
import json
import os


class GeneratingFunction:
    """
    Encodes the coefficients of a generating function of the form:
    F = sum_(g,n) h^(g-1) F_(g,n)[l_1,...,l_k] x_(l_n) ... x_(l_n)
    """

    def __init__(self, F_dict={}, name="F"):
        self.dict = F_dict
        self.name = name

        self.max_chi = max(
            (2 * g - 2 + n for (g, n) in self.dict),
            default=0,
        )

        self.vars = None
        self.poly = None

    def save_to_txt(self, dirname=""):
        for (g, n) in self.dict:
            self.save_gn_to_txt(g, n, dirname)

    def save_gn_to_txt(self, g, n, dirname=""):
        F_gn_str = "F[{}, {}] = {}".format(
            g, n, poly_dict_to_str(self.dict[(g, n)], False)
        )
        filename = "{}{}F_g={},n={}.txt".format(
            dirname, "/" if dirname != "" else "", g, n
        )

        with open(filename, "w") as outfile:
            outfile.write(F_gn_str)
            print("F polynomial output in plain text at: " + str(outfile.name))

    # Only implemented for the case when all extraneous variables are evaluated.
    def save_to_JSON(self, dirname=""):
        for (g, n) in self.dict:
            self.save_gn_to_JSON(g, n, dirname)

    def save_gn_to_JSON(self, g, n, dirname):
        F_data = {}
        for term, w in self.dict[(g, n)].items():
            val = w

            for var_name, _, power in term:
                if var_name == "x":
                    val *= factorial(power)

            F_data[str(term)] = str(val)

        filename = "{}{}F_g={},n={}.json".format(
            dirname, "/" if dirname != "" else "", g, n
        )

        with open(filename, "w") as outfile:
            json.dump(F_data, outfile, indent=4)
            print("Coefficients of F stored as JSON at: " + str(outfile.name))
            return outfile.name

    def load_from_JSON(self, chi, dirname):

        for x in range(1, chi + 1):
            for g in range((x + 3) // 2):
                n = x + 2 - 2 * g
                filename = "{}{}F_g={},n={}.json".format(
                    dirname, "/" if dirname != "" else "", g, n
                )
                self.load_gn_from_JSON(g, n, filename)

        return self.dict

    def load_gn_from_JSON(self, g, n, filename):
        with open(filename, "r") as infile:
            data = json.load(infile)
            self.dict[(g, n)] = {}
            for s_term, s_val in data.items():
                term = eval(s_term)
                val = sage_eval(s_val)

                for var_name, _, power in term:
                    if var_name == "x":
                        val /= factorial(power)

                self.dict[(g, n)][term] = val

        return self.dict

    def copy(self):
        return GeneratingFunction(deepcopy(self.dict), self.name)

    def set_name(self, name):
        self.name = name
        return self.name

    def get_poly(self, genus=None, num_sources=None):
        self.poly = {(g, n): 0 for (g, n) in self.dict}

        self.vars = {"h": var("h", latex_name="\\hbar")}

        for (g, n) in self.dict:
            for term, w in self.dict[(g, n)].items():
                t = Integer(1)
                for var_name, indices, power in term:
                    full_var_name = var_name + "_" + "_".join(str(i) for i in indices)

                    if full_var_name not in self.vars:
                        self.vars[full_var_name] = var(full_var_name)

                    t *= self.vars[full_var_name] ** power

                self.poly[(g, n)] += w * t

        if genus == None:
            return sum(
                self.poly[(g, n)] * (self.vars["h"] ** (g - 1)) for (g, n) in self.poly
            )
        elif num_sources == None:
            return sum(self.poly[(g, n)] for (g, n) in self.poly if g == genus)
        else:
            return self.poly[(genus, num_sources)]

    def evaluate(self, evaluation_dict):
        evaluated_dict = {}
        for (g, n), term_dict in self.dict.items():
            evaluated_dict[(g, n)] = {}

            for term, w in term_dict.items():
                evaluated_w = w
                evaluated_term = []

                for var_name, indices, power in term:
                    criterion = None

                    if var_name in evaluation_dict:
                        if type(evaluation_dict[var_name]) == tuple:
                            criterion, evaluator = evaluation_dict[var_name]
                        else:
                            evaluator = evaluation_dict[var_name]

                        if criterion != None and not criterion(*indices):
                            evaluated_term.append((var_name, indices, power))
                        else:
                            if type(evaluator) == list:
                                offset_indices = tuple(index - 1 for index in indices)
                                evaluated_w *= (
                                    get_nested_element(evaluator, offset_indices)
                                    ** power
                                )
                            elif callable(evaluator):
                                evaluated_w *= evaluator(*indices) ** power
                            else:
                                raise Exception(
                                    "Invalid coefficient evaluator. Evaluator must be a list or function."
                                )
                    else:
                        evaluated_term.append((var_name, indices, power))

                evaluated_term = tuple(evaluated_term)

                if evaluated_w != 0:
                    add_to_dict(evaluated_dict[(g, n)], evaluated_term, evaluated_w)

            if evaluated_dict[(g, n)] == {}:
                del evaluated_dict[(g, n)]

        self.dict = evaluated_dict

    def relabel_vars(self, relabeling_dict):
        relabeled_dict = {}

        for (g, n), term_dict in self.dict.items():
            relabeled_dict[(g, n)] = {}

            for term, w in term_dict.items():
                relabeled_term = []

                for var_name, indices, power in term:
                    if var_name in relabeling_dict:
                        relabeled_var_name, relabeled_indices = relabeling_dict[
                            var_name
                        ](*indices)
                        relabeled_term.append(
                            (relabeled_var_name, relabeled_indices, power)
                        )
                    else:
                        relabeled_term.append((var_name, indices, power))

                relabeled_term = tuple(relabeled_term)

                add_to_dict(relabeled_dict[(g, n)], relabeled_term, w)

        self.dict = relabeled_dict

    def scale_vars(self, scaling_dict):
        scaled_dict = {}

        for (g, n), term_dict in self.dict.items():
            scaled_dict[(g, n)] = {}

            for term, w in term_dict.items():
                scaling_factor = Integer(1)
                for var_name, indices, power in term:
                    if var_name in scaling_dict:
                        scaling_factor *= scaling_dict[var_name](*indices) ** power

                add_to_dict(scaled_dict[(g, n)], term, w * scaling_factor)

        self.dict = scaled_dict

    def change_vars(self, change_of_var_dict):
        """
        {
            var_name : lambda indices : ((new_var_name, new_indices, new_power), scaling_factor)
        }
        e.g.
        {
            "t" : lambda a : (("x", (2*a+1,), 1), r_factorial(2*a+1, 2))
        }
        """
        new_dict = {}

        for (g, n), term_dict in self.dict.items():
            new_dict[(g, n)] = {}

            for term, w in term_dict.items():
                new_term = []
                total_scaling_factor = Integer(1)

                for var_name, indices, power in term:
                    if var_name in change_of_var_dict:
                        (
                            new_var_name,
                            new_indices,
                            new_power,
                        ), scaling_factor = change_of_var_dict[var_name](*indices)
                        total_scaling_factor *= scaling_factor**power
                        new_term.append((new_var_name, new_indices, new_power * power))
                    else:
                        new_term.append((var_name, indices, power))

                add_to_dict(new_dict[(g, n)], tuple(new_term), w * total_scaling_factor)

        self.dict = new_dict

    def sort_terms(self):
        sorted_dict = {}

        for (g, n), term_dict in self.dict.items():
            sorted_dict[g, n] = {}

            for term, w in term_dict.items():
                sorted_term = tuple(sorted(term))
                sorted_dict[g, n][sorted_term] = w

        self.dict = sorted_dict

    def __call__(self, g, n):
        return self.get_poly(g, n)

    def __eq__(self, other):
        # Sorts terms first!
        self.sort_terms()
        other.sort_terms()

        for (g, n) in self.dict:
            if (g, n) not in other.dict:
                return False
            else:
                if self.dict[(g, n)].keys() != other.dict[(g, n)].keys():
                    return False
                else:
                    for term in self.dict[(g, n)]:
                        if self.dict[(g, n)][term] != other.dict[(g, n)][term]:
                            return False
        return True

    def __str__(self):
        F_str = ""
        for (g, n) in self.dict:
            F_str += "{}[{}, {}] = {}\n".format(
                self.name, g, n, poly_dict_to_str(self.dict[(g, n)])
            )

        return F_str[:-1]

    def __repr__(self):
        return str(self)


def poly_dict_to_str(term_dict, include_newlines=False):
    F_str = ""

    for term, w in term_dict.items():
        if w != 0:
            str_w = str(w)

            if str_w[0] == "-" and len(F_str) > 0:
                F_str = F_str[:-2] + "- "
                str_w = str_w[1:]
            F_str += (
                (str_w + " " if str_w != "1" else "")
                + term_to_str(term)
                + ("\n" if include_newlines else "")
                + " + "
            )

    F_str = F_str[:-3].rstrip(" ")

    return F_str


def term_to_str(term):
    num = ""
    denom = ""
    add_parens = False

    for var_name, indices, power in term:
        if type(indices[0]) == tuple:
            indices_str = "".join(
                str(list(sub_indices)).replace(" ", "") for sub_indices in indices
            )
        else:
            indices_str = str(list(indices)).replace(" ", "")

        item = "{}{}{} ".format(
            var_name, indices_str, "^" + str(abs(power)) if abs(power) > 1 else ""
        )
        if power > 0:
            num += item
        else:
            if denom != "":
                add_parens = True
            denom += item

    denom = denom[:-1]
    num = num[:-1]

    if add_parens:
        denom = "(" + denom + ")"

    term_string = num + (" / " + denom if denom != "" else "")

    return term_string


def load_F_for_rs_curve(r, s, chi, include_permutations=True):
    F = GeneratingFunction()

    F.load_from_JSON(
        chi, os.path.join(os.path.dirname(__file__), "data", f"F_r={r},s={s}")
    )

    F_val_dict = {}
    for (g, n) in F.dict:
        F_val_dict[(g, n)] = {}
        for term, w in F.dict[(g, n)].items():
            x_indices = []
            for var_name, indices, power in term:
                if var_name != "x":
                    raise ("Unevaluated terms error")
                else:
                    index = indices[0]
                    for _ in range(power):
                        x_indices.append(index)
            F_val_dict[(g, n)][tuple(x_indices)] = w * tuple_automorphism_number(x_indices)

    if include_permutations:
        for (g, n) in F_val_dict:
            for b in F_val_dict[(g, n)].copy():
                aut_factor = tuple_automorphism_number(b)
                for b_perm in multiset_permutations(b):
                    F_val_dict[(g, n)][tuple(b_perm)] = F_val_dict[(g, n)][b] * aut_factor

    return F_val_dict
