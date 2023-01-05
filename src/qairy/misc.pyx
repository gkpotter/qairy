from sage.misc.all import prod
from sage.rings.all import Integer
from sage.arith.all import xgcd
from sage.misc.all import cached_function


def r_factorial(n, r):
    return Integer(prod([*range(n, 0, -r)]))


def add_to_dict(d, k, v):
    if k not in d:
        d[k] = v
    else:
        d[k] += v


def get_nested_element(nested_lists, keys):
    element = nested_lists
    for key in keys:
        element = element[key]

    return element


@cached_function
def compute_rs_qairy_index(r, s, q):
    """
    Finds a solution to q = r * k + (s - r) * (i - 1).
    """

    _, u, v = xgcd(r, s)
    v *= q
    u *= q

    while v > r - 1:
        v -= r
        u += s

    while v < 0:
        v += r
        u -= s

    i = v + 1
    k = u + v

    return (i, k)


def delta(a, b):
    return Integer(a == b)
