# Import qairy!
from qairy import *

# Compute the invariants and generating function of the Airy curve
# for 2*g - 2 + n <= 3:
C_Airy = AiryCurve()
w_Airy = C_Airy.invariants(3)
F_Airy = C_Airy.generating_function(3)

# Compute the generating function of (r,s)-curves:
RSCurve(4, 3).generating_function(3)

# The (2,3)-curve is the Airy curve:
assert RSCurve(2, 3).generating_function(7) == C_Airy.generating_function(7)

# Compute the generating function for an arbitrary spectral curve
# given a ramification profile:
C_weierstrass = LocalSpectralCurve(((2, 3), (2, 3), (2, 3)), name="weierstrass")
F_weierstrass = C_weierstrass.generating_function(3)

# Change, scale, and evaluate variables by using:
r = 4
Phi = RSCurve(r, r + 1).generating_function(4)
Phi.change_vars(
    {"x": lambda l: (("t", ((l - 1) // r, (l - 1) % r), 1), 1 / r_factorial(l, r))}
)
Phi.scale_vars(
    {"t": lambda a, m: 1 / (-1 / r) ** (((r - 2) * (1 - a) - 3 * m) / (2 * r + 2))}
)
Phi.evaluate({"t": (lambda a, m: a > 0, lambda a, m: 0)})

# A local spectral curve with a single simple ramification point and
# local data of F01[l] = delta(l, 3) and F02[l1,l2] = 0.
C_simple = LocalSpectralCurve(((2, 3),))
F_simple = C_simple.generating_function(3)

F_simple.evaluate(
    {"F01": lambda a, l: -1 if l == 3 else 0, "F02": lambda a1, a2, l1, l2: 0}
)
F_simple.change_vars({"x": lambda a, l: (("x", (l,), 1), 1)})

assert F_simple == RSCurve(2, 3).generating_function(3)

# Save to JSON (which can be re-imported):
LocalSpectralCurve(((2, 3), (2, 3))).generating_function(2).save_to_JSON(
    "./computations/F_rams=((2,3),(2,3))"
)

# Load from JSON:
GeneratingFunction().load_from_JSON(2, "./computations/F_rams=((2,3),(2,3))")

# Save to plain text (compatible with Mathematica):
LocalSpectralCurve(((2, 3), (2, 3))).generating_function(2).save_to_txt(
    "./computations/F_rams=((2,3),(2,3))"
)
