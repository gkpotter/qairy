# qairy
qairy is a SageMath package written in Python and C for computing topological recursion via higher quantum Airy structures. It can calculate the invariants of topological recursion for a given ramification profile and can export to LaTeX, plain text, or Mathematica-compatible expressions.

## Example
Computing the invariants for a curve with one ramificaiton point of type (2,3) and another of type (4,5) can be done as follows.
```  
sage: C = SpectralCurve([(2,3),(4,5)]);
sage: w = C.invariants(3);
sage: w[(1,2)]
25/24 F02[1,1,1,3] dX[1,1] dX[1,1] / F01[1,3]^2 + 3/2 F02[1,1,1,1] dX[1,3] dX[1,1] / F01[1,3]^2 + 3/2 F02[1,1,1,1] dX[1,1] dX[1,3] / F01[1,3]^2 - 13/8 F01[1,5] F02[1,1,1,1] dX[1,1] dX[1,1] / F01[1,3]^3 + 3 F02[2,2,1,7] dX[2,1] dX[2,1] / F01[2,5]^2 + ...
```
## Release date

qairy will be released by the end of 2022.
