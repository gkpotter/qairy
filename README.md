# qairy
qairy is a [SageMath](https://www.sagemath.org/) package written in [Cython](https://cython.org/) (which combines C and Python) for computing topological recursion via higher quantum Airy structures. It can calculate the invariants and generating function from topological recursion for a local spectral curve with arbitrary local data.

## Installation
SageMath is a prerequisite, and instructions for installing it can be found [here](https://doc.sagemath.org/html/en/installation/). Then, qairy can be installed using pip (the Python package manager), via:
```
sage -pip install git+https://github.com/gkpotter/qairy --user
```
Alternatively, this entire repository can be downloaded, which includes examples and computations, and then qairy can be installed via:
```
cd /path/to/qairy
sage -pip install . --user
```
The `--user` flag is optional, and installs qairy in the home folder instead of in the SageMath folder.
## Example
To start a SageMath session, use the command `sage` from the command line. Then, qairy can be imported via:
```
sage: from qairy import *
```
The generating function from topological recursion for a local spectral curve with two ramification points of types (2,3) and (4,5) respectively, can be computed via:
```  
sage: C = SpectralCurve(((2,3),(4,5)))
sage: F = C.generating_function(1)
sage: F
F_rams=((2,3),(4,5))[0, 3] = -1/6 x[1,1]^3 / F01[1,3] - 3/2 x[2,1]^2 x[2,3] / F01[2,5] - 2 x[2,1] x[2,2]^2 / F01[2,5] + 1/2 F01[2,7] x[2,1]^3 / F01[2,5]^2 + 2 F01[2,6] x[2,1]^2 x[2,2] / F01[2,5]^2 - 2/3 F01[2,6]^2 x[2,1]^3 / F01[2,5]^3
F_rams=((2,3),(4,5))[1, 1] = -1/2 F02[1,1,1,1] x[1,1] / F01[1,3] - F02[2,2,1,3] x[2,1] / F01[2,5] - 1/2 F02[2,2,2,2] x[2,1] / F01[2,5] - 2 F02[2,2,1,2] x[2,2] / F01[2,5] - 3/2 F02[2,2,1,1] x[2,3] / F01[2,5] + 2 F01[2,6] F02[2,2,1,2] x[2,1] / F01[2,5]^2 + 3/2 F01[2,7] F02[2,2,1,1] x[2,1] / F01[2,5]^2 + 2 F01[2,6] F02[2,2,1,1] x[2,2] / F01[2,5]^2 - 2 F01[2,6]^2 F02[2,2,1,1] x[2,1] / F01[2,5]^3 - 1/8 x[1,3] / F01[1,3] + 1/8 F01[1,5] x[1,1] / F01[1,3]^2 - 5/8 x[2,5] / F01[2,5] + 5/8 F01[2,9] x[2,1] / F01[2,5]^2 - 3/8 F01[2,7] x[2,3] / F01[2,5]^2 + 3/8 F01[2,7]^2 x[2,1] / F01[2,5]^3 + 2 F01[2,6] F01[2,7] x[2,2] / F01[2,5]^3 + F01[2,6]^2 x[2,3] / F01[2,5]^3 - 3 F01[2,6]^2 F01[2,7] x[2,1] / F01[2,5]^4 - 2 F01[2,6]^3 x[2,2] / F01[2,5]^4 + 2 F01[2,6]^4 x[2,1] / F01[2,5]^5
```
This computation can be saved to JSON for use later via:
```
sage: F.save_to_JSON()
Coefficients of F stored as JSON at: F_g=0,n=3.json
Coefficients of F stored as JSON at: F_g=1,n=1.json
```
It can also be saved to a [Mathematica](https://www.wolfram.com/mathematica/)-compatible plain text expression via:
```
sage: F.save_to_txt()
F polynomial output in plain text at: F_g=0,n=3.txt
F polynomial output in plain text at: F_g=1,n=1.txt
```
There are a variety of other examples in [examples.sage](examples.sage).
