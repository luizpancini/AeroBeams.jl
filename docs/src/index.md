# AeroBeams.jl

Documentation for [AeroBeams](https://github.com/luizpancini/AeroBeams.jl).

## Overview
AeroBeams is a finite-element implementation of the geometrically exact beam theory of [[1]](@ref References), augmented with aerodynamic formulations in order to solve aeroelastic problems. The structural part of the code was developed based on the works of [[2]](@ref References), [[3]](@ref References) and [[4]](@ref References), whereas the aerodynamic part follows [[5]](@ref References), [[6]](@ref References) and [[7]](@ref References).

## Installation

To install AeroBeams, simply go to the package manager mode in the Julia REPL by typing ], and then
```julia-repl
pkg> add AeroBeams
```

## References
[1] Hodges, D. H. "Nonlinear Composite Beam Theory". 2006. American Institute of Aeronautics and Astronautics. [10.2514/4.866821](https://doi.org/10.2514/4.866821)

[2] Hodges, D. H., Shang, X. and Cesnik, C. E. S. "Finite element solution of nonlinear intrinsic equations for curved composite beams". 1996. Journal of the American Helicopter Society. [10.2514/6.1995-1174](https://doi.org/10.2514/6.1995-1174)

[3] Yu, W. and Blair, M. "GEBT: A general-purpose nonlinear analysis tool for composite beams". 2012. Composite Structures. [10.1016/j.compstruct.2012.04.007](https://doi.org/10.1016/j.compstruct.2012.04.007)

[4] Wang, Q. and Yu, W. "Geometrically nonlinear analysis of composite beams using Wiener-Milenković parameters". 2017. Journal of Renewable and Sustainable Energy. [10.1063/1.4985091](https://doi.org/10.1063/1.4985091)

[5] Leishman, J. G. "Principles of Helicopter Aerodynamics". 2006. Cambridge University Press.

[6] Peters, D. A., Karunamoorthy, S. and Cao, W. "Finite state induced flow models. I: Two-dimensional thin airfoil". 1995. Journal of Aircraft. [10.2514/3.46718](https://doi.org/10.2514/3.46718)

[7] dos Santos, L. G. P. and Marques, F. D. "Improvements on the Beddoes-Leishman dynamic stall model for low speed applications". 2021. Journal of Fluids and Structures. [10.1016/j.jfluidstructs.2021.103375](https://doi.org/10.1016/j.jfluidstructs.2021.103375)

