# shadoks-CGSHOP2026
Reconfiguration of triangulations with parallel flips

Shadoks solvers of the CG:SHOP 2026 competition.

 - exact: Produces exact solutions
 - solver: Produces initial solution
 - improver: Improves an existing solution
 - makelb: Creates the "lb.txt" file that solver and improver use

Linux static precompiled binaries are included. To compile, run
```
cmake .
make
```

Dependencies like EvalMAXSat, CaDiCal, and cxxopts are included. Check description.pdf for a description of the algorithms and heuristics used.
