# acg-alp-ldpc

Implementation of [Adaptive Cut Generation Algorithm](https://ieeexplore.ieee.org/document/6218777) for linear programming decoding of LDPC codes.

Comparison with Belief Propagation.

### Run

All experiments (need ```glpk``` installed):

```
make run
```

To install ```glpk``` this instruction can be used:

- Download ```glpk-5.0.tar.gz``` from here: http://ftp.gnu.org/gnu/glpk/.
- Run command ```./configure && make``` from the downloaded directory.
- Move ```src/glpk.h``` and ```src/.libs/libglpk.a``` to directory with program.
