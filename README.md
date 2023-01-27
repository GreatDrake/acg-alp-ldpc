# acg-alp-ldpc

Implementation of LDPC codes decoding algorithms:

* Belief Propagation
* [Adaptive Cut Generation Algorithm](https://ieeexplore.ieee.org/document/6218777)
* [QP-ADMM Algorithm](https://arxiv.org/abs/1910.12712)

Also implementation of the optimization of check matrix for QP-ADMM algorithm.

### Installation

To run Linear Programming based methods, ```glpk``` library needs to be installed. Installation instructions (works only on mac):

```
git clone https://github.com/GreatDrake/acg-alp-ldpc.git && cd acg-alp-ldpc
wget http://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz
tar -xf glpk-5.0.tar.gz && cd glpk-5.0
./configure && make
mv src/glpk.h ../ && mv src/.libs/libglpk.a ../ && cd ../
```

On Ubuntu one can evaluate Belief Propagation and QP-ADMM algorithms without installing ```glpk```. For that, after installation above run:

```
sed -i '1 s/^/\/\//' main.cpp
```

### Run experiments

#### Algorithms

To evaluate algorithms on matrix at ```data/optimalH.txt```, run: 

```
make run
```

This command will output algorithm results to stdout and ```report.csv```.

#### Matrix optimization

To run local decsent algorithm for check matrix optimization for QP-ADMM algorithm (see its parameters in ```optimize_H.cpp```), run:

```
make optimize
```

This command will overwrite file ```data/optimalH.txt```. Adjust constant ```THREADS_NUM``` in ```optimize_H.cpp``` before run.

#### QP-ADMM parameters optimization

To run grid-search on QP-ADMM parameters for matrix at ```data/optimalH.txt```, run:

```
make run_qpadmm_params
```

This command will output optimal ```alpha``` and ```mu``` to stdout.
