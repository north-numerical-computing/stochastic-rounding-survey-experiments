# Numerical MATLAB experiments with stochastic rounding
This repository contains the source code for reproducing the results in [Sec. 8, 1]. The following scripts for generating the data presented in the survey are made available:
   * [`test_sum.m`](./test_sum.m) - Experiment with the harmonic sum [Sec. 8a, 1].
   * [`test_matvec.m`](./test_matvec.m) - Experiment with matrix-vector multiplication [Sec. 8a, 1].
   * [`ODE_test.m`](./ODE_tests.m) - Solution of an exponential decay ODE with Euler's method [Sec. 8d, 1].
   * [`unit_circle_ODE.m`](./unit_circle_ODE.m) - Solution of a unit circle ODE with Euler's method [Sec. 8d, 1].
   
The script [`run_tests.m`](./run_tests.m) performs all of the tests from the paper [1] in one run.

These scripts rely on the [chop](https://github.com/higham/chop) library for implementing low precision arithmetics with stochastic rounding and therefore it should be downloaded and placed on the MATLAB search path. These experiments were developed and run on MATLAB version 2021b.

### References

 [1] M. Croci, M. Fasi, N. J. Higham, T. Mary, and M. Mikaitis. [*Stochastic Rounding: Implementation, Error Analysis, and Applications*](https://royalsocietypublishing.org/doi/10.1098/rsos.211631). R. Soc. Open Sci., 9:3. Mar. 2022.

### License

This software is distributed under the terms of the 2-clause BSD software license (see [LICENSE](./LICENSE)).

The MATLAB function `chop` is distributed under the terms of the [BSD 2-Clause "Simplified" License](https://raw.githubusercontent.com/higham/chop/master/license.txt).

The CPFloat C library is distributed under the [GNU Lesser General Public License, Version 2.1 or later](https://raw.githubusercontent.com/mfasi/cpfloat/master/LICENSES/LGPL-2.1-or-later.txt).
