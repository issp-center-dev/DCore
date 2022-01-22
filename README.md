[![Build Status](https://travis-ci.org/issp-center-dev/DCore.svg?branch=master)](https://travis-ci.org/issp-center-dev/DCore)
# DCore

DCore is aimed at model calculations and ab-initio calculations by the dynamical mean-field theory (DMFT). This package consists of programs with text-based and hdf5-based interface. These programs enable users to perform DMFT calculations and analyze results without writing computer code. ALPS and TRIQS impurity solvers are supported.

# Install

```
> pip3 install dcore
```

## Documentation

https://issp-center-dev.github.io/DCore/index.html

See the the link above for documentations including installation, tutorial, reference manual, and FAQ/Troubleshooting.

## Related paper

Technical details are described in the following paper:

- *"DCore: Integrated DMFT software for correlated electrons"*,  
  H. Shinaoka, J. Otsuki, M. Kawamura, N. Takemori, K. Yoshimi,  
  [SciPost Phys. 10, 117 (2021)](https://scipost.org/10.21468/SciPostPhys.10.5.117)

## Run tests (only for developers)

```
> pytest tests/non-mpi/*/*.py
> mpirun -np 2 pytest tests/mpi/*/*.py
```
