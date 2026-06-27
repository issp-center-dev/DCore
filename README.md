<div align="center">
  <img src="doc/_static/logo_dcore1.png" alt="DCore logo" width="300"><br>
</div>

<div align="center">

[![Tests](https://github.com/issp-center-dev/DCore/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/issp-center-dev/DCore/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/dcore.svg)](https://pypi.org/project/dcore/)
[![Python versions](https://img.shields.io/pypi/pyversions/dcore.svg)](https://pypi.org/project/dcore/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

</div>

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
