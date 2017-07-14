
## `standard.py`

### Usage

``` bash
$ pytriqs standard.py input-file
```

### Input file format

Exsample
```
t = 1.0
seedname = test
U = 4.0
lattice = square
model = t2g
nelec = 1.0
nk = 8
```

Parameters

- `lattice` (default `chain`): Lattice. Chosen from `chain`, `square`, `cubic`, and `bethe`
- `model` (default `single`): Structure of a shell. Chosen from `single`, `eg`, `t2g`, and `full-d`.
- `nelec` (default `1.0`): Number of electrons per uunit cell.
- `nk` (default `8`): Number of *k* points along single axis.
- `seedname` (default `pydmft`): The title of system
- `t` (default: `1.0`): Nearest neighbor hopping
- `t'` (default `0.0`): Second nearest hopping
- `U` (default `0.0`)`: Coulomb parameter U
- `J` (default `0.0`)`: Coulomb parameter J

### Output files

#### `seedname`

General-Hk formatted temporary file.

#### `seedname.h5`

HDF5 formatted file which contains the SumkDFT data and the interaction.
