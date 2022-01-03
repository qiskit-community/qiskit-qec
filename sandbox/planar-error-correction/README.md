# Software for planar, fault-tolerant quantum error correction

This repository contains Python and C++ software for studying planar quantum codes with circuit-level noise using Qiskit. The software has objects and methods to:

* define families of codes
* construct fault-tolerant circuits
* construct decoding graphs
* process syndrome data
* interact with Qiskit backends/simulators
* enumerate and propagate faults in circuits

## Quantum codes, circuits, and layouts

* rotated subsystem surface code (three layouts, including heavy-hexagon lattice)
* heavy-hexagon compass code (heavy-hexagon lattice)

## Syndrome processing algorithms

* matching decoder using networkx
* matching decoder using pymatching


## Installing plerco (planar-error-correction)
`plerco` is currently not available pypi so you have to build from source. 

### Conda Environment
Cool kids don't put everything in their local Python. Be a cool kid.  Use conda to make a new virtual environment into which you will install `plerco`/whatever else you need. 
0. Install anaconda: https://www.anaconda.com/
1. `conda create -y -n <conda-name> python=3.9`

### 1. Cloning the repo
We use `--recurse-submodules` to make sure `pybind11` gets downloaded as a submodule
```buildoutcfg
git clone --recurse-submodules git://github.com/foo/bar.git
cd planar-error-correction
```

### 2. Building from source
```buildoutcfg
pip install . 
```

### 3. Jupyter Notebook (Optional)
Feel free to run the `example.ipynb` found in the `example` folder. 

### Installing Notes
* If you want to re-install `plerco` after installing in editable mode, please run `./rm_install.sh` to remove any leftover installation files. 
* If you wish to be able to edit the python without reinstalling each time, run `pip install -e .`

## Running Simulation Runner

After pip installing `plerco` in your conda environment, simulation parameters can be written into configuration files and run.
See ``multi.yaml`` and ``runyaml.py`` in the ``simulation_runner`` folder for an example.

```
usage: runyaml [-h] [--log LOG] [--prefix PREFIX] yaml_file

Run a set of simulations.

positional arguments:
  yaml_file        input YAML 1.0 file

optional arguments:
  -h, --help       show this help message and exit
  --log LOG        debug logging to file LOG (default: None)
  --prefix PREFIX  filename prefix for output files (default: )
```

