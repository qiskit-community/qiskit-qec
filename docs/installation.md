Installation guide
==================

## Content
* [Setting up python environment](#setting-up-python-environment)
* [Installing Qiskit QEC](#installing)
* [Installing optional dependencies](#installing-optional-dependencies)

### Setting up python environment

Next, we assume that [Python](https://www.python.org/) 3.6 or higher is installed.  It is recommended to use a Python virtual environment that is dedicated to working with `qrao`.  The steps in the remainder of this tutorial will assume that this environment is activated using either method.


#### Option 1: conda (recommended)

Install [conda](https://docs.conda.io/en/latest/).

The following commands create and activate a conda virtual environment named `qiskit_qec_env`:

```sh
$ conda create -n qiskit_qec_env python=3
$ conda activate qiskit_qec_env
```

Any time you open a new shell and wish to work with Qiskit QEC, you will need to activate it using the second line above.


#### Option 2: venv (included in Python)

You can create and activate a virtual environment with the following commands:

```sh
$ python3 -m venv venv
$ source venv/bin/activate
```

The first command creates a virtual environment in the `venv` folder of the current directory.  We recommend using this name, as it will be ignored by git (i.e., we have added it to `.gitignore`).

Any time you open a new shell and wish to work with Qiskit QEC, you will need to activate it using the second line above.  [If you prefer that the virtual environment be activated automatically any time you enter the directory, you may wish to look into using [direnv](https://direnv.net/) or a similar tool.]

### Installing

1. Clone repo

```shell
git clone https://github.com/qiskit-community/qiskit-qec
cd qiskit-qec
```

2. Install from source using the following command.

```shell
pip install -e .
```

Note that some alternative commands to install from source, or to install directly from github, may not work properly. So if you have problems, make sure to use the exact command shown above.

### Installing optional dependencies

Installing dev dependencies (tox, pylint, etc).

```shell
pip install -r requirements-dev.txt
```

You will probably need jupyter notebook to work on project.

```shell
conda install jupyter
```
