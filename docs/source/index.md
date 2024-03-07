% eesunhong documentation master file, created by
% sphinx-quickstart on Mon Jan  1 22:48:08 2024.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# Welcome to eesunhong's documentation!

## Installation

To install `eesunhong` use

```shell
pip install eesunhong
```

Although not required, as is generally good practice for any development project, we highly recommend creating a separate virtual environment for each distinct project. For example, via Conda, creating a virtual environment for a project using `eesunhong` might look like

```
conda create -n eesunhong_env python=3.11
```

Then before working, be sure to activate your environment with

```shell
conda activate eesunhong_env
```

Then install `eesunhong` within this environment.

```{toctree}
:caption: 'Contents:'
:maxdepth: 2

guides/guide_for_fortran_developers
```

# Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
