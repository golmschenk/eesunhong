Guide for Fortran developers
============================

Setting up the developer environment
------------------------------------

The only prerequisites before these steps are that you have git installed (which is included on many systems by default) and some version of Conda installed.

#. Clone the repository and enter the cloned project directory.

    .. code:: sh

        git clone https://github.com/golmschenk/eesunhong
        cd eesunhong

#. Create a Conda virtual environment with the correct dependencies installed.

    .. code:: sh

        conda env create --file=environments.yml

#. Activate the Conda environment.

    .. code:: sh

        conda activate eesunhong_env

In the future, whenever working on developing ``eesunhong``, activate the existing environment with the above line. When you are done working, you can either deactivate the environment with ``conda deactivate`` or just close the terminal.

Building
--------

``eesunhong`` uses CMake for building. CMake is similar to Make, but allows builds to be much more toolchain agnostic (compiler, operating system, etc) as well as providing many other useful features (such as automatically downloading dependencies). All commands expected from inside the project directory.

To build, first configure the build with:

    .. code:: sh

        cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

Then build with:

    .. code:: sh

        cmake --build build --config Debug

Just the second build command is required for subsequent builds, and will only rebuild changed components.

To clean the build, simply delete the build directory:

    .. code:: sh

        rm -rf build

To build an optimized version, first clean the build, then run the configure and build commands again, but this time replacing ``Debug`` with ``Release``.

To compile using ``ifort``, make sure ``ifort`` is available in the path, then run

    .. code:: sh

        export FC=ifort

prior to the build commands.

