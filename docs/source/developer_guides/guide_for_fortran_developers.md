# Guide for Fortran developers

## Setting up the developer environment

The only prerequisites before these steps are that you have git installed (which is included on many systems by default) and some version of Conda installed.

1. Clone the repository and enter the cloned project directory.

   ```sh
   git clone https://github.com/golmschenk/eesunhong
   cd eesunhong
   ```

2. Create a Conda virtual environment with the correct dependencies installed.

   ```sh
   conda env create --file=environments.yml
   ```

3. Activate the Conda environment.

   ```sh
   conda activate eesunhong_env
   ```

In the future, whenever working on developing `eesunhong`, activate the existing environment with the above line. When you are done working, you can either deactivate the environment with `conda deactivate` or just close the terminal.

## Building

`eesunhong` uses CMake for building. CMake is similar to Make, but allows builds to be much more toolchain agnostic (compiler, operating system, etc) as well as providing many other useful features (such as automatically downloading dependencies). All commands expected from inside the project directory.

To build, first configure the build with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```

Then build with:

```sh
cmake --build build --config Debug
```

Just the second build command is required for subsequent builds, and will only rebuild changed components.

To clean the build, simply delete the build directory:

```sh
rm -rf build
```

To build an optimized version, first clean the build, then run the configure and build commands again, but this time replacing `Debug` with `Release`.

To compile using `ifort`, make sure `ifort` is available in the path, then run

```sh
export FC=ifort
```

prior to the build commands.

## Adding source files

The file specifying which source files to compile, `CMakeLists.txt`, may look somewhat confusing to those unfamiliar to CMake, especially due to the complexity added by importing remote GitHub dependencies and building both a Python library and the main binary executable. For Fortran developers, the main piece you will likely need to change is to add additional source files to be compiled. In most cases, you should add your new source files to the statement starting with `add_executable(eesunhong_main`. If the source files will include functions or routines that should be callable from Python, the should generally be included in the statement starting with `add_library(eesunhong_internal`.
