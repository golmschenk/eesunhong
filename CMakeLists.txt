cmake_minimum_required(VERSION 3.19)

project(eesunhong LANGUAGES Fortran)

# Skip the build if the user is manually building.
if ("$ENV{EESUNHONG_ONLY_INSTALL_PYTHON}" STREQUAL "1")
    message(STATUS "Manual build set, skipping build.")
    return()
endif ()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_TESTING OFF)
include(FetchContent)

# Set build directories
if (SKBUILD)
    set(library_directory "${SKBUILD_PLATLIB_DIR}")
    set(binary_directory "${SKBUILD_PLATLIB_DIR}")
else ()  # Calling CMake directly instead of through scikit-build means this is a developer build.
    set(library_directory ".")
    set(binary_directory ".")
    set(CMAKE_INSTALL_PREFIX .)
endif ()

# Enable Python.
find_package(Python COMPONENTS Interpreter Development REQUIRED)
list(APPEND CMAKE_INSTALL_RPATH ${Python_LIBRARY_DIRS})
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}")
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
list(APPEND CMAKE_INSTALL_RPATH "${library_directory}")
set(CMAKE_MACOSX_RPATH TRUE)
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Enable Rust.
FetchContent_Declare(
        Corrosion
        GIT_REPOSITORY https://github.com/corrosion-rs/corrosion.git
        GIT_TAG v0.4
)
FetchContent_MakeAvailable(Corrosion)

# Build Fortran stdlib.
FetchContent_Declare(
        fortran_stdlib
        GIT_REPOSITORY https://github.com/fortran-lang/stdlib.git
        GIT_TAG fb4ca801f0c8e0ed09f9d137c620676fa348ebdd  # v0.2.1
)
FetchContent_MakeAvailable(fortran_stdlib)
set(BUILD_TESTING ON)

# Build eesunhong.

## Rust build.
corrosion_import_crate(MANIFEST_PATH Cargo.toml)
# Required to make sure Rust uses the correct path to libpython at runtime.
corrosion_add_target_local_rustflags(eesunhong_rust_library "-Clink-args=-Wl,-rpath,${Python_LIBRARY_DIRS}")

## Main Fortran build.
add_library(polyroots OBJECT third_party/polyroots-fortran/polyroots_cmplx_roots_gen.f90)
add_library(roots OBJECT third_party/roots-fortran/root_module.F90)
add_library(eesunhong_recipes_replacements OBJECT src/eesunhong_recipes_replacements.f90)
target_link_libraries(eesunhong_recipes_replacements PUBLIC fortran_stdlib)

add_custom_target(dummy_to_force_build_order)
add_dependencies(dummy_to_force_build_order eesunhong_recipes_replacements)

add_library(eesunhong_fortran_library SHARED $<TARGET_OBJECTS:eesunhong_recipes_replacements> $<TARGET_OBJECTS:roots>)
target_link_libraries(eesunhong_fortran_library PUBLIC fortran_stdlib)

add_library(eesunhong_complete_static STATIC $<TARGET_OBJECTS:eesunhong_recipes_replacements> $<TARGET_OBJECTS:polyroots> $<TARGET_OBJECTS:roots>)
add_dependencies(eesunhong_complete_static dummy_to_force_build_order)
target_link_libraries(eesunhong_complete_static PUBLIC fortran_stdlib)

add_executable(eesunhong_main src/main.f third_party/minuit/minuit_94a_dblb.f src/fcnrvg4_Ctpar.f src/bilens.f src/critical.f
        src/microcurve_rvg4Ctpar.f src/hexadec_only.f src/geo_par.f src/eesunhong_real_complex_conversion.f90 src/eesunhong_astropy_interface.f90)
target_link_libraries(eesunhong_main PUBLIC eesunhong_complete_static eesunhong_rust_library Python::Python)

# Install eesunhong.
install(TARGETS eesunhong_fortran_library DESTINATION ${library_directory})
install(TARGETS eesunhong_rust_library DESTINATION ${library_directory})
install(TARGETS eesunhong_main DESTINATION ${binary_directory})
