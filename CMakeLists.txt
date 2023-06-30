cmake_minimum_required(VERSION 3.17)

project(src/eesunhong LANGUAGES Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(BUILD_TESTING OFF)
include(FetchContent)
FetchContent_Declare(
        fortran_stdlib
        GIT_REPOSITORY https://github.com/fortran-lang/stdlib.git
        GIT_TAG fb4ca801f0c8e0ed09f9d137c620676fa348ebdd  # v0.2.1
)
FetchContent_MakeAvailable(fortran_stdlib)
set(BUILD_TESTING ON)

add_library(polyroots OBJECT third_party/polyroots-fortran/polyroots_cmplx_roots_gen.f90)
add_library(roots OBJECT third_party/roots-fortran/root_module.F90)
add_library(eesunhong_recipes_replacements OBJECT src/eesunhong_recipes_replacements.f90)
target_link_libraries(eesunhong_recipes_replacements PUBLIC fortran_stdlib)

add_library(eesunhong_fortran_library SHARED $<TARGET_OBJECTS:eesunhong_recipes_replacements> $<TARGET_OBJECTS:roots>)
target_link_libraries(eesunhong_fortran_library PUBLIC fortran_stdlib)

add_library(eesunhong_complete_static STATIC $<TARGET_OBJECTS:eesunhong_recipes_replacements> $<TARGET_OBJECTS:polyroots> $<TARGET_OBJECTS:roots>)
target_link_libraries(eesunhong_complete_static PUBLIC fortran_stdlib)

add_executable(eesunhong_main src/main.f third_party/minuit/minuit_94a_dblb.f src/fcnrvg4_Ctpar.f src/bilens.f src/critical.f
        src/microcurve_rvg4Ctpar.f src/hexadec_only.f src/geo_par.f src/eesunhong_real_complex_conversion.f90)
target_link_libraries(eesunhong_main PUBLIC eesunhong_complete_static)

if(SKBUILD)
    set(library_directory "${SKBUILD_PLATLIB_DIR}")
    set(binary_directory "${SKBUILD_PLATLIB_DIR}")
else()  # Calling CMake directly instead of through scikit-build means this is a developer build.
    set(library_directory ".")
    set(binary_directory ".")
    set(CMAKE_INSTALL_PREFIX .)
endif()

install(TARGETS eesunhong_fortran_library DESTINATION ${library_directory})
install(TARGETS eesunhong_main DESTINATION ${binary_directory})
