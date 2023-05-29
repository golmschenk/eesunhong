cmake_minimum_required(VERSION 3.17)
project(eesunhong LANGUAGES Fortran)

set(CMAKE_LANG_DISABLE C)
set(CMAKE_LANG_DISABLE CXX)

find_package(Python
        COMPONENTS Interpreter Development.Module
        REQUIRED)

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
enable_testing()

add_library(polyroots OBJECT third_party/polyroots-fortran/polyroots_cmplx_roots_gen.f90)
add_library(roots OBJECT third_party/roots-fortran/root_module.F90)
add_library(eesunhong_recipes_replacements OBJECT eesunhong/eesunhong_recipes_replacements.f90)
target_link_libraries(eesunhong_recipes_replacements PUBLIC fortran_stdlib)

add_library(eesunhong_fortran_library SHARED $<TARGET_OBJECTS:roots> $<TARGET_OBJECTS:eesunhong_recipes_replacements>)
target_link_libraries(eesunhong_fortran_library PUBLIC fortran_stdlib)

add_library(eesunhong_complete_static STATIC $<TARGET_OBJECTS:polyroots> $<TARGET_OBJECTS:roots> $<TARGET_OBJECTS:eesunhong_recipes_replacements>)

add_executable(eesunhong_run eesunhong/main.f third_party/minuit/minuit_94a_dblb.f eesunhong/fcnrvg4_Ctpar.f eesunhong/bilens.f eesunhong/critical.f
        eesunhong/microcurve_rvg4Ctpar.f eesunhong/hexadec_only.f eesunhong/geo_par.f)
target_link_libraries(eesunhong_run PUBLIC fortran_stdlib)
target_link_libraries(eesunhong_run PUBLIC eesunhong_complete_static)

if(SKBUILD)
    set(library_directory "${SKBUILD_PLATLIB_DIR}/lib")
    set(binary_directory "${SKBUILD_SCRIPTS_DIR}")
else()  # Calling CMake directly instead of through scikit-build means this is a developer build.
    set(library_directory "lib")
    set(binary_directory "bin")
    set(CMAKE_INSTALL_PREFIX .)
endif()

install(TARGETS eesunhong_fortran_library DESTINATION ${library_directory})
install(TARGETS eesunhong_run DESTINATION ${binary_directory})
