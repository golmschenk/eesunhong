cmake_minimum_required(VERSION 3.17)
project(eesunhong)
enable_language(Fortran)

set(BUILD_TESTING OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(FetchContent)
FetchContent_Declare(
        fortran_stdlib
        GIT_REPOSITORY https://github.com/fortran-lang/stdlib.git
        GIT_TAG fb4ca801f0c8e0ed09f9d137c620676fa348ebdd  # v0.2.1
)
FetchContent_MakeAvailable(fortran_stdlib)
set(BUILD_TESTING ON)
enable_testing()

add_library(polyroots third_party/polyroots-fortran/polyroots_cmplx_roots_gen.f90)
add_library(roots third_party/roots-fortran/root_module.F90)
add_library(eesunhong SHARED eesunhong/eesunhong_recipes_replacements.f90)
target_link_libraries(eesunhong PUBLIC fortran_stdlib)
target_link_libraries(eesunhong PRIVATE roots)

add_executable(minuit_all_rvg4Ctpar.xO eesunhong/main.f third_party/minuit/minuit_94a_dblb.f eesunhong/fcnrvg4_Ctpar.f eesunhong/bilens.f eesunhong/critical.f
        eesunhong/microcurve_rvg4Ctpar.f eesunhong/hexadec_only.f eesunhong/geo_par.f)
target_link_libraries(minuit_all_rvg4Ctpar.xO PRIVATE polyroots)
target_link_libraries(minuit_all_rvg4Ctpar.xO PRIVATE eesunhong)

add_custom_command(TARGET minuit_all_rvg4Ctpar.xO
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:minuit_all_rvg4Ctpar.xO> ${CMAKE_SOURCE_DIR}/build/$<TARGET_FILE_NAME:minuit_all_rvg4Ctpar.xO>)

get_target_property(library_output_file_name eesunhong LIBRARY_OUTPUT_NAME)
add_custom_command(TARGET eesunhong
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:eesunhong> ${CMAKE_SOURCE_DIR}/build/$<TARGET_FILE_NAME:eesunhong>)

find_package(Pytest)
pytest_discover_tests(
        PythonTest
        LIBRARY_PATH_PREPEND $<TARGET_FILE_DIR:eesunhong>
        PYTHON_PATH_PREPEND $<TARGET_FILE_DIR:eesunhong>
        DEPENDS eesunhong
)
