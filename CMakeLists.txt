cmake_minimum_required(VERSION 3.17)

project(eesunhong LANGUAGES Fortran CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(cmake/skipBuildBranch.cmake)
include(cmake/fetchStdlib.cmake)
include(cmake/fetchVbbl.cmake)

add_library(polyroots OBJECT third_party/polyroots-fortran/polyroots_cmplx_roots_gen.f90)
add_library(roots OBJECT third_party/roots-fortran/root_module.F90)

add_library(eesunhong_internal OBJECT src/eesunhong_recipes_replacements.f90 src/geo_par.f src/eesunhong_geo_par_c_wrapper.f90)
target_link_libraries(eesunhong_internal PUBLIC fortran_stdlib)

add_library(eesunhong_fortran_library SHARED $<TARGET_OBJECTS:eesunhong_internal> $<TARGET_OBJECTS:roots> $<TARGET_OBJECTS:VBBinaryLensing>)
target_link_libraries(eesunhong_fortran_library PUBLIC fortran_stdlib)
# Without this was defaulting to CXX which was resulting in warnings with compact unwinding.
set_target_properties(eesunhong_fortran_library PROPERTIES LINKER_LANGUAGE Fortran)

add_library(eesunhong_complete_static STATIC $<TARGET_OBJECTS:eesunhong_internal> $<TARGET_OBJECTS:polyroots> $<TARGET_OBJECTS:roots>)
target_link_libraries(eesunhong_complete_static PUBLIC fortran_stdlib)

add_executable(eesunhong_main src/main.f third_party/minuit/minuit_94a_dblb.f src/fcnrvg4_Ctpar.f src/bilens.f src/critical.f
        src/microcurve_rvg4Ctpar.f src/hexadec_only.f src/eesunhong_real_complex_conversion.f90 src/eesunhong_vbbl_interface.f90)
target_link_libraries(eesunhong_main PUBLIC eesunhong_complete_static VBBinaryLensing_static)
# Without this was defaulting to CXX which was resulting in warnings with compact unwinding.
set_target_properties(eesunhong_main PROPERTIES LINKER_LANGUAGE Fortran)

include(cmake/install.cmake)
