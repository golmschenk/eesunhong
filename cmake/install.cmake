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
