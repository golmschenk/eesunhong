if("$ENV{SKIP_CMAKE_BUILD}" STREQUAL "1")
    message(STATUS "Manual build set, skipping build.")
    return()
endif()