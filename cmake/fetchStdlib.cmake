set(BUILD_TESTING OFF)
include(FetchContent)
FetchContent_Declare(
        fortran_stdlib
        GIT_REPOSITORY https://github.com/fortran-lang/stdlib.git
        GIT_TAG fb4ca801f0c8e0ed09f9d137c620676fa348ebdd  # v0.2.1
)
FetchContent_MakeAvailable(fortran_stdlib)
set(BUILD_TESTING ON)