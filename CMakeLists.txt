cmake_minimum_required(VERSION 3.15...3.27)
project(GeoPH LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

if (NOT SKBUILD)
    message(WARNING "\
    This CMake file is meant to be executed using 'scikit-build-core'.
    Running it directly will almost certainly not produce the desired
    result. If you are a user trying to install this package, use the
    command below, which will install all necessary build dependencies,
    compile the package in an isolated environment, and then install it.
    =====================================================================
    $ pip install .")
endif()

if (CMAKE_VERSION VERSION_LESS 3.18)
    set(DEV_MODULE Development)
else()
    set(DEV_MODULE Development.Module)
endif()

find_package(Python 3.8
    REQUIRED COMPONENTS Interpreter ${DEV_MODULE}
    OPTIONAL_COMPONENTS Development.SABIModule)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
    # Name of the extension
    geoph_impl

    # Target the stable ABI for Python 3.12+, which reduces
    # the number of binary wheels that must be built. This
    # does nothing on older Python versions
    STABLE_ABI

    # Source code goes here
    python/pyGeoPH.cpp
)

if (WIN32)
    message("Windows: setting GMP/MPFR path")
    set(GMP_INCLUDE_DIR "deps/auxiliary/gmp/include")
    get_filename_component(GMP_LIBRARIES "deps/auxiliary/gmp/lib/libgmp-10.lib" ABSOLUTE)
    set(MPFR_INCLUDE_DIR "deps/auxiliary/gmp/include")
    get_filename_component(MPFR_LIBRARIES "deps/auxiliary/gmp/lib/libmpfr-4.lib" ABSOLUTE)

    file(DOWNLOAD
            "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
            "${PROJECT_SOURCE_DIR}/deps/eigen-3.4.0.tar.gz"
            EXPECTED_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72)
    file(ARCHIVE_EXTRACT
            INPUT "${PROJECT_SOURCE_DIR}/deps/eigen-3.4.0.tar.gz"
            DESTINATION "${PROJECT_SOURCE_DIR}/deps")
    set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/deps/eigen-3.4.0")
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    find_package(Eigen3 REQUIRED NO_MODULE)
    target_link_libraries (geoph_impl PUBLIC Eigen3::Eigen)
endif()

find_package(CGAL REQUIRED)
target_include_directories (geoph_impl PUBLIC ${CGAL_INCLUDE_DIR})
target_link_libraries (geoph_impl PUBLIC ${CGAL_LIBRARIES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_compile_definitions(ENABLE_OPENMP)
endif()

find_package(TBB)
if (TBB_FOUND)
    message("TBB FOUND")
    target_link_libraries (geoph_impl PUBLIC TBB::tbb TBB::tbbmalloc)
    add_compile_definitions(ENABLE_TBB)
    add_compile_definitions(CGAL_LINKED_WITH_TBB)
    add_compile_definitions(TBB_PREVIEW_GLOBAL_CONTROL)
endif()

if (SKBUILD)
    install(TARGETS geoph_impl LIBRARY DESTINATION geoph)
endif()

#if (NOT SKBUILD)
#    add_executable(test src/test.cpp)
#    target_include_directories (test PUBLIC ${CGAL_INCLUDE_DIR})
#    target_link_libraries (test PUBLIC ${CGAL_LIBRARIES})
#    target_include_directories (test PUBLIC ${EIGEN3_INCLUDE_DIR})
#    target_link_libraries (test PUBLIC ${EIGEN3_LIBRARIES})
#    target_link_libraries (test PUBLIC TBB::tbb TBB::tbbmalloc)
#endif()