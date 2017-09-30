cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(qfitlib Fortran)
option(ENABLE_QFIT_TOOLS "Build standalone QFITLIB tools." OFF)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

option(ENABLE_MPI               "Enable MPI"                    OFF)
option(ENABLE_MPIF              "Use mpif.h instead of mpi.mod" OFF)
option(ENABLE_64BIT_INTEGERS    "Enable 64-bit integers"        OFF)
option(ENABLE_STATIC_LINKING    "Enable static linking"         OFF)

if(ENABLE_64BIT_INTEGERS)
    message("-- Using 64-bit integers")
    if(ENABLE_MPI)
        set(ENABLE_MPIF ON)
        message("-- WARNING: Using mpif.h instead of mpi.mod, which can cause conversion")
        message("            issues if MPI is not built with 64-bit integers.")
    endif()
endif()

if(ENABLE_MPI AND NOT DEFINED HOST_PROGRAM)
    find_package(MPI REQUIRED)
    if(MPI_Fortran_FOUND)
        if(NOT ${MPI_Fortran_COMPILER} MATCHES ${CMAKE_Fortran_COMPILER})
            message("-- WARNING: Detected MPI Fortran compiler does not match Fortran compiler")
#            set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
        endif()
        include_directories(${MPI_Fortran_INCLUDE_PATH})
        add_definitions(-DVAR_MPI)
        if(NOT ENABLE_MPIF)
            add_definitions(-DUSE_MPI_MOD_F90)
        endif()
    # the following use of MPI_ variables is deprecated
    elseif(MPI_FOUND)
        include_directories(${MPI_INCLUDE_PATH})
        add_definitions(-DVAR_MPI)
        if(NOT ENABLE_MPIF)
            add_definitions(-DUSE_MPI_MOD_F90)
        endif()
    else()
        message(FATAL_ERROR "-- MPI Fortran not found")
    endif()
elseif(ENABLE_MPI AND DEFINED HOST_PROGRAM)
    add_definitions(-DVAR_MPI)
    if(NOT ENABLE_MPIF)
        add_definitions(-DUSE_MPI_MOD_F90)
    endif()
endif()

if(NOT DEFINED HOST_PROGRAM)
    message("-- No external host program defined. Locating MATH libraries.")

    #find_package( BLAS REQUIRED )
    find_package( LAPACK REQUIRED )

    set(EXTERNAL_LIBS
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})
else()
    message("-- External host program is: ${HOST_PROGRAM}")
    if("${HOST_PROGRAM}" STREQUAL "DALTON")
        add_definitions(-DPRG_DALTON)
    endif()
endif()

if("${INTEGRAL_LIBRARY}" STREQUAL "GEN1INT")
    add_definitions(-D${INTEGRAL_LIBRARY})
    message("-- Using Gen1Int integral library.")
else()
    message("-- No integral library enabled.")
endif()

if(DEFINED PARENT_INCLUDE_DIR)
    include_directories(${PARENT_INCLUDE_DIR})
endif()

if(NOT DEFINED PARENT_MODULE_DIR)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
else()
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PARENT_MODULE_DIR})
endif()
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

if(DEFINED PARENT_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${PARENT_INSTALL_PREFIX})
endif()

if(NOT DEFINED CMAKE_BUILD_TYPE)
    if(DEFINED PARENT_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE ${PARENT_BUILD_TYPE})
    else()
        set(CMAKE_BUILD_TYPE release CACHE STRING "-- Release build (options are: none release debug profile)" FORCE)
    endif()
else()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "")
        set(CMAKE_BUILD_TYPE release)
    endif()
endif()
message("-- build type is ${CMAKE_BUILD_TYPE}")

if(NOT CMAKE_Fortran_MODULE_DIRECTORY)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/modules)
endif()

include(ConfigCompilerFlags)

if(ENABLE_QFIT_TOOLS)
    message("-- Building QFITLIB standalone tools.")

    set(TOOLS_GRID
        ${PROJECT_SOURCE_DIR}/tools/fortran/grid.F90
    )

    set(TOOLS_MATRICES
        ${PROJECT_SOURCE_DIR}/tools/fortran/matrices.F90
    )

    set(TOOLS_XYZ2GRID
        ${PROJECT_SOURCE_DIR}/tools/fortran/xyz2grid.F90
    )

    add_executable(grid ${TOOLS_GRID})
    target_link_libraries (grid qfitlib ${EXTERNAL_LIBS})

    add_executable(xyz2grid ${TOOLS_XYZ2GRID})
    target_link_libraries (xyz2grid qfitlib ${EXTERNAL_LIBS})

    add_executable(matrices ${TOOLS_MATRICES})
    target_link_libraries (matrices qfitlib ${EXTERNAL_LIBS})
endif(ENABLE_QFIT_TOOLS)

find_package( Doxygen )
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in
                   ${CMAKE_BINARY_DIR}/Doxyfile)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
    message("-- Doxygen found. Run `make doc` to make documentation")
else()
    message("-- Doxygen not found. Run `sudo apt install doxygen doxygen-latex graphviz`")
endif(DOXYGEN_FOUND)

set(SOURCES_LIBRARY
    ${PROJECT_SOURCE_DIR}/src/qfit.F90
    ${PROJECT_SOURCE_DIR}/src/connolly.F90
    ${PROJECT_SOURCE_DIR}/src/linear_solver.F90
    ${PROJECT_SOURCE_DIR}/src/input_readers.F90
    ${PROJECT_SOURCE_DIR}/src/global_variables.F90
    ${PROJECT_SOURCE_DIR}/src/precision.F90
    ${PROJECT_SOURCE_DIR}/src/integrals.F90
    ${PROJECT_SOURCE_DIR}/src/utilities.F90
    ${PROJECT_SOURCE_DIR}/src/io.F90
    ${PROJECT_SOURCE_DIR}/src/printpkg.F
    )

add_library(qfitlib STATIC ${SOURCES_LIBRARY})
install(TARGETS qfitlib ARCHIVE DESTINATION lib)
