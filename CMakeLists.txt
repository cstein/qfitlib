cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(qfitlib Fortran)
option(ENABLE_QFIT_TOOLS "Build standalone QFITLIB tools." OFF)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

include(ConfigParentSettings)

if(NOT DEFINED HOST_PROGRAM)
    message("-- No external host program defined. Locating MATH libraries.")

    #find_package( BLAS REQUIRED )
    find_package( LAPACK REQUIRED )

    set(EXTERNAL_LIBS
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES})
endif()

if("${INTEGRAL_LIBRARY}" STREQUAL "GEN1INT")
    add_definitions(-D${INTEGRAL_LIBRARY})
    message("-- Using Gen1Int integral library.")
else()
    message("-- No integral library enabled.")
endif()

if(NOT CMAKE_Fortran_MODULE_DIRECTORY)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/modules)
endif()

include_directories(
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

include(ConfigCompilerFlags)

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

add_library(qfitlib ${SOURCES_LIBRARY})

if(ENABLE_QFIT_TOOLS)
    message("-- Building QFITLIB standalone tools.")

    set(TOOLS_GRID
        ${PROJECT_SOURCE_DIR}/tools/grid.F90
    )

    set(TOOLS_MATRICES
        ${PROJECT_SOURCE_DIR}/tools/matrices.F90
    )

    set(TOOLS_XYZ2GRID
        ${PROJECT_SOURCE_DIR}/tools/xyz2grid.F90
    )

    add_executable(grid ${TOOLS_GRID})
    target_link_libraries (grid qfitlib ${EXTERNAL_LIBS})

    add_executable(xyz2grid ${TOOLS_XYZ2GRID})
    target_link_libraries (xyz2grid qfitlib ${EXTERNAL_LIBS})

    add_executable(matrices ${TOOLS_MATRICES})
    target_link_libraries (matrices qfitlib ${EXTERNAL_LIBS})
endif(ENABLE_QFIT_TOOLS)

install(TARGETS qfitlib ARCHIVE DESTINATION lib)
