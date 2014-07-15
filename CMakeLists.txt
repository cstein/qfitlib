cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(qfitlib Fortran)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

include(ConfigParentSettings)

if("${INTEGRAL_LIBRARY}" STREQUAL "GEN1INT")
    add_definitions(-D${INTEGRAL_LIBRARY})
    message("-- Using Gen1Int integral library.")
else()
    message(FATAL_ERROR "-- No integral library enabled.")
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

#option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
#option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)

include(ConfigCompilerFlags)
#include(ConfigMath)

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
    )

#set(TOOLS_GRID
#    ${PROJECT_SOURCE_DIR}/tools/grid.F90
#    )

#set(TOOLS_SOLVER
#    ${PROJECT_SOURCE_DIR}/tools/solver.F90
#    )

add_library(qfitlib ${SOURCES_LIBRARY})
#add_executable(grid ${TOOLS_GRID})
#add_executable(solver ${TOOLS_SOLVER})
#target_link_libraries (grid qfitlib ${EXTERNAL_LIBS})
#target_link_libraries (solver qfitlib ${EXTERNAL_LIBS})

install(TARGETS qfitlib ARCHIVE DESTINATION lib)
