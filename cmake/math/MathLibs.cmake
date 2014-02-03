if(ENABLE_64BIT_INTEGERS)
    set(MATH_LIB_SEARCH_ORDER MKL ACML)
    if(NOT ENABLE_BUILTIN_BLAS AND NOT ENABLE_BUILTIN_LAPACK)
        if(NOT HAVE_BLAS AND NOT HAVE_LAPACK)
            message(STATUS "Since you specified 64bit integers the math lib search order is (only) ${MATH_LIB_SEARCH_ORDER}")
            message(STATUS "This is because apart from MKL and ACML default math library installations are built for 32bit integers")
            message(STATUS "If you know that the library you want to use provides 64bit integers, you can select the library")
            message(STATUS "with -D BLAS_TYPE=X or -D LAPACK_TYPE X (X: MKL ESSL ATLAS ACML SYSTEM_NATIVE)")
            message(STATUS "or by redefining MATH_LIB_SEARCH_ORDER")
        endif()
    endif()
else()
    set(MATH_LIB_SEARCH_ORDER MKL ESSL ATLAS ACML SYSTEM_NATIVE)
    if(NOT ENABLE_BUILTIN_BLAS AND NOT ENABLE_BUILTIN_LAPACK)
        if(NOT HAVE_BLAS AND NOT HAVE_LAPACK)
            message(STATUS "Math lib search order is ${MATH_LIB_SEARCH_ORDER}")
            message(STATUS "You can select a specific type by defining for instance -D BLAS_TYPE=ATLAS or -D LAPACK_TYPE=ACML")
            message(STATUS "or by redefining MATH_LIB_SEARCH_ORDER")
        endif()
    endif()
endif()

#-------------------------------------------------------------------------------
# SYSTEM_NATIVE

set(SYSTEM_NATIVE_BLAS_INCLUDE_PATH_SUFFIXES)
set(SYSTEM_NATIVE_LAPACK_INCLUDE_PATH_SUFFIXES)

set(SYSTEM_NATIVE_BLAS_HEADERS   cblas.h)
set(SYSTEM_NATIVE_LAPACK_HEADERS clapack.h)

set(SYSTEM_NATIVE_BLAS_LIBRARY_PATH_SUFFIXES)
set(SYSTEM_NATIVE_LAPACK_LIBRARY_PATH_SUFFIXES)

set(SYSTEM_NATIVE_BLAS_LIBS   blas)
set(SYSTEM_NATIVE_LAPACK_LIBS lapack)

#-------------------------------------------------------------------------------
# ESSL

set(ESSL_BLAS_INCLUDE_PATH_SUFFIXES)
set(ESSL_LAPACK_INCLUDE_PATH_SUFFIXES)

set(ESSL_BLAS_HEADERS)
set(ESSL_LAPACK_HEADERS)

set(ESSL_BLAS_LIBRARY_PATH_SUFFIXES)
set(ESSL_LAPACK_LIBRARY_PATH_SUFFIXES)

set(ESSL_BLAS_LIBS   essl)
set(ESSL_LAPACK_LIBS essl)

#-------------------------------------------------------------------------------
# ACML

set(ACML_BLAS_INCLUDE_PATH_SUFFIXES)
set(ACML_LAPACK_INCLUDE_PATH_SUFFIXES)

set(ACML_BLAS_HEADERS   cblas.h)
set(ACML_LAPACK_HEADERS clapack.h)

set(ACML_BLAS_LIBRARY_PATH_SUFFIXES   libso)
set(ACML_LAPACK_LIBRARY_PATH_SUFFIXES libso)

set(ACML_BLAS_LIBS   acml)
set(ACML_LAPACK_LIBS acml)

#-------------------------------------------------------------------------------
# ATLAS

set(ATLAS_BLAS_INCLUDE_PATH_SUFFIXES   atlas)
set(ATLAS_LAPACK_INCLUDE_PATH_SUFFIXES atlas)

set(ATLAS_BLAS_HEADERS   cblas.h)
set(ATLAS_LAPACK_HEADERS clapack.h)

set(ATLAS_BLAS_LIBRARY_PATH_SUFFIXES   atlas atlas-base atlas-base/atlas atlas-sse3)
set(ATLAS_LAPACK_LIBRARY_PATH_SUFFIXES atlas atlas-base atlas-base/atlas atlas-sse3)

set(ATLAS_BLAS_LIBS   f77blas cblas atlas)
set(ATLAS_LAPACK_LIBS atlas lapack)

#-------------------------------------------------------------------------------
# MKL

set(MKL_BLAS_INCLUDE_PATH_SUFFIXES)
set(MKL_LAPACK_INCLUDE_PATH_SUFFIXES)

set(MKL_BLAS_HEADERS   mkl_cblas.h)
set(MKL_LAPACK_HEADERS mkl_clapack.h)

if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(MKL_BLAS_LIBRARY_PATH_SUFFIXES   intel64 em64t)
    set(MKL_LAPACK_LIBRARY_PATH_SUFFIXES intel64 em64t)
else()
    set(MKL_BLAS_LIBRARY_PATH_SUFFIXES   ia32 32)
    set(MKL_LAPACK_LIBRARY_PATH_SUFFIXES ia32 32)
endif()

if(ENABLE_THREADED_MKL)
    set(_thread_lib)
    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        set(_thread_lib mkl_intel_thread)
    endif()
    if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
        set(_thread_lib mkl_pgi_thread)
    endif()
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(_thread_lib mkl_gnu_thread)
    endif()
else()
    set(_thread_lib mkl_sequential)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(_compiler_mkl_interface mkl_intel)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(_compiler_mkl_interface mkl_intel)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(_compiler_mkl_interface mkl_gf)
endif()
if(DEFINED MKL_FLAG)
    if(NOT DEFINED _compiler_mkl_interface)
         message(FATAL_ERROR "compiler MKL interface not defined for your compiler")
    endif()
endif()

set(_lib_suffix)
if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    if(ENABLE_64BIT_INTEGERS)
        set(_lib_suffix _ilp64)
    else()
        set(_lib_suffix _lp64)
    endif()
endif()

if(ENABLE_SCALAPACK)
    set(_scalapack_lib      mkl_scalapack${_lib_suffix})
    set(_blacs_intelmpi_lib mkl_blacs_intelmpi${_lib_suffix})
    set(_blacs_sgimpt_lib   mkl_blacs_sgimpt${_lib_suffix})
else()
    set(_scalapack_lib)
    set(_blacs_intelmpi_lib)
    set(_blacs_sgimpt_lib)
endif()

# first try this MKL BLAS combination with SGI MPT
set(MKL_BLAS_LIBS  ${_scalapack_lib} ${_compiler_mkl_interface}${_lib_suffix} ${_thread_lib} mkl_core ${_blacs_sgimpt_lib}   guide pthread m)
# now with Intel MPI
set(MKL_BLAS_LIBS2 ${_scalapack_lib} ${_compiler_mkl_interface}${_lib_suffix} ${_thread_lib} mkl_core ${_blacs_intelmpi_lib} guide pthread m)
# newer MKL BLAS versions do not have libguide
set(MKL_BLAS_LIBS3 ${_scalapack_lib} ${_compiler_mkl_interface}${_lib_suffix} ${_thread_lib} mkl_core ${_blacs_sgimpt_lib}         pthread m)
# now with Intel MPI
set(MKL_BLAS_LIBS4 ${_scalapack_lib} ${_compiler_mkl_interface}${_lib_suffix} ${_thread_lib} mkl_core ${_blacs_intelmpi_lib}       pthread m)
# ancient MKL BLAS
set(MKL_BLAS_LIBS5 mkl guide m)

set(MKL_LAPACK_LIBS mkl_lapack95${_lib_suffix} ${_compiler_mkl_interface}${_lib_suffix})

# older MKL LAPACK
set(MKL_LAPACK_LIBS2 mkl_lapack)

unset(_lib_suffix)
unset(_thread_lib)
unset(_compiler_mkl_interface)
unset(_scalapack_lib)
unset(_blacs_intelmpi_lib)
unset(_blacs_sgimpt_lib)
