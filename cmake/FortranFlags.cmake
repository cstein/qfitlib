#if(NOT DEFINED DEFAULT_Fortran_FLAGS_SET)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    set(CMAKE_Fortran_FLAGS         "-fbacktrace -cpp")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize")
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
    endif()
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fbounds-check -Wall -Wunderflow -Wextra"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
    set(CMAKE_Fortran_FLAGS         "-cpp -fno-second-underscore -ftrace=full")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fsloppy-char")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Wall -fbounds-check"
            )
    endif()
    if(ENABLE_STATIC_LINKING)
        message("-- Static linking is not available using g95 compilers")
    endif()
    if(ENABLE_CODE_COVERAGE)
        message("-- Code coverage is not available using g95 compilers")
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-fpp -assume byterecl -traceback")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    if(DEFINED MKL_FLAG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_FLAG}")
    endif()
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static-libgcc -static-intel"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -check all -debug all -fpstkchk -ftrapuv"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        message("-- Code coverage is not available using Intel compilers")
    endif()
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(CMAKE_Fortran_FLAGS         "-mcmodel=medium -Mpreprocess")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast -Munroll")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8 -i8storage"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        message("-- Bounds checking is not available using PGI compilers")
    endif()
    if(ENABLE_CODE_COVERAGE)
        message("-- Code coverage is not available using PGI compilers")
    endif()
    if(ENABLE_STATIC_LINKING)
        message("-- Static linking is not available using PGI compilers")
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname -qfree -qlanglvl=extended -qinit=f90ptr")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-qstrict -O3")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -q64"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -C"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        message("-- Code coverage is not available using XL compilers")
    endif()
    if(ENABLE_STATIC_LINKING)
        message("-- Static linking is not available using XL compilers")
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Cray) 
    add_definitions(-DVAR_CRAY)
    set(CMAKE_Fortran_FLAGS         "-eZ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE " ")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -s integer64"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -R bps"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        message("-- Code coverage is not available using Cray compilers")
    endif()
    if(ENABLE_STATIC_LINKING)
        message("-- Static linking is not available using Cray compilers")
    endif()
endif()

if(DEFINED EXTRA_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_Fortran_FLAGS}")
endif()

#save_compiler_flags(Fortran)

#endif()
