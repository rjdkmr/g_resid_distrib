
# Checking GROMACS PAth
if (NOT DEFINED GMX_PATH AND (NOT DEFINED GMX_LIB OR NOT DEFINED GMX_INCLUDE) )
message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs \nOR\n -DGMX_LIB=/path/to/gromacs/lib -DGMX_INCLUDE=/path/to/gromacs/include")
endif()


if(DEFINED GMX_PATH)
	set(GMX_INCLUDE "${GMX_PATH}/include")
	set(GMX_LIB "${GMX_PATH}/lib")
endif()

include_directories("${GMX_INCLUDE}/gromacs")


# Searching FFTW library
find_library(LIBFFTW fftw3f HINTS ${FFTW_LIB})
if(DEFINED LIBFFTW-NOTFOUND)
message(FATAL_ERROR "\nFFTW library file libfftw3f.so or libfftw3f.a not found at standard location\nUSE:\n -DFFTW_LIB=/path/to/fftw3/lib")
endif()

# Searching GROMACS Library
find_library(LIBGMX gmx HINTS ${GMX_LIB})
find_library(LIBGMXANA gmxana HINTS ${GMX_LIB})
find_library(LIBMD md HINTS ${GMX_LIB})

# For Library
#install(TARGETS gmxhole
#		EXPORT gmxholeTargets
#		RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
#		LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
#		PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/g_resid_distrib"
#		COMPONENT dev)

add_executable(g_resid_distrib  g_resid_distrib.hpp g_residue_distribution.cpp g_resid_distrib_aux.cpp)


target_link_libraries(g_resid_distrib m ${LIBFFTW} ${LIBGMX} ${LIBGMXANA} ${LIBMD} pthread dl)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

install(TARGETS g_resid_distrib
  # IMPORTANT: Add the bar executable to the "export-set"
        EXPORT g_resid_distribTargets
	  	RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
