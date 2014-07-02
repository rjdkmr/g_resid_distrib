cmake_minimum_required(VERSION 2.8)
project(g_resid_distrib CXX)

set(g_resid_distrib_MAJOR_VERSION 1)
set(g_resid_distrib_MINOR_VERSION 0)
set(g_resid_distrib_PATCH_VERSION 0)
set(g_resid_distrib_VERSION
  ${FOOBAR_MAJOR_VERSION}.${FOOBAR_MINOR_VERSION}.${FOOBAR_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
	set(DEF_INSTALL_CMAKE_DIR CMake)
else()
set(DEF_INSTALL_CMAKE_DIR lib/CMake/g_resid_distrib)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
	set(var INSTALL_${p}_DIR)
	if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	endif()
endforeach()

# set up include-directories
include_directories(
	"${g_resid_distrib_SOURCE_DIR}"   # to find src/*.h
	"${g_resid_distrib_BINARY_DIR}")  # to find src/config.h
	 
# Add sub-directories
add_subdirectory(src)


# ============== Eclipse CDT Error Handling======
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_VERBOSE_MAKEFILE ON)

# ===============================
 
# Add all targets to the build-tree export set
export(TARGETS g_resid_distrib FILE "${PROJECT_BINARY_DIR}/g_resid_distribTargets.cmake")
    
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE g_resid_distrib)


# Create the g_resid_distribConfig.cmake and g_resid_distribConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")

# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(g_resid_distribConfig.cmake.in "${PROJECT_BINARY_DIR}/g_resid_distribConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${g_resid_distrib_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(g_resid_distribConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/g_resid_distribConfig.cmake" @ONLY)
# ... for both
configure_file(g_resid_distribConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/g_resid_distribConfigVersion.cmake" @ONLY)

# Install the g_resid_distribConfig.cmake and g_resid_distribConfigVersion.cmake
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/g_resid_distribConfig.cmake" 
              "${PROJECT_BINARY_DIR}/g_resid_distribConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
	   
# Install the export set for use with the install-tree
install(EXPORT g_resid_distribTargets DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
