CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
PROJECT(cdmft_ed Fortran)
SET(VERSION 0.15)

MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)


# CMakeGraphVizOptions.cmake
set(GRAPHVIZ_CUSTOM_TARGETS TRUE)


# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Check that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
OPTION(USE_MPI "Set MPI pre-compilation flag" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
OPTION(WITH_BLAS_LAPACK "Use internal Blas/Lapack" OFF)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)



############################################################
# Look for MPI if needed
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MpiConfig.cmake)






############################################################
# Look for required libraries: scifor, dmft_tools
# (includes Blas/Lapack/Scalapack/Quadpack/Minpack/FFTpack/Arpack/...) 
############################################################
# Find Scifor
FIND_PACKAGE(PkgConfig)
PKG_SEARCH_MODULE(SCIFOR REQUIRED scifor)
IF(${SCIFOR_FOUND})
  MESSAGE(STATUS "SCIFOR found at: ${SCIFOR_PREFIX}")
ELSE()
  MESSAGE(FATAL_ERROR "
Can not find SciFortran in this system. Check SciFor installation. Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_scifor_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})


# Find Dmft_tools
PKG_SEARCH_MODULE(DMFT_TOOLS REQUIRED dmft_tools)
IF(${DMFT_TOOLS_FOUND})
  MESSAGE(STATUS "DMFT_TOOLS found at: ${DMFT_TOOLS_PREFIX}")
ELSE()
  MESSAGE( FATAL_ERROR "
Can not find DMFT_TOOLS in this system. Check DMFT_TOOLS installation. Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_dmft_tools_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${DMFT_TOOLS_INCLUDE_DIRS})



############################################################
# Configure Fortran Compiler Main options
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/FortranCompilerConfig.cmake)





############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/GitConfig.cmake)
# IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
#   FIND_PACKAGE(Git)
#   IF(GIT_FOUND)
#     EXECUTE_PROCESS(
#       COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
#       WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
#       OUTPUT_VARIABLE GIT_SHA1
#       ERROR_QUIET
#       OUTPUT_STRIP_TRAILING_WHITESPACE)
#     MESSAGE( STATUS "Git version: ${GIT_SHA1}" )
#     EXECUTE_PROCESS(
#       COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
#       WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
#       OUTPUT_VARIABLE GIT_BRANCH
#       ERROR_QUIET
#       OUTPUT_STRIP_TRAILING_WHITESPACE)
#     SET(GIT_BRANCH_ "_${GIT_BRANCH}")
#     MESSAGE( STATUS "Git branch: ${GIT_BRANCH}" )
#     IF( "${GIT_BRANCH}" MATCHES "master")
#       UNSET(GIT_BRANCH)
#       UNSET(GIT_BRANCH_)
#     ENDIF()
#   ELSE(GIT_FOUND)
#     SET(GIT_SHA1 0)
#   ENDIF(GIT_FOUND)
# ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
# SET(REVISION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/revision.inc)
# FILE(WRITE  ${REVISION_FILE} "character(len=41),parameter,public :: revision = \"${GIT_SHA1}\"\n")


############################################################
# Define the actual files and folders that make up the build
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MainConfig.cmake)




############################################################
############################################################
############################################################
# Define some directories
SET(LIB_SRC ${CMAKE_SOURCE_DIR})
SET(DRIVERS ${CMAKE_SOURCE_DIR}/drivers)
SET(LIB_TMP_INC ${CMAKE_BINARY_DIR}/include)

#define PREFIX
SET(PREFIX  "$ENV{HOME}/.bin" CACHE PATH "Prefix prepended to install directories")
SET(CMAKE_INSTALL_PREFIX "${PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})

# define the executable to compiler
SET(EXE cdn_hm_2dsquare CACHE STRING "Driver code to compile")
# SET(CMAKE_INSTALL_EXE "${EXE}" CACHE INTERNAL "Driver code to compile" FORCE)


#INIT:
CONFIGURE_FILE(${LIB_SRC}/revision.in ${LIB_SRC}/ED_VERSION.f90 @ONLY)
FILE(GLOB ED_INIT_SRC
  ED_VERSION.f90
  ED_INPUT_VARS.f90
  ED_SPARSE_MATRIX.f90
  ED_SPARSE_MAP.f90
  ED_VARS_GLOBAL.f90
  ED_AUX_FUNX.f90
  ED_SECTOR.f90
  ED_SETUP.f90
  ED_EIGENSPACE.f90
)
ADD_LIBRARY(ED_INIT_LIB OBJECT  ${ED_INIT_SRC})


#BATH
ADD_SUBDIRECTORY(${LIB_SRC}/ED_BATH)
ADD_DEPENDENCIES(ED_BATH_LIB
  ED_INIT_LIB
)

#NORMAL
ADD_SUBDIRECTORY(${LIB_SRC}/ED_NORMAL)
ADD_DEPENDENCIES(ED_NORMAL_LIB
  ED_INIT_LIB
  ED_BATH_LIB
)


# #SUPERC
# ADD_SUBDIRECTORY(${LIB_SRC}/ED_SUPERC)
# ADD_DEPENDENCIES(ED_SUPERC_LIB
#   ED_INIT_LIB
#   ED_BATH_LIB
# )


#FIT
ADD_SUBDIRECTORY(${LIB_SRC}/ED_FIT)
ADD_DEPENDENCIES(ED_FIT_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  # ED_SUPERC_LIB
)


#HAMILTONIAN WRAPPER
FILE(GLOB ED_HAMILTONIAN_SRC
  ED_HAMILTONIAN.f90
)
ADD_LIBRARY(ED_HAMILTONIAN_LIB OBJECT  ${ED_HAMILTONIAN_SRC})
ADD_DEPENDENCIES(ED_HAMILTONIAN_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  # ED_SUPERC_LIB
)


#DIAG WRAPPER
FILE(GLOB ED_DIAG_SRC
  ED_DIAG.f90
)
ADD_LIBRARY(ED_DIAG_LIB OBJECT  ${ED_DIAG_SRC})
ADD_DEPENDENCIES(ED_DIAG_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  # ED_SUPERC_LIB
  ED_HAMILTONIAN_LIB
)


#GF WRAPPER
FILE(GLOB ED_GF_SRC
  ED_GREENS_FUNCTIONS.f90
)
ADD_LIBRARY(ED_GF_LIB OBJECT  ${ED_GF_SRC})
ADD_DEPENDENCIES(ED_GF_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  # ED_SUPERC_LIB
  ED_HAMILTONIAN_LIB
  ED_DIAG_LIB
)


#OBSERVABLES WRAPPER
FILE(GLOB ED_OBSERVABLES_SRC
  ED_OBSERVABLES.f90  
)
ADD_LIBRARY(ED_OBSERVABLES_LIB OBJECT  ${ED_OBSERVABLES_SRC})
ADD_DEPENDENCIES(ED_OBSERVABLES_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  # ED_SUPERC_LIB
  ED_HAMILTONIAN_LIB
  ED_DIAG_LIB
)


#IO
ADD_SUBDIRECTORY(${LIB_SRC}/ED_IO)
ADD_DEPENDENCIES(ED_IO_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_GF_LIB
)



#FINAL GLOBAL WRAPPER
FILE(GLOB ED_GLOB_SRC
  ED_MAIN.f90
  CDMFT_ED.f90
)
ADD_LIBRARY(ED_GLOB_LIB OBJECT  ${ED_GLOB_SRC})


ADD_DEPENDENCIES(ED_GLOB_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  # ED_SUPERC_LIB
  ED_HAMILTONIAN_LIB
  ED_DIAG_LIB
  ED_GF_LIB
  ED_OBSERVABLES_LIB
  ED_IO_LIB
  ED_FIT_LIB
)



SET(CDMFT_ED_OBJS
  $<TARGET_OBJECTS:ED_INIT_LIB>
  $<TARGET_OBJECTS:ED_BATH_LIB>
  $<TARGET_OBJECTS:ED_NORMAL_LIB>
  # $<TARGET_OBJECTS:ED_SUPERC_LIB>
  $<TARGET_OBJECTS:ED_HAMILTONIAN_LIB>
  $<TARGET_OBJECTS:ED_DIAG_LIB>
  $<TARGET_OBJECTS:ED_GF_LIB>
  $<TARGET_OBJECTS:ED_OBSERVABLES_LIB>
  $<TARGET_OBJECTS:ED_IO_LIB>
  $<TARGET_OBJECTS:ED_FIT_LIB>
  $<TARGET_OBJECTS:ED_GLOB_LIB>)



# SET(CDMFT_ED_OBJS ${CDMFT_ED_LOCAL_SRC})

ADD_EXECUTABLE(cdmft_ed ${DRIVERS}/${EXE}.f90 ${CDMFT_ED_OBJS})


target_LINK_LIBRARIES(cdmft_ed ${DMFT_TOOLS_LDFLAGS} ${SCIFOR_LDFLAGS})
SET_TARGET_PROPERTIES(cdmft_ed PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}" OUTPUT_NAME "${EXE}_${GIT_BRANCH}")
MESSAGE(STATUS "Executable created at: ${Yellow}${CMAKE_INSTALL_PREFIX}/${EXE}_${GIT_BRANCH}${ColourReset}")

ADD_CUSTOM_TARGET(
  completion
  COMMAND scifor_completion.sh ${DRIVERS}/${EXE}.f90
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
