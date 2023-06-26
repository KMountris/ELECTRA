#--------------------------------------------------------------
cmake_minimum_required(VERSION 3.14.0 FATAL_ERROR)

project(ELECTRA
  VERSION 0.5.3
  LANGUAGES CXX
)

# Set project's CMake modules path
list(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
  ${CMAKE_SOURCE_DIR}/cmake/modules/
)

#------- Compilation Options -------#
option(BUILD_SHARED_LIBS "Build ${PROJECT_NAME} as shared library" ON)
option(BUILD_APPS "Build all ${PROJECT_NAME} applications" ON)
option(BUILD_TESTS "Build all ${PROJECT_NAME} tests" OFF)
option(BUILD_DOC "Build ${PROJECT_NAME} documentation" ON)
option(${PROJECT_NAME}_WITH_MPI "Build ${PROJECT_NAME} with MPI parallelization backend" OFF)
option(${PROJECT_NAME}_WITH_CUDA "Build ${PROJECT_NAME} with CUDA GPU parallelization backend" OFF)
#-----------------------------------#


# Set the install prefix under the supported architectures.
if(MSVC)
    set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/ELECTRA" CACHE PATH "CMake prefix")
else()
    set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "CMake prefix")
endif()

# General set up
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type to Release mode.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "ELECTRA Build type. Options: Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
message(STATUS "${PROJECT_NAME} build type set to ${CMAKE_BUILD_TYPE}")


# Find dependencies of ELECTRA.
# -------------------------

# Find IMP.
find_package(IMP REQUIRED 0.0.1)
include_directories(${IMP_INCLUDE_DIRS})

# Find CLOUDEA.
find_package(CLOUDEA REQUIRED)
include_directories(${CLOUDEA_INCLUDE_DIRS})
link_directories(${CLOUDEA_LIBRARY_DIR})

# Find Boost.
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find Armadillo
find_package(Armadillo REQUIRED NO_MODULE)

# Find TinyXML2
find_package(TinyXML2 REQUIRED)
include_directories(${TINYXML2_INCLUDE_DIRS})

# Find CGAL.
find_package(CGAL REQUIRED COMPONENTS Core)
include_directories(${CGAL_INCLUDE_DIRS})

# Find CUDA
if (${PROJECT_NAME}_WITH_CUDA)
  find_package(CUDA REQUIRED)
  include(CudaHelper)
endif()

# Find MPI
# if(${ELECTRA_USE_MPI})
#     find_package(MPI REQUIRED)
# endif()
# -------------------------

# Set build and installation directories
include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# Set the install prefix for the system's architecture.
if (WIN32)
  if (CMAKE_GENERATOR_PLATFORM EQUAL x64)
    set(INSTALL_PREFIX "C:/Program Files/${PROJECT_NAME}")
  else()
    set(INSTALL_PREFIX "C:/Program Files (x86)/${PROJECT_NAME}")
  endif(CMAKE_GENERATOR_PLATFORM EQUAL x64)
else()
  set(INSTALL_PREFIX "/usr/local")
endif (WIN32)
set(CMAKE_INSTALL_PREFIX INSTALL_PREFIX CACHE PATH "Installation directory prefix for ${PROJECT_NAME}")

message(STATUS "${PROJECT_NAME} will be installed in ${CMAKE_INSTALL_PREFIX}")

# Allow overriding the installation directories
set(VERSIONED_PROJECT_NAME ${PROJECT_NAME}-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR}/${VERSIONED_PROJECT_NAME} CACHE PATH "Installation directory for ${PROJECT_NAME} libraries")
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR}/${VERSIONED_PROJECT_NAME} CACHE PATH "Installation directory for ${PROJECT_NAME} applications")
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/${VERSIONED_PROJECT_NAME} CACHE PATH "Installation directory for ${PROJECT_NAME} header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKEDIR CMake/${VERSIONED_PROJECT_NAME})
else()
  set(DEF_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${VERSIONED_PROJECT_NAME})
endif()
set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for ${PROJECT_NAME} CMake files")

# Make relative paths, if any, absolute.
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()
#------------------------------------------------------#

#--- Report Installation directories to user
foreach(p LIB BIN INCLUDE CMAKE)
  file(TO_NATIVE_PATH ${INSTALL_${p}DIR} _path )
  message(STATUS "Installing ${p} components in ${_path}")
  unset(_path)
endforeach()

# # Setup execution version definition
# IF(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
#     ADD_DEFINITIONS(-DEXECUTION_MODEL="MPI_CPU")
# ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES Windows)
#     ADD_DEFINITIONS(-DCMAKEMACROSAMPLE="SINGLE_CPU")
# ENDIF()

# To choose system
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(STATUS "Configuring on/for Linux")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  message(STATUS "Configuring on/for macOS")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  message(STATUS "Configuring on/for Windows")
else()
  message(STATUS "Configuring on/for ${CMAKE_SYSTEM_NAME}")
endif()

# Set compile options for supported compilers
if(MSVC)
    set(ELECTRA_COMPILE_OPTIONS "/W4" "$<$<CONFIG:RELEASE>:/O2>" "-msse2")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(ELECTRA_COMPILE_OPTIONS "-Wall" "-Wextra" "-Wdouble-promotion" "-ansi" "-pedantic" "-fexceptions"
                            "-fsignaling-nans" "-fPIC" "-stdlib=libc++" "$<$<CONFIG:RELEASE>:-O2>" "-ftree-vectorize" "-msse2" "-frounding-math" "-ftrapping-math" "-ffast-math")
else()
    set(ELECTRA_COMPILE_OPTIONS "-Wall" "-Wextra" "-Wdouble-promotion" "-Wc++14-compat" "-Werror=format-security"
                                "-pedantic" "-fsignaling-nans" "-m64" "-funroll-loops"
                                "-march=native" "-mfpmath=sse" "-mieee-fp"
                                "-fexceptions" "-fPIC" "-shared" "$<$<CONFIG:RELEASE>:-O2>" "-ftree-vectorize" "-msse2" "-frounding-math" "-ftrapping-math" "-ffast-math" "-pthread")
endif()

add_compile_options(${ELECTRA_COMPILE_OPTIONS})

# Compile the ELECTRA library components.
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

# Uninstall ELECTRA custom command.
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}_uninstall.cmake.in
  ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_uninstall.cmake  IMMEDIATE  @ONLY
)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_uninstall.cmake)

# Build ELECTRA applications.
if(BUILD_APPS)
  add_subdirectory(${CMAKE_SOURCE_DIR}/apps)
endif()

# Add ELECTRA documentation in the compilation.
if(BUILD_DOC)
  add_subdirectory(${CMAKE_SOURCE_DIR}/docs)
endif()

# Add ELECTRA tests in the compilation.
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(${CMAKE_SOURCE_DIR}/tests)
endif()


# Report details for ELECTRA compilation and installation.
message("=================================================================================================")
message(STATUS "$BUILD_PATH=${CMAKE_BINARY_DIR}")
if(BUILD_DOC)
message(STATUS "$BUILD_DOC_PATH=${CMAKE_BINARY_DIR}/doc")
endif()
message(STATUS "$INSTALL_PATH=${CMAKE_INSTALL_PREFIX}")
message(STATUS "$COMPILATION OPTIONS:")
message(STATUS "                         make            |  build ${PROJECT_NAME} libary in ${PROJECT_NAME}_BUILD_PATH")
if(BUILD_DOC)
message(STATUS "                         make doc        |  build ${PROJECT_NAME} documentation in ${PROJECT_NAME}_BUILD_DOC_PATH")
endif()
message(STATUS "                         make install    |  install ${PROJECT_NAME} library in ${PROJECT_NAME}_INSTALL_PATH")
message(STATUS "                         make uninstall  |  uninstall ${PROJECT_NAME} library from ${PROJECT_NAME}_INSTALL_PATH")
if(BUILD_TESTS)
message(STATUS "                         make test       |  execute ${PROJECT_NAME} unitary tests")
endif()
message("=================================================================================================")
