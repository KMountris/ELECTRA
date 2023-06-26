set(CUDA_ARCH "" CACHE STRING "Target CUDA Architectures multiple are allowed")

# # We use *.cu.cc as the default as most tool do not understand cu as CUDA.
# file(GLOB_RECURSE source_list "*.cu.cc")
# foreach(child ${source_list})
#   set_source_files_properties(${child} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
# endforeach()


if(CUDA_FOUND)
  message(STATUS "Found CUDA Toolkit v.${CUDA_VERSION}")

  if(CUDA_VERSION_MAJOR LESS 8)
    message(FATAL_ERROR "${PROJECT_NAME} requires CUDA Toolkit v.8.x or more recent. Please disable ${PROJECT_NAME}_WITH_CUDA option.")
  endif()

  message(STATUS "${PROJECT_NAME} build with CUDA support")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_CUDA ")
  include_directories(${CUDA_INCLUDE_DIRS})

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 --expt-relaxed-constexpr -DWITH_CUDA ")

  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}")
  endif()

  if(CUDA_ARCH STREQUAL "")
    # good defaults for CUDA Toolkit 8.x
    if(CUDA_VERSION_MAJOR MATCHES 8)
      set(CUDA_ARCH "35 37 52 60")
    endif()

    # good defaults for CUDA Toolkit 9.x
    if(CUDA_VERSION_MAJOR MATCHES 9)
      set(CUDA_ARCH "35 52 60 70")
    endif()

    # good defaults for CUDA Toolkit 10.x
    if(CUDA_VERSION_MAJOR MATCHES 10)
      set(CUDA_ARCH "35 52 60 70")
    endif()
  endif()

  # str replace ' ' with ;
  STRING(REGEX REPLACE " " ";" CUDA_ARCH ${CUDA_ARCH})

  # set the compiler flags for each NV target
  foreach(target ${CUDA_ARCH})
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode=arch=compute_${target},code=\\\"sm_${target},compute_${target}\\\")
  endforeach()

endif()