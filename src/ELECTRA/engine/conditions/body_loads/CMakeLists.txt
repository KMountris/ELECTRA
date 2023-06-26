#----------------------------------------------
set(LIB_NAME BodyLoads)
set(LIB_FOLDER conditions/body_loads)

# Module source files.
set(SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/stimulus.cpp 
)

# Module header files.
set(HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/body_loads.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/stimulus.hpp
)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(${LIB_NAME} OBJECT ${SOURCES} ${HEADERS})
target_link_libraries(${LIB_NAME} PUBLIC Eigen3::Eigen)

target_include_directories(${LIB_NAME} PUBLIC
    $<INSTALL_INTERFACE:${${PROJECT_NAME}_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

set_target_properties(${LIB_NAME}
    PROPERTIES
        LINKER_LANGUAGE CXX
        PUBLIC_HEADER "${HEADERS}"
)

install(TARGETS
    PUBLIC_HEADER
        DESTINATION ${${PROJECT_NAME}_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${LIB_FOLDER}
        COMPONENT ${PROJECT_NAME}_Development
)
