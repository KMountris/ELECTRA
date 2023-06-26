#--------------------------------------------------------------
set(LIB_NAME ${PROJECT_NAME})
string(TOLOWER ${PROJECT_NAME} LC_PROJECT_NAME)

# The headers collection of the IMP library.
set(HEADERS ${CMAKE_SOURCE_DIR}/src/${PROJECT_NAME}/${PROJECT_NAME})

# Subdirectories of ELECTRA.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/conditions)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/conduction_system)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/electrophysiology)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/exporters)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/fibers)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/materials)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/physics)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/post_process)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/engine/utilities)

set(${LIB_NAME}_LIBS
    ${LC_PROJECT_NAME}Conditions
    ${LC_PROJECT_NAME}ConductionSystem
    ${LC_PROJECT_NAME}Electrophysiology
    ${LC_PROJECT_NAME}Exporters
    ${LC_PROJECT_NAME}Fibers
    ${LC_PROJECT_NAME}Materials
    ${LC_PROJECT_NAME}Physics
    ${LC_PROJECT_NAME}PostProcess
    ${LC_PROJECT_NAME}Utilities
)


#-------- Build library --------
add_library(${LIB_NAME} "")
target_sources(${LIB_NAME} PRIVATE ${HEADERS})

target_link_libraries(${LIB_NAME} PUBLIC ${${LIB_NAME}_LIBS})
target_link_libraries(${LIB_NAME} PUBLIC IMP::IMP CLOUDEA::CLOUDEA)
target_link_libraries(${LIB_NAME} PUBLIC Eigen3::Eigen armadillo Boost::boost tinyxml2 CGAL::CGAL CGAL::CGAL_Core)

if(${ELECTRA_WITH_MPI})
    target_link_libraries(${LIB_NAME} PUBLIC MPI::MPI_CXX)
endif()

target_include_directories(${LIB_NAME} PUBLIC
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

set_target_properties(${LIB_NAME}
    PROPERTIES
        LINKER_LANGUAGE CXX
        POSITION_INDEPENDENT_CODE 1
        DEBUG_POSTFIX "_d"
        PUBLIC_HEADER "${HEADERS}"
        $<$<BUILD_SHARED_LIBS>:CXX_VISIBILITY_PRESET hidden>
        $<$<BUILD_SHARED_LIBS>:VISIBILITY_INLINES_HIDDEN 1>
        $<$<BUILD_SHARED_LIBS>:SOVERSION ${PROJECT_VERSION}>
        $<$<BUILD_SHARED_LIBS>:OUTPUT_NAME ${PROJECT_NAME}>
        $<$<BUILD_SHARED_LIBS>:MACOSX_RPATH ON>
        $<$<NOT:$<BUILD_SHARED_LIBS>>:ARCHIVE_OUTPUT_NAME ${PROJECT_NAME}>
)
#--------------------------------------

install(TARGETS
        ${LIB_NAME}
    EXPORT
        ${LIB_NAME}Targets
    RUNTIME
        DESTINATION ${INSTALL_BINDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
    ARCHIVE
        DESTINATION ${INSTALL_LIBDIR}
        COMPONENT ${PROJECT_NAME}_Development
    PUBLIC_HEADER
        DESTINATION ${INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        COMPONENT ${PROJECT_NAME}_Development
)

install(EXPORT ${LIB_NAME}Targets
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION ${INSTALL_CMAKEDIR}
    COMPONENT ${PROJECT_NAME}_Development
)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION
        ${INSTALL_CMAKEDIR}
    PATH_VARS
        INSTALL_INCLUDEDIR
        INSTALL_LIBDIR
)

write_basic_package_version_file(${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION
        ${PROJECT_VERSION}
    COMPATIBILITY
        SameMajorVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION
        ${INSTALL_CMAKEDIR}
)
