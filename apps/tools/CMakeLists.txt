#----------------------------------------------
set(LIB_NAME Tools)
set(LIB_FOLDER tools)

# Module source files.
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/config_electrophys.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_stimuli.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_units.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/parser.cpp
)

# Module header files.
SET(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/parser.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/parser.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_approximation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_approximation.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_conduct_sys.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_conduct_sys.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_electric.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_electric.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_electrophys.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_fibers.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_fibers.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_geo.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_geo.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_output.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_output.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_physics.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_physics.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_post_process.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_post_process.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_sim.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_sim.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_stimuli.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/config_units.hpp
)

#-------- Build library --------
add_library(${LIB_NAME} "")
target_sources(${LIB_NAME} PRIVATE ${SOURCES})
target_compile_definitions(${LIB_NAME} PUBLIC ELECTRA_VERSION="${PROJECT_VERSION}")
target_link_libraries(${LIB_NAME} PRIVATE ${PROJECT_NAME} nlohmann_json::nlohmann_json)

include(GenerateExportHeader)
generate_export_header(${LIB_NAME}
    BASE_NAME "${LIB_NAME}"
    EXPORT_MACRO_NAME "${LIB_NAME}_EXPORT"
    EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/apps/${LIB_FOLDER}/${LIB_NAME}Export.h"
    DEPRECATED_MACRO_NAME "${LIB_NAME}_DEPRECATED"
    NO_EXPORT_MACRO_NAME "${LIB_NAME}_NO_EXPORT"
    STATIC_DEFINE "${LIB_NAME}_STATIC_DEFINE"
    NO_DEPRECATED_MACRO_NAME "${LIB_NAME}_NO_DEPRECATED"
    DEFINE_NO_DEPRECATED
)

target_include_directories(${LIB_NAME} PUBLIC
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/apps>
)

set_target_properties(${LIB_NAME}
    PROPERTIES
        LINKER_LANGUAGE CXX
        POSITION_INDEPENDENT_CODE 1
        DEBUG_POSTFIX "_d"
        PUBLIC_HEADER "${HEADERS};${CMAKE_BINARY_DIR}/apps/${LIB_FOLDER}/${LIB_NAME}Export.h"
        $<$<BUILD_SHARED_LIBS>:CXX_VISIBILITY_PRESET hidden>
        $<$<BUILD_SHARED_LIBS>:VISIBILITY_INLINES_HIDDEN 1>
        $<$<BUILD_SHARED_LIBS>:SOVERSION ${PROJECT_VERSION}>
        $<$<BUILD_SHARED_LIBS>:OUTPUT_NAME ${LIB_NAME}>
        $<$<BUILD_SHARED_LIBS>:MACOSX_RPATH ON>
        $<$<NOT:$<BUILD_SHARED_LIBS>>:ARCHIVE_OUTPUT_NAME ${LIB_NAME}>
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
        DESTINATION ${INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${LIB_FOLDER}
        COMPONENT ${PROJECT_NAME}_Development
)

# install(EXPORT ${LIB_NAME}Targets
#     NAMESPACE "${PROJECT_NAME}::"
#     DESTINATION ${INSTALL_CMAKEDIR}
#     COMPONENT ${PROJECT_NAME}_Development
# )