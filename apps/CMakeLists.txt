####
# Build ELECTRA apps

# Options for building ELECTRA apps
option(BUILD_APPS_WITH_EXTERNAL_JSON "Use an external build of the NLOHMANN JSON library" OFF)
option(BUILD_APPS_WITH_EXTERNAL_TERMCOLOR "Use an external build of the TERMCOLOR library" OFF)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tools)

# Third party libraries for building ELECTRA apps
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/nlohmann ${CMAKE_BINARY_DIR}/thirdparty/nlohmann)
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/termcolor ${CMAKE_BINARY_DIR}/thirdparty/termcolor)


# Build ElectraSim app
add_executable(ElectraSim electra_sim.cpp)
target_compile_definitions(ElectraSim PRIVATE ELECTRA_VERSION="${PROJECT_VERSION}")
target_link_libraries(ElectraSim PRIVATE ${PROJECT_NAME} Tools nlohmann_json::nlohmann_json termcolor::termcolor pthread)
target_include_directories(ElectraSim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)

# Build ElectraPre app
add_executable(ElectraPre electra_pre.cpp)
target_compile_definitions(ElectraPre PRIVATE ELECTRA_VERSION="${PROJECT_VERSION}")
target_link_libraries(ElectraPre PRIVATE ${PROJECT_NAME} Tools nlohmann_json::nlohmann_json termcolor::termcolor pthread)
target_include_directories(ElectraPre PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)

# Build ElectraPost app
add_executable(ElectraPost electra_post.cpp)
target_compile_definitions(ElectraPost PRIVATE ELECTRA_VERSION="${PROJECT_VERSION}")
target_link_libraries(ElectraPost PRIVATE ${PROJECT_NAME} Tools nlohmann_json::nlohmann_json termcolor::termcolor pthread)
target_include_directories(ElectraPost PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)

# Build single cell dev app
add_executable(SingleCelldev  electra_cell_dev.cpp)
target_link_libraries(SingleCelldev PRIVATE ${PROJECT_NAME})
target_include_directories(SingleCelldev PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)

# Install apps in targets
install(TARGETS
    ElectraSim
    ElectraPre
    ElectraPost
    SingleCelldev
    RUNTIME
        DESTINATION ${INSTALL_BINDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
)