# Build TERMCOLOR thirdparty library

if(BUILD_APPS_WITH_EXTERNAL_TERMCOLOR)
    find_package(termcolor REQUIRED)
else()
    include(FetchContent)

    FetchContent_Declare(termcolor
        GIT_REPOSITORY https://github.com/ikalnytskyi/termcolor
        GIT_TAG v2.0.0
    )

    FetchContent_GetProperties(termcolor)
    if(NOT termcolor_POPULATED)
        FetchContent_Populate(termcolor)
        FetchContent_MakeAvailable(termcolor)
        add_subdirectory(${termcolor_SOURCE_DIR} ${termcolor_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
endif() 