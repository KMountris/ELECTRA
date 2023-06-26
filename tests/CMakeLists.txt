#----------------------------------------------

# Add definitions for compilation in Apple Architecture.
IF (APPLE)
    ADD_DEFINITIONS(-DGTEST_USE_OWN_TR1_TUPLE)
    ADD_DEFINITIONS(-D__GLIBCXX__)
ENDIF (APPLE)

# Find Google Test dependency.
FIND_PACKAGE(GTest  REQUIRED)
INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
