#-----------------------------------------

# Check if Doxygen is installed in the OS.
FIND_PACKAGE(Doxygen)
# # Set input and output files.
# SET(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/electra_doxyfile.in)
# SET(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/electra_doxyfile)
# # Request to configure the Doxygen documentation file.
# CONFIGURE_FILE(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
# # Add custom target "doc" for documentation compilation.
# ADD_CUSTOM_TARGET(doc COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#                       COMMENT "Generating the ELECTRA Doxygen documentation..."
#                       VERBATIM )
# Find all the public headers
GET_TARGET_PROPERTY(${PROJECT_NAME}_PUBLIC_HEADER_DIR ${PROJECT_NAME} INTERFACE_INCLUDE_DIRECTORIES)

FILE(GLOB_RECURSE ${PROJECT_NAME}_PUBLIC_HEADERS ${${PROJECT_NAME}_PUBLIC_HEADER_DIR}/*.hpp)

SET(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/electra/docs)
SET(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
SET(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
SET(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in)
SET(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/doxyfile)

#Replace variables inside @@ with the current values
CONFIGURE_FILE(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# FILE(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
ADD_CUSTOM_COMMAND(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${${PROJECT_NAME}_PUBLIC_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating the ${PROJECT_NAME} Doxygen documentation..."
)
ADD_CUSTOM_TARGET(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})


FIND_PACKAGE(Sphinx)

SET(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
SET(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
SET(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND ${SPHINX_EXECUTABLE} -b html
                     # Tell Breathe where to find the Doxygen output
                     -Dbreathe_projects.${PROJECT_NAME}=${DOXYGEN_OUTPUT_DIR}/xml
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                   # Other docs files you want to track should go here (or in some variable)
                   ${SPHINX_SOURCE}/index.rst
                   ${DOXYGEN_INDEX_FILE}
                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                   COMMENT "Generating the ${PROJECT_NAME} Sphinx documentation..."
)

# Nice named target so we can run the job easily
ADD_CUSTOM_TARGET(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})