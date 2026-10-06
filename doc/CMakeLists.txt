# ---- Build documentation ---------------------------------------------------
# Find and configure Doxygen
find_package(Doxygen 1.8.15 REQUIRED)

set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/doc/Doxygen.in)
set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/doc/Doxygen)
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

if (RODIN_USE_MCSS)
  set(DOXYGEN_MCSS_IN ${CMAKE_SOURCE_DIR}/doc/Doxygen.mcss.in)
  set(DOXYGEN_MCSS_OUT ${CMAKE_BINARY_DIR}/doc/Doxygen.mcss)
  configure_file(${DOXYGEN_MCSS_IN} ${DOXYGEN_MCSS_OUT} @ONLY)

  set(MCSS_CONF_IN ${CMAKE_SOURCE_DIR}/doc/conf.py.in)
  set(MCSS_CONF_OUT ${CMAKE_BINARY_DIR}/doc/conf.py)
  configure_file(${MCSS_CONF_IN} ${MCSS_CONF_OUT} @ONLY)

  file(COPY favicon.ico DESTINATION ${CMAKE_BINARY_DIR}/doc)

  rodin_add_mcss("${MCSS_CONF_OUT}")
else()
  add_custom_target(RodinDoxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generate API documentation with Doxygen"
    VERBATIM
    )
endif()
