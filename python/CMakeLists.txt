# where will the python end up in triqs?
set(python_destination dcore)

# site_customize for build
set(package_name "dcore")
configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# Install path
install(DIRECTORY ${CMAKE_SOURCE_DIR}/python/ DESTINATION
${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/dcore FILES_MATCHING PATTERN *.py)

# Copy all Python files for ctest into a build directory
file(GLOB_RECURSE all_py_files RELATIVE ${CMAKE_SOURCE_DIR}/python *.py)
add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/py_source.tar
    COMMAND ${CMAKE_COMMAND} -E tar cf ${CMAKE_BINARY_DIR}/py_source.tar ${all_py_files} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/python/dcore
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR}/python/dcore tar xf ${CMAKE_BINARY_DIR}/py_source.tar
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/py_source.tar
    DEPENDS ${all_py_files}
)

add_custom_target(copy_tar_py ALL DEPENDS ${CMAKE_BINARY_DIR}/py_source.tar)

# add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${DCORE_SITE_PACKAGES}/dcore)
