# where will the python end up in triqs?
set(python_destination pytriqs/applications/pydmft)

# site_customize for build
set(package_name "pytriqs.applications.pydmft")
configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# make a local pytriqs copy
triqs_prepare_local_pytriqs(${python_destination})

# add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
#install(FILES ${CMAKE_SOURCE_DIR}/bin/pydmft.py DESTINATION ${TRIQS_PATH}/bin)
