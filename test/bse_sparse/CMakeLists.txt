set(CONFIG_FILES sparse_fit.ini.ref sparse_fit.ini.Lambda100)
foreach(config_file ${CONFIG_FILES})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${config_file}.in ${CMAKE_CURRENT_BINARY_DIR}/${config_file} @ONLY)
endforeach()

add_python_test(bse_sparse)
