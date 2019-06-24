FILE(COPY sparse_fit.ini.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY sparse_fit.ini.Lambda10 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_python_test(bse_sparse)
