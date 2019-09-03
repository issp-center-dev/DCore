FILE(COPY external.h5.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY external_ref.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_python_test(pre_external)
