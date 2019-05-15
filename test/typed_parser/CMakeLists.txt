FILE(COPY parser.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_python_test(typed_parser)
