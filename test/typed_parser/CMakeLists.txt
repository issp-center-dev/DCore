FILE(COPY parser.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

triqs_add_python_test(typed_parser)
