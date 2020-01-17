set(TEST_FILES dmft.ini)
foreach(test_file ${TEST_FILES})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${test_file} ${CMAKE_CURRENT_BINARY_DIR}/${test_file} COPYONLY)
endforeach()

# COPY REFERENCE DATA
file(GLOB REF_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/ref ${CMAKE_CURRENT_SOURCE_DIR}/ref/*)
foreach(test_file ${REF_FILES})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ref/${test_file} ${CMAKE_CURRENT_BINARY_DIR}/ref/${test_file} COPYONLY)
endforeach()

add_python_test(chain_hubbardI_so)
