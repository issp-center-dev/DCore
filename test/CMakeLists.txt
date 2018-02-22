# load triqs helper to set up tests
find_package(TriqsTest)

add_subdirectory(typed_parser)
add_subdirectory(openmx)
add_subdirectory(respack)
add_subdirectory(pre_preset)
add_subdirectory(pre_wannier)
add_subdirectory(pre_wannier_so)
add_subdirectory(main_chain_hubbardI)
add_subdirectory(check_chain_hubbardI)
add_subdirectory(post_chain_hubbardI)
