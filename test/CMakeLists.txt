# load triqs helper to set up tests
find_package(TriqsTest)

# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# Copy other files
FILE(COPY parser.in DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY nis_hr.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#FILE(COPY SrVO3.pmat SrVO3.struct SrVO3.outputs SrVO3.oubwin SrVO3.ctqmcout SrVO3.symqmc SrVO3.sympar SrVO3.parproj hk_convert_hamiltonian.hk LaVO3-Pnma_hr.dat LaVO3-Pnma.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

triqs_add_python_test(parser_test)
triqs_add_python_test(dmft_solver)
triqs_add_python_test(standard_test)
triqs_add_python_test(wannier_test)
triqs_add_python_test(U_matrix_cubic_test)
#triqs_add_python_test(w90_convert)
#triqs_add_python_test(sumkdft_basic)
#triqs_add_python_test(srvo3_Gloc)
#triqs_add_python_test(srvo3_transp)
#triqs_add_python_test(sigma_from_file)
#triqs_add_python_test(blockstructure)
