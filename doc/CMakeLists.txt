find_package(Sphinx REQUIRED)
message(STATUS "Found sphinx at ${SPHINX_EXECUTABLE}")

find_package (Git)
if (GIT_FOUND)
      message("git found: ${GIT_EXECUTABLE} in version ${GIT_VERSION_STRING}")
endif(GIT_FOUND)

# Get branch name and substitute to GIT_BRANCH_NAME
execute_process(
    COMMAND ${GIT_EXECUTABLE} 
    --git-dir=${CMAKE_SOURCE_DIR}/.git --work-tree=${CMAKE_SOURCE_DIR} rev-parse --abbrev-ref HEAD
    RESULT_VARIABLE result
    OUTPUT_VARIABLE GIT_BRANCH_NAME)
string(STRIP ${GIT_BRANCH_NAME} GIT_BRANCH_NAME)
message(STATUS "git branch name is ${GIT_BRANCH_NAME}")

# generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

# all rst files of the documentation
file(GLOB_RECURSE doc_sources *.rst EXCLUDE *.swp)

# Make a local copy of doc/* as doc_source
file(GLOB_RECURSE all_doc_files RELATIVE ${CMAKE_SOURCE_DIR}/doc * EXCLUDE *.swp)
add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/doc_source.tar DEPENDS ${all_doc_files} COMMAND ${CMAKE_COMMAND} -E tar cf ${CMAKE_BINARY_DIR}/doc_source.tar ${all_doc_files} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc)
add_custom_target(make_doc_local_copy_dir ALL COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/doc_source)
add_custom_target(copy_tar_doc ALL DEPENDS ${CMAKE_BINARY_DIR}/doc_source.tar make_doc_local_copy_dir COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_BINARY_DIR}/doc_source.tar WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc_source)
add_custom_target(make_option_tables ALL DEPENDS copy_tar_doc COMMAND ${PYTHON_INTERPRETER} ${CMAKE_SOURCE_DIR}/python/option_tables.py ${CMAKE_BINARY_DIR}/doc_source/reference)
add_custom_target(make_doc_local_copy ALL DEPENDS make_option_tables)

# create documentation target
set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${doc_sources} 
    COMMAND ${SPHINX_EXECUTABLE} ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -b html ${CMAKE_BINARY_DIR}/doc_source html)

#option(DOC_SPELLING "Activate spellcheck for documentations" OFF)
#if (DOC_SPELLING)
    #set(sphinx_spelling_out ${CMAKE_CURRENT_BINARY_DIR}/spelling/output.txt)
    #add_custom_command(OUTPUT ${sphinx_spelling_out} DEPENDS ${doc_sources} spelling_wordlist.txt
                   #COMMAND ${CMAKE_BINARY_DIR}/build_pytriqs ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -b spelling ${CMAKE_CURRENT_SOURCE_DIR} spelling)
#endif ()
#add_custom_target(doc_sphinx ALL DEPENDS make_doc_local_copy ${sphinx_top} ${sphinx_spelling_out})
add_custom_target(doc_sphinx ALL DEPENDS make_doc_local_copy ${sphinx_top})

# install
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/dcore
 FILES_MATCHING 
 PATTERN "*.html" 
 PATTERN "*.png"
 PATTERN "*.js"
 PATTERN "_*"
 PATTERN "*.jpg" 
 PATTERN "*.gif" 
 PATTERN "*.xsl" 
 PATTERN "*.css"
 PATTERN "*.pdf"
 PATTERN "*.py"
 PATTERN "*.txt"
 PATTERN "*.inv"
 PATTERN "*.bib"
 PATTERN "*.ini"
 PATTERN "*.in"
 PATTERN "*.win"
 PATTERN "*.dat"
 )
