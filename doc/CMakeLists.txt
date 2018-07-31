# generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

# all rst files of the documentation
file(GLOB_RECURSE doc_sources *.rst)

# create documentation target
set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${doc_sources} py_copy
                   COMMAND ${CMAKE_BINARY_DIR}/build_pytriqs ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -b html ${CMAKE_CURRENT_SOURCE_DIR} html)
option(DOC_SPELLING "Activate spellcheck for documentations" OFF)
if (DOC_SPELLING)
    set(sphinx_spelling_out ${CMAKE_CURRENT_BINARY_DIR}/spelling/output.txt)
    add_custom_command(OUTPUT ${sphinx_spelling_out} DEPENDS ${doc_sources} py_copy spelling_wordlist.txt
                   COMMAND ${CMAKE_BINARY_DIR}/build_pytriqs ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -b spelling ${CMAKE_CURRENT_SOURCE_DIR} spelling)
endif ()
add_custom_target(doc_sphinx ALL DEPENDS ${sphinx_top} ${sphinx_spelling_out})

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
 )
