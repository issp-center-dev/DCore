# list all files to be installed...

set(FILES openmx2dcore.py respack2wan90.py)

foreach(fname IN LISTS FILES)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${fname}  DESTINATION  bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )
    #configure_file(${fname} ${CMAKE_BINARY_DIR}/python/dcore/${fname} @ONLY)
endforeach()
