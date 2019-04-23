set(EXEC_NAMES dcore_pre dcore dcore_post dcore_check dcore_bse)

foreach(exec IN LISTS EXEC_NAMES)
    message(STATUS ${exec})
    configure_file(${exec}.bash.in ${exec})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${exec} DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endforeach()
