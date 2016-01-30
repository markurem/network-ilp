# find python
find_package(PythonInterp)
if (NOT PYTHON_EXECUTABLE)
    set(BUILD_PYTHON_BINDINGS OFF)
	message(STATUS "No python executable found.")
endif()

configure_file(setup.py.tpl setup.py)

install(DIRECTORY netopt DESTINATION share/netopt/python)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup.py 
		DESTINATION share/netopt/python)

if (PYTHON_EXECUTABLE)
    install(CODE "execute_process(
        	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py install
        	WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}\")")
endif()
