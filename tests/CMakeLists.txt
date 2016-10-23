
file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include_directories(..)
set(CMAKE_CXX_FLAGS "-Wno-null-dereference -g -O0 -std=c++14")

foreach(_src_file_path ${TEST_SOURCES})
    message(STATUS ${_src_file_path})
    string(REPLACE ".cpp" "" test_case ${_src_file_path})
    string(REGEX MATCH "([^/]*)$" test_name ${test_case})

    set(EXECUTABLE_NAME "${test_name}.out")
    add_executable(${EXECUTABLE_NAME} ${_src_file_path})

    set_target_properties(${EXECUTABLE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
endforeach()
