
cmake_minimum_required(VERSION 3.13)

project(ProgrammingUtilities)

add_subdirectory(src)
add_subdirectory(include)
#add_subdirectory(pybindings)

add_subdirectory(tests)

target_include_directories(extramath PUBLIC include)
