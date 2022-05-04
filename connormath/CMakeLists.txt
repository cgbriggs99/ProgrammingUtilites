
project(ProgrammingUtilities)

add_subdirectory(roots)
add_subdirectory(utils)

add_subdirectory(tests)

target_include_directories(connormath PUBLIC include)
