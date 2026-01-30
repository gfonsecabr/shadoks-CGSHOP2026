cmake_minimum_required(VERSION 3.12)
set(projectName shadoks)
message("cmake for ${projectName}")
project(${projectName})

enable_testing()

# List of local libraries
set(localLibs
    MaLib
    cadical
    EvalMaxSAT
    glucose
)

# List of external libraries
set(externLibs
    z
)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
set(CMAKE_EXE_LINKER_FLAGS "-static")

include_directories(lib/eigen)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add local libraries
FOREACH(localLib ${localLibs})
    if(NOT TARGET ${localLib})
        add_subdirectory(lib/${localLib})
    endif()
ENDFOREACH(localLib)

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIR})

include_directories(lib)

add_executable(solver solver.cpp)
add_executable(improver improver.cpp)
add_executable(exact exact.cpp)
add_executable(makelb makelb.cpp)

target_link_libraries(solver
    ${localLibs}
    ${externLibs}
)

target_link_libraries(improver
    ${localLibs}
    ${externLibs}
)

target_link_libraries(exact
    ${localLibs}
    ${externLibs}
)

# install(TARGETS ${projectName}
#     DESTINATION ${CMAKE_INSTALL_BINDIR}
# )


