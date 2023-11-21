cmake_minimum_required(VERSION 3.14)

project(FMHD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++)

file(GLOB_RECURSE src CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -msse4.2 -mavx2")
#add_compile_options(-O2 -fopenmp -msse4.2 -mavx)
add_library(minhash_lib ${src})
target_include_directories(minhash_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
set_target_properties(minhash_lib PROPERTIES OUTPUT_NAME "minhash")

add_executable(fmhd ${src})
target_include_directories(fmhd PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(fmhd-dev ${src})
target_include_directories(fmhd-dev PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_definitions(fmhd-dev PUBLIC "-DDEVELOPER_BUILD")

if (CYGWIN OR MINGW)
    target_link_libraries(fmhd pthread)
    target_link_libraries(fmhd -static-libgcc)
    target_link_libraries(fmhd -static-libstdc++)
    target_link_libraries(fmhd-dev pthread)
    target_link_libraries(fmhd-dev -static-libgcc)
    target_link_libraries(fmhd-dev -static-libstdc++)
else()
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(fmhd PUBLIC OpenMP::OpenMP_CXX )
        target_link_libraries(fmhd-dev PUBLIC OpenMP::OpenMP_CXX )
    endif()
endif()
