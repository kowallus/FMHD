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

add_executable(minhash_executable
        ${src}
)
target_include_directories(minhash_executable PUBLIC ${PROJECT_SOURCE_DIR}/include)

if (CYGWIN OR MINGW)
    target_link_libraries(minhash_executable pthread)
    target_link_libraries(minhash_executable -static-libgcc)
    target_link_libraries(minhash_executable -static-libstdc++)
else()
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(minhash_executable
                PUBLIC OpenMP::OpenMP_CXX
        )
    endif()
endif()
