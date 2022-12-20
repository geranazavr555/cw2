cmake_minimum_required(VERSION 3.22)
project(cw2)

set(CMAKE_CXX_STANDARD 20)

#add_compile_options("-fopencilk" "-O0" "-g" "-fsanitize=address" "-fno-omit-frame-pointer")
#add_link_options("-fopencilk" "-fsanitize=address")
#
#add_compile_options("-fopencilk" "-O0" "-g")
#add_link_options("-fopencilk")

#add_compile_options("-fopencilk" "-O3" "-g")
#add_link_options("-fopencilk")
#
#add_executable(cw2 main.cpp)

find_package(Threads REQUIRED)
find_package(Parlay REQUIRED)

#add_compile_options("-fopencilk" "-DPARLAY_OPENCILK" "-O3")
#add_link_options("-fopencilk")

add_compile_options("-O3")

#
#add_compile_options("-O0" "-g" "-fsanitize=address" "-fno-omit-frame-pointer")
#add_link_options("-fsanitize=address")

add_executable(cw2 main2.cpp)
target_link_libraries(cw2 PUBLIC Parlay::parlay)
