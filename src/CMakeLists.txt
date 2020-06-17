cmake_minimum_required(VERSION 3.8)
add_executable(compass "main.cpp" "data.hpp" "fit.hpp" "circle.hpp")

target_link_libraries(compass
    PUBLIC
        Eigen3::Eigen
)
