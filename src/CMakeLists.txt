cmake_minimum_required(VERSION 3.8)

file(GLOB fit_dir
    "*hpp"
)

add_executable(compass "main.cpp" "data.hpp" "fit.hpp" "circle.hpp" "util.hpp" ${fit_dir})

target_link_libraries(compass
    PUBLIC
        Eigen3::Eigen
)
