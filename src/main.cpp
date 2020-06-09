#include <iostream>
#include <Python.h>
#include "data.hpp"

int main() {
    auto pair = Data::RandomNormalPair<double>();
    std::cout << pair.first << std::endl;
    std::cout << pair.second << std::endl;
}
