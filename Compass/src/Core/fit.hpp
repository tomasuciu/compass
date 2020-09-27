#ifndef FIT_HPP
#define FIT_HPP

#include <vector>
#include <iostream>
#include <typeinfo>
#include <type_traits>

#include "Compass/src/Util/util.hpp"
//#include "Compass/src/Algebraic/Taubin.hpp"
#include "circle.hpp"

#include "FitBase.hpp"

namespace compass {

using DataMatrixD = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;

// TODO: phase out usage
using DataMatrix = Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>; // 2 X N



} // end compass

#endif
