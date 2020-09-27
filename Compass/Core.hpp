#ifndef COMPASS_CORE_MODULE_HPP
#define COMPASS_CORE_MODULE_HPP

#include <random>
#include <type_traits>

#include <eigen-master/Eigen/SVD>
#include <eigen-master/Eigen/QR>
#include <eigen-master/Eigen/Cholesky>

// Forward declarations speed up compile time significantly and reduce the usage of circular dependencies
#include "src/Core/ForwardDeclarations.hpp"
#include "src/Core/StaticAssert.hpp"

// includes convenience functions and matrix typedefs
#include "src/Util/util.hpp"

#include "src/Core/circle.hpp"

// All algorithms derive from FitBase;
#include "src/Core/FitBase.hpp"




#endif
