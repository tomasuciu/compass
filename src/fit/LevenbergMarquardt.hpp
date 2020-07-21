#ifndef LEVENBERG_MARQUARDT_HPP
#define LEVENBERG_MARQUARDT_HPP

#include "../fit.hpp"

namespace compass {

class LevenbergMarquardtFull : public GeometricFit<LevenbergMarquardtFull> {
    public:
        LevenbergMarquardtFull() : GeometricFit<LevenbergMarquardtFull>() {
        }
};

class LevenbergMarquardtReduced : public GeometricFit<LevenbergMarquardtReduced> {};

}

#endif
