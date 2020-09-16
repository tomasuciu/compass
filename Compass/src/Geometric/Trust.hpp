#ifndef TRUST_HPP
#define TRUST_HPP

#include "Compass/src/Core/fit.hpp"

namespace compass {
template <class A>
class Trust : public GeometricFit<Trust<A>, A> {
    friend class GeometricFit<Trust<A>, A>;
    typedef GeometricFit<Trust<A>, A> Base;

    public:
        Trust<A>() : Base() {}
        Trust<A>(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        Trust& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double>& initialGuess) {

            return *this;
        }
    private:
        void computeStep() {}

};

}

#endif /* TRUST_HPP */
