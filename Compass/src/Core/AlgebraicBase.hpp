#ifndef ALGEBRAIC_BASE_HPP
#define ALGEBRAIC_BASE_HPP

//#include "FitBase.hpp"

namespace compass {

template<typename Derived>
class AlgebraicFit : public FitBase<Derived> {
    public:
        void fit(const Eigen::Ref<const DataMatrixD>& data) {
            this->mean = data.colwise().mean();
            this->derived().compute(data);
        }

    protected:
        AlgebraicFit() : FitBase<Derived>() {}
        AlgebraicFit(const Eigen::Ref<const DataMatrixD>& data) : FitBase<Derived>(data) {}

        void computeCircleParameters(Eigen::Ref<const Eigen::VectorXd> v, bool recenter=true) {
            auto a = -v(1) / v(0) / 2.0;
            auto b = -v(2) / v(0) / 2.0;
            auto radius = std::sqrt(v(1) * v(1) + v(2) * v(2) - 4*v(0) * v(3)) / std::abs(v(0)) / 2.0;

            if (recenter) {
                this->derived().circle.setParameters(a + this->mean(0), b + this->mean(1), radius);
            } else {
                this->derived().circle.setParameters(a, b, radius);
            }
        }
};
}

#endif /* ALGEBRAIC_BASE_HPP */
