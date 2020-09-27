#ifndef SPECIALIZED_BASE_HPP
#define SPECIALIZED_BASE_HPP

#include "FitBase.hpp"
#include <vector>

namespace compass {
template <typename Derived>
class SpecializedFit : public FitBase<Derived> {
    typedef FitBase<Derived> Base;

    public:
        void fit(const Eigen::Ref<const DataMatrixD>& data) {
            this->mean = data.colwise().mean();
            this->derived().compute(data);
        }

    protected:
        SpecializedFit() : Base() {}
        SpecializedFit(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}
        SpecializedFit(const Eigen::Ref<const DataMatrixD>& data, Eigen::Ref<const Eigen::Vector2d>& pole) : Base(data, pole) {}

};

template <typename Derived>
class SpecializedFitWithPole : public FitBase<Derived> {
    typedef FitBase<Derived> Base;

    public:
        void fit(const Eigen::Ref<const DataMatrixD>& data, const Eigen::Ref<const Eigen::Vector2d>& pole) {
            this->mean = data.colwise().mean();
            this->derived().compute(data, pole);
        }

        void fit(const Eigen::Ref<const DataMatrixD>& data, const std::vector<double>& pole) {}

    protected:
        SpecializedFitWithPole() : Base() {}
        SpecializedFitWithPole(const Eigen::Ref<const DataMatrixD>& data, const Eigen::Ref<const Eigen::Vector2d>& pole) : Base(data, pole) {}
};

template <typename Derived>
class SpecializedFitRandomPole : public FitBase<Derived> {
    // selects one of the data points and uses that as the pole
};

}
#endif /* SPECIALIZED_BASE_HPP */
