#ifndef FIT_HPP
#define FIT_HPP

#include <vector>
#include <iostream>
#include <typeinfo>
#include <type_traits>

#include "Compass/src/Util/util.hpp"
#include "circle.hpp"

namespace compass {

using DataMatrixD = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;

// TODO: phase out usage
using DataMatrix = Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>; // 2 X N

template <typename Derived>
struct FitCRTP {
    protected:
        Derived& derived() { return static_cast<Derived&>(*this); }
        Derived const& derived() const { return static_cast<Derived const&>(*this); }

        FitCRTP() = default;
        ~FitCRTP() = default;
};

template <typename Derived>
class FitBase : public FitCRTP<Derived> {
    public:
        [[nodiscard]] Circle<double> getCircle() {
            return this->derived().circle;
        }

        [[nodiscard]] Eigen::Vector3<double> getVector() {
            return this->derived().circle.getVector();
        }

    protected:
        Eigen::RowVector2<double> mean;
        Circle<double> circle;

        FitBase() {}

        explicit FitBase(const Eigen::Ref<const DataMatrixD>& data) {
            this->derived().fit(data);
        }

        //TODO: Overload for in-place operation
        explicit FitBase(Eigen::Ref<const DataMatrixD>& data) {
            this->derived().fit(data);
        }

        // used for GeometricFitWithGuess; initial guess is precomputed
        explicit FitBase(const Eigen::Ref<const DataMatrixD>& data, Circle<double> initialGuess) {
            this->derived().fit(data, initialGuess);
        }

        // SpecializedFitWithPole, const
        explicit FitBase(const Eigen::Ref<const DataMatrixD>& data, const Eigen::Ref<const Eigen::Vector2<double>>& pole) {
            this->derived().fit(data, pole);
        }
};


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

// TODO: Implement functionality for geometric fits where the initial guess has been precomputed
template <typename Derived>
class GeometricFitWithGuess : public FitBase<Derived>{
    typedef FitBase<Derived> Base;

    public:
        void fit(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialGuess) {
            this->mean = data.colwise().mean();
            this->derived().compute(data, initialGuess);
        }
    protected:
        GeometricFitWithGuess() : Base() {}
        GeometricFitWithGuess(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialGuess) : Base(data, initialGuess) {}

};

template <typename Derived>
class GeometricFitGuessIndiscriminate : public FitBase<Derived> {
    // specifically created for algorithms that are gauranteed to converge with any initial guess.
    // Notably : Chernov-Lesort and Chernov-Houssam
};

template <typename Derived, class A,
         class = std::enable_if_t<std::is_base_of<AlgebraicFit<A>, A>::value>>
class GeometricFit: public FitBase<Derived> {
    typedef FitBase<Derived> Base;

    public:
        void fit(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialGuess) {
            this->mean = data.colwise().mean();
            this->derived().compute(data, initialGuess);
        }

        void fit(const Eigen::Ref<const DataMatrixD>& data) {
            Circle<double> initalGuess = A(data).getCircle();
            fit(data, initalGuess);
        }

    protected:
        GeometricFit() : Base(){}
        GeometricFit(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

        // used in cases where initial guess is precomputed
        GeometricFit(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initialguess) : Base(data, initialguess) {}

};

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

} // end compass

#endif
