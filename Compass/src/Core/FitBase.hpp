#ifndef FIT_BASE_HPP
#define FIT_BASE_HPP

//#include "circle.hpp"

namespace compass {

//using DataMatrixD = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;

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

}

#endif /* FIT_BASE_HPP */
