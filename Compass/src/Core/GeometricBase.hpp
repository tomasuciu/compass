#ifndef GEOMETRIC_BASE_HPP
#define GEOMETRIC_BASE_HPP

//#include "FitBase.hpp"
//#include "AlgebraicBase.hpp"
#include <type_traits>

namespace compass {

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

template<typename Derived, class A, bool Enable>
class GeometricFit;

//template<typename Derived, class A>
//class GeometricFit<Derived, A, (std::is_base_of<AlgebraicFit<A>, A>::value)> {};

template <typename Derived, class A>
class GeometricFit<Derived, A, true> : public FitBase<Derived> {
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

template<typename Derived, class A>
class GeometricFit<Derived, A, false> {
    static_assert(False<A>{}, "Invalid parameter for algorithm of type GeometricFit");
};


//template<typename Derived, class A>
//using GeometricFitImpl = GeometricFit<Derived, A, (std::is_base_of<AlgebraicFit<A>, A>::value)>;


}

#endif /* GEOMETRIC_BASE_HPP */
