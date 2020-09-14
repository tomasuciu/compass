#ifndef FIT_HPP
#define FIT_HPP

#include <vector>
#include <optional>
#include <any>
#include <variant>
#include <numeric>
#include <algorithm>
#include <random>

#include "util.hpp"
#include "data.hpp"
#include "circle.hpp"

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

namespace compass {
using DataMatrixD = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;

template <typename Derived>
struct FitCRTP {
    Derived& derived() { return static_cast<Derived&>(*this); }
    Derived const& derived() const { return static_cast<Derived const&>(*this); }
};

template <typename Derived>
class FitBase : public FitCRTP<Derived> {
    friend class FitCRTP<Derived>;

    public:
        [[nodiscard]] inline Circle<double> getCircle() {
            std::cout << "CRTP getCircle() called" << std::endl;
            return this->derived().circle;
        }

    protected:
        Eigen::RowVector2<double> mean;
        Circle<double> circle;

        // used with FitBase::compute(data)
        FitBase() {}

        explicit FitBase(const Eigen::Ref<const DataMatrixD>& data) {
            this->derived().fit(data);
        }

        //TODO: Overload for in-place operation
        explicit FitBase(Eigen::Ref<const DataMatrixD>& data) {
            std::cout << "Calling non-const FitBase constructor" << std::endl;
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

// TODO: Consider implementing generalized algebraic fits.
template<typename Derived>
class AlgebraicFit : public FitBase<Derived> {
    friend class FitBase<Derived>;

    public:
        void fit(const Eigen::Ref<const DataMatrixD>& data) {
            this->mean = data.colwise().mean();
            this->derived().compute(data);
        }

    protected:
        AlgebraicFit() : FitBase<Derived>() {}
        AlgebraicFit(const Eigen::Ref<const DataMatrixD>& data) : FitBase<Derived>(data) {
            std::cout << "calling parameterized algebraicfit constructor" << std::endl;
        }

        void computeCircleParams(Eigen::Ref<const Eigen::MatrixXd> solution,
                Eigen::Ref<const Eigen::RowVectorXd> mean) {

            auto a = -solution(0) / 2.0 + mean(0);
            auto b = -solution(1) / 2.0 + mean(1);
            auto r = std::sqrt(std::pow(solution(0), 2) + std::pow(solution(1), 2));

            this->derived().circle.setParameters(a, b, r);
        }

        void computeCircleParams(Eigen::Ref<const Eigen::MatrixXd> solution) {
            computeCircleParams(solution, Eigen::RowVector2d::Zero());
        }
};

// TODO: Implement functionality for geometric fits where the initial guess has been precomputed
template <typename Derived>
class GeometricFitWithGuess : public FitBase<Derived>{
    friend class FitBase<Derived>;
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
         class = std::enable_if_t<std::is_base_of_v<AlgebraicFit<A>, A>>>
class GeometricFit: public FitBase<Derived> {

    friend class FitBase<Derived>;
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
    friend class FitBase<Derived>;
    typedef FitBase<Derived> Base;

    public:
        void fit(const Eigen::Ref<const DataMatrixD>& data) {
            this->mean = data.colwise().mean();
            this->derived().compute(data);
        }

    protected:
        SpecializedFit() : Base() {}
        SpecializedFit(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}
        //SpecializedFit(const Eigen::Ref<const DataMatrixD>& data, Eigen::Ref<const Eigen::Vector2d>& pole) : Base(data, pole) {}

};

template <typename Derived>
class SpecializedFitWithPole : public FitBase<Derived> {
    friend class FitBase<Derived>;
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


template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>
struct Algebraic {

using ExtendedDesignMatrix = Eigen::Matrix<T, Eigen::Dynamic, 4>; // N x 4
using DesignMatrix = Eigen::Matrix<T, Eigen::Dynamic, 3>; // N x 3

// TODO: phase out usage
using DataMatrix = Eigen::Matrix<T, 2, Eigen::Dynamic, Eigen::RowMajor>; // 2 X N
using DataMatrix3 = Eigen::Matrix<T, 3, Eigen::Dynamic>;

//TODO: Debug Nystrom functionality; not currently working
static void TaubinNystromSVD(const DataMatrix& data) {
    DataMatrix3 C(3, data.cols());
    Eigen::Vector2<T> mean{ data.row(0).mean(), data.row(1).mean() };
    auto centered = data.colwise() - mean;
    auto centeredSquared = centered.cwiseProduct(centered);
    auto Mzz = centeredSquared.colwise().sum();
    auto Zmean = Mzz.mean();
    auto Zi = Mzz.array() - Zmean;
    auto Zval = Zi/(2.0 * std::sqrt(Zmean));

    C << Zval, centered;

    // generates a random gaussian matrix; abstract away eventually
    std::default_random_engine generator;
    std::normal_distribution<double> normal{};
    auto gaussian = [&] (double) {return normal(generator);};
    int l = 10;
    Eigen::MatrixXd omega = Eigen::MatrixXd::NullaryExpr(C.rows(), l, gaussian);

    Eigen::MatrixXd CTC = C * C.transpose();
    Eigen::MatrixXd X = CTC.colPivHouseholderQr().solve(omega);

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
    qr.compute(X);

    Eigen::MatrixXd Q = qr.matrixQ();

    Eigen::MatrixXd Y = CTC.colPivHouseholderQr().solve(Q);

    Eigen::MatrixXd Z = Q.transpose() * Y;


    Eigen::LLT<Eigen::MatrixXd> cholesky = Z.llt();
    Eigen::MatrixXd G = cholesky.matrixL();
    Eigen::MatrixXd K = Y * G.inverse();

    Eigen::BDCSVD<Eigen::MatrixXd> svd(K, Eigen::ComputeFullU | Eigen::ComputeFullV);

    auto V = svd.matrixV();
    std::cout << "Nystrom SVD: \n" << V << std::endl;
    auto v = V.col(0);

    exit(1);
    auto A = V.col(2);

    auto A0 = A(0)/(2.0*std::sqrt(Zmean));
    auto A_4 = -Zmean*A0;

    auto a = -A(1)/A0/2.0 + mean(0);
    auto b = -A(2)/A0/2.0 + mean(1);
    auto radius = std::sqrt(std::pow(A(1), 2) + std::pow(A(2), 2) - 4*A0 * A_4)/std::abs(A0)/2.0;

    std::cout << a << ", " << b << " | radius : " << radius << std::endl;
    exit(1);

    std::cout << C.transpose() << std::endl;
    std::cout << CTC << std::endl;
    std::cout << omega << "\n\n";
    std::cout << X << std::endl;
    std::cout << Q << "\n\n" << std::endl;
    std::cout << Y << "\n\n" << std::endl;
    std::cout << Z << std::endl;
    std::cout << K << std::endl;
}


};

};

#endif
