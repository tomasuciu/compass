#ifndef FIT_HPP
#define FIT_HPP

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
    public:
        FitBase& fit (const Eigen::Ref<const DataMatrixD>& data) {
            std::cout << "Algebraic fit() called" << std::endl;
            this->derived().fit(data);
            return *this;
        }

        FitBase& fit (const Eigen::Ref<const DataMatrixD>& data, Circle<double> guess) {
            std::cout << "Geometric fit() called" << std::endl;
            this->derived().fit(data, guess);
            return *this;
        }

        [[nodiscard]] inline Circle<double> getCircle() {
            std::cout << "CRTP getCircle() called" << std::endl;
            return this->derived().circle;
        }

    protected:
        Eigen::RowVector2<double> mean;
        Circle<double> circle;

        FitBase() {}
        FitBase(const Eigen::Ref<const DataMatrixD>& data){
            this->mean = data.colwise().mean();

            fit(data);
        }
};


// TODO: Consider implementing generalized algebraic fits.
template<typename Derived>
class AlgebraicFit : public FitBase<Derived> {
    protected:
        AlgebraicFit() : FitBase<Derived>() {}
        AlgebraicFit(const Eigen::Ref<const DataMatrixD>& data) : FitBase<Derived>(data) {}

    protected:
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

template <typename Derived>
class GeometricFit : public FitBase<Derived> {
    protected:
        GeometricFit() {}
        GeometricFit(Eigen::Ref<DataMatrixD> data, std::any algebraicMethod) {
            //std::variant<AlgebraicFit<Kasa>> var;
            // use std::visit and CRTP to compute the initial guess; pass the class as an argument and then
            // fit. After the Circle<double> object is obtained, pass that into the constructor below as such:
            // GeometricFit(data, algebraicMethodGuess);

            //compute the intial guess using one of the algebraic fit methods, Taubin by default
        }
        GeometricFit(Eigen::Ref<DataMatrixD> data, Circle<double> guess) {
            this->derived().fit(data, guess);
        }
};

/*
    static void Trust(){}
    static void Landau(){}
    static void ChernovLesort(){}
*/

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>
struct Algebraic {

using ExtendedDesignMatrix = Eigen::Matrix<T, Eigen::Dynamic, 4>; // N x 4
using DesignMatrix = Eigen::Matrix<T, Eigen::Dynamic, 3>; // N x 3
using DataMatrix = Eigen::Matrix<T, 2, Eigen::Dynamic, Eigen::RowMajor>; // 2 X N
using DataMatrix3 = Eigen::Matrix<T, 3, Eigen::Dynamic>;

static void KukushMarkovskyHuffel(const DataMatrix& data) {}

static void PrattNewton(const Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>& data) {
    auto mean = Eigen::Vector2<T>(data.row(0).mean(), data.row(1).mean());
    auto colwiseData = data.colwise();

    std::cout << mean << std::endl;
    auto Mxx=0.0, Myy=0.0, Mxy=0.0, Mxz=0.0, Myz=0.0, Mzz=0.0;

    std::for_each(colwiseData.begin(), colwiseData.end(), [&](const auto &column) {
        auto Xi = column(0) - mean(0);
        auto Yi = column(1) - mean(1);
        auto Zi = std::pow(Xi, 2) + std::pow(Yi, 2);

        Mxx += std::pow(Xi, 2);
        Myy += std::pow(Yi, 2);
        Mzz += std::pow(Zi, 2);
        Mxy += Xi*Yi;
        Mxz += Xi*Zi;
        Myz += Yi*Zi;
    });

    int n = data.cols();

    Mxx /= n;
    Myy /= n;
    Mxy /= n;
    Mxz /= n;
    Myz /= n;
    Mzz /= n;

    // computing coefficients of characteristic polynomial

    double Mz = Mxx + Myy;
    double Cov_xy = Mxx*Myy - Mxy*Mxy;
    double Var_z = Mzz - Mz*Mz;

    double A2 = 4.0*Cov_xy - 3.0*Mz*Mz - Mzz;
    double A1 = Var_z*Mz + 4.0*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    double A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    double A22 = A2 + A2;

    double x, y, Dy, xnew, ynew;
    int iter, IterMAX=99;

    // finds root using newton's method
    for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)
        {
            Dy = A1 + x*(A22 + 16.*x*x);
            xnew = x - y/Dy;
            if ((xnew == x)||(!std::isfinite(xnew))) break;
            ynew = A0 + xnew*(A1 + xnew*(A2 + 4.0*xnew*xnew));
            if (abs(ynew)>=abs(y))  break;
            x = xnew;  y = ynew;
        }

    auto DET = x*x - x*Mz + Cov_xy;
    auto Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2.0;
    auto Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2.0;

    std::cout << "X: " << Xcenter + mean(0) << std::endl;
    std::cout << "Y: " << Ycenter + mean(1) << std::endl;
    std::cout << "Radius: " << sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz + x + x) << std::endl;
}

static Circle<T> Nievergelt(const DataMatrix& data) {
    using GolubDataMatrix = Eigen::Matrix<T, 3, Eigen::Dynamic>;
    GolubDataMatrix designMat(3, data.cols());

    Eigen::Vector2<T> mean{data.row(0).mean(), data.row(1).mean()};
    auto centered = data.colwise() - mean;

    auto centeredSquared = centered.cwiseProduct(centered);
    auto Mzz = centeredSquared.colwise().sum();

    auto Zmean = Mzz.mean();
    auto Zcentered = Mzz.array() - Zmean;

    designMat << Zcentered, centered;
    Eigen::BDCSVD<Eigen::MatrixX<T>> svd(designMat.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);

    auto V = svd.matrixV();
    auto A = V.col(2);
    auto A_4 = -Zmean * A(0);

    auto a = -A(1)/A(0)/2.0 + mean(0);
    auto b = -A(2)/A(0)/2.0 + mean(1);
    auto radius = std::sqrt(std::pow(A(1), 2) + std::pow(A(2), 2) - 4*A(0)*A_4)/std::abs(A(0))/2.0;

    return Circle<T>(a, b, radius);
}

static void TaubinNewton(const DataMatrix& data) {
    Eigen::Vector2<T> mean{data.row(0).mean(), data.row(1).mean()};
    /*auto centered = data.colwise() - mean;
    auto centeredSquared = centered.cwiseProduct(centered);
    auto Z = centeredSquared.colwise().sum();*/
    auto colwise = data.colwise();

    double Mxy = 0, Mxx = 0, Myy = 0, Mxz = 0, Myz = 0, Mzz = 0;
    std::for_each(colwise.begin(), colwise.end(), [&](const auto &column){
        auto Xi = column(0) - mean(0);
        auto Yi = column(1) - mean(1);
        auto Zi = std::pow(Xi, 2) + std::pow(Yi, 2);

        Mxy += Xi * Yi;
        Mxx += std::pow(Xi, 2);
        Myy += std::pow(Yi, 2);
        Mxz += Xi * Zi;
        Myz += Yi * Zi;
        Mzz += std::pow(Zi, 2);
    });

    int n = data.cols();

    Mxy /= n;
    Mxx /= n;
    Myy /= n;
    Mxz /= n;
    Myz /= n;
    Mzz /= n;

    // computing coefficients of characteristic polynomial
    double Mz = Mxx + Myy;
    double Cov_xy = Mxx*Myy - std::pow(Mxy, 2);

    double Var_z = Mzz - std::pow(Mz, 2);
    double A3 = 4*Mz;
    double A2 = -3*Mz*Mz - Mzz;
    double A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    double A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    double A22 = A2 + A2;
    double A33 = A3 + A3 + A3;

// finding the root of the characteristic polynomial using Newton's method starting at x = 0;
    double x, y, Dy, xnew, ynew;
    int i, iter;

    for (x=0.,y=A0,iter=0; iter < 99; iter++)  // usually, 4-6 iterations are enough
    {
    	Dy = A1 + x*(A22 + A33*x);
        xnew = x - y/Dy;
        if ((xnew == x)||(!std::isfinite(xnew))) break;
        ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
        if (abs(ynew)>=abs(y))  break;
        x = xnew;  y = ynew;
    }

// computing paramters of the fitting circle

    double DET = x*x - x*Mz + Cov_xy;
    double Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2.0;
    double Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2.0;


    double a = Xcenter + mean(0);
    double b = Ycenter + mean(1);
    double r = std::sqrt(std::pow(Xcenter, 2) + std::pow(Ycenter, 2) + Mz);
    std::cout << "(" << a << "," << b << ")" << std::endl;
    std::cout << "radius : " << r << std::endl;
}


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


static void GanderGolubStrebel(const DataMatrix& data()){}

};

};
#endif
