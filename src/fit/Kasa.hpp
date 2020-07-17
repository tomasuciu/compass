#ifndef KASA_HPP
#define KASA_HPP

//TODO: Temporary, eventually specify /I flag with CMake
#include "../fit.hpp"
namespace compass {

class Kasa : public AlgebraicFit<Kasa> {
    public:
        Kasa() : AlgebraicFit<Kasa>() {}
        Kasa(const DataMatrixD& data) : AlgebraicFit<Kasa>(data) {}

        Kasa& fit (const DataMatrixD& data) {
            Eigen::RowVector2<double> mean {data.col(0).mean(), data.col(1).mean()};
            Eigen::MatrixX<double> centered = data.rowwise() - mean;

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            DesignMatrix mat = (DesignMatrix(centered.rows(), 3)
                    << centered, Eigen::VectorXd::Ones(centered.rows())).finished();

            Eigen::MatrixX<double> smean = (mat.transpose() * mat).array() / data.rows();
            clamp(smean);

            centered.conservativeResize(centered.rows(), centered.cols()+1);
            centered.col(centered.cols()-1) = Eigen::VectorXd::Zero(centered.rows());

            Eigen::MatrixXd XYZ = centered.array().colwise() * Z.array();
            Eigen::Vector3d rhs = (XYZ.colwise().sum() / -XYZ.rows());

            Eigen::LDLT<Eigen::MatrixXd> chol = smean.ldlt();
            auto sol= chol.solve(rhs.transpose());

            auto B = -1*sol(0)/2.0;
            auto C = -1*sol(1)/2.0;

            std::cout << B + mean(0) << ", " << C + mean(1) << std::endl;

            return *this;
        }
    private:
};

class KasaConsistent : public AlgebraicFit<KasaConsistent> {};


}
#endif
