#ifndef PRATT_HPP
#define PRATT_HPP

#include <Eigen/SVD>
#include "../fit.hpp"

namespace compass {

class PrattSVD : public AlgebraicFit<PrattSVD> {
    public:
        PrattSVD() : AlgebraicFit<PrattSVD>() {}
        PrattSVD(Eigen::Ref<DataMatrixD> data) : AlgebraicFit<PrattSVD>(data) {}

        PrattSVD& fit(Eigen::Ref<DataMatrixD> data) {
            this->mean = center<double>(data);
            Eigen::MatrixX<double> centered = data;

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            DesignMatrix mat = (DesignMatrix(centered.rows(), 3)
                    << centered, Eigen::VectorXd::Ones(centered.rows())).finished();

            Eigen::BDCSVD<Eigen::MatrixX<double>> svd{mat, Eigen::ComputeThinU | Eigen::ComputeThinV};

            return *this;
        }
};

class PrattNewton : public AlgebraicFit<PrattNewton> {};

// Possibly using regula falsi with Illinois modification (i.e., introduces quasi-regularization parameters to
// scale function values according to iteration history)
class PrattRobust : public AlgebraicFit<PrattRobust> {};
}

#endif
