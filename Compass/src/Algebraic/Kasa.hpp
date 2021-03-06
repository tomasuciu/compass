#ifndef KASA_HPP
#define KASA_HPP

namespace compass {

class Kasa : public AlgebraicFit<Kasa> {
    friend class AlgebraicFit< Kasa >;

    public:
        Kasa() : AlgebraicFit<Kasa>() {}
        Kasa(const Eigen::Ref<const DataMatrixD>& data) : AlgebraicFit<Kasa>(data) {}

    protected:
        Kasa& compute (const Eigen::Ref<const DataMatrixD>& data) {
            Eigen::MatrixX<double> centered = data.rowwise() - mean;
            const std::size_t N = centered.rows();

            Eigen::VectorX<double> Z = (centered.array().square()).rowwise().sum();
            DesignMatrix mat = (DesignMatrix(N, 3)
                    << centered, Eigen::VectorXd::Ones(N)).finished();

            Eigen::MatrixX<double> smean = (mat.transpose() * mat).array() / N;
            clamp(smean);

            centered.conservativeResize(centered.rows(), centered.cols()+1);
            centered.col(centered.cols()-1) = Eigen::VectorXd::Zero(centered.rows());

            Eigen::MatrixXd XYZ = centered.array().colwise() * Z.array();
            Eigen::Vector3d rhs = (XYZ.colwise().sum() / -XYZ.rows());

            Eigen::LDLT<Eigen::MatrixXd> chol = smean.ldlt();
            auto sol= chol.solve(rhs);

            computeCircleParams(sol, mean);
            return *this;
        }

    private:
        //TODO: fix radius calculation for regular Kasa fit
        void computeCircleParams(Eigen::Ref<const Eigen::RowVectorXd> solVector,
                Eigen::Ref<const Eigen::RowVectorXd> mean) {
            double B = -solVector(0) / 2.0;
            double C = -solVector(1) / 2.0;
            double radius = std::sqrt(B*B + C*C + solVector(2));
            circle.setParameters(B + mean(0), C + mean(1), radius);
        }
};

class KasaConsistent : public AlgebraicFit<KasaConsistent> {
    friend class AlgebraicFit<KasaConsistent>;

    public:
        typedef AlgebraicFit<KasaConsistent> Base;

        KasaConsistent() : Base() {}
        KasaConsistent(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        KasaConsistent& compute (const Eigen::Ref<const DataMatrixD>& data) {
            size_t n = data.rows();

            // uncentered matrices
            Eigen::VectorX<double> Z = (data.array().square()).rowwise().sum();
            ExtendedDesignMatrix mat = (ExtendedDesignMatrix(n, 4)
                    << Z, data, Eigen::VectorXd::Ones(n)).finished();

            DesignMatrix uncentered = (DesignMatrix(n, 3)
                    << data, Eigen::VectorXd::Ones(n)).finished();

            Eigen::MatrixX<double> M0 = mat.transpose() * mat;

            Eigen::Matrix4<double> M1 = (Eigen::Matrix4<double>() <<
                    8 * M0(0, 3), 4 * M0(1, 3), 4* M0(2, 3), 2 * n, 4 * M0(1, 3), n, 0, 0,
                    4 * M0(2, 3), 0, n, 0, 2 * n, 0, 0, 0).finished();

            Eigen::Matrix4<double> M2 = (Eigen::Matrix4<double>() <<
                    8 * n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).finished();

            Eigen::MatrixX<double> scatter = data.transpose() * data;

            Eigen::EigenSolver<Eigen::MatrixX<double>> solver;
            solver.compute(scatter, false);
            double Vmax = solver.eigenvalues().real().minCoeff();
            double Vmin = 0.0;

            float epsilon = 0.00001 * Vmax;

            double V = 0.0;
            while (Vmax - Vmin > epsilon) {
                V = (Vmin + Vmax)/2.0;
                Eigen::MatrixX<double> M = M0 - M1*V + V*V*M2;

                auto eig = M.eigenvalues();
                if (eig.real().minCoeff() > 0) {
                    Vmin = V;
                } else {
                    Vmax = V;
                }
            }

            Eigen::MatrixX<double> matN = (Eigen::MatrixX<double>(3, 3)
                    << n, 0, 0, 0, n, 0, 0, 0, 0).finished();

            Eigen::MatrixX<double> K = (uncentered.transpose() * uncentered) - (V * matN);

            Eigen::Vector3<double> scale = (Eigen::Vector3<double>() << 0, 0, 2*V*n).finished();
            Eigen::Vector3<double> N = uncentered.transpose() * (Z.array() - 4*V).matrix() + scale;

            Eigen::LDLT<Eigen::MatrixXd> chol = K.ldlt();
            auto sol = chol.solve(N);

            auto a = sol(0)/2.0;
            auto b = sol(1)/2.0;

            std::cout << "(" << a << ", " << b << ")" << std::endl;
            auto rad = std::sqrt(std::pow(a, 2) + std::pow(b, 2) + sol(2));
            std::cout << "rad: " << rad << std::endl;

            return *this;
        }
};

}
#endif
