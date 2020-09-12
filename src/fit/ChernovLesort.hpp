#ifndef CHERNOV_LESORT_HPP
#define CHERNOV_LESORT_HPP
#include "../fit.hpp"
#include <math.h>

namespace compass {
template <class A>
class ChernovLesort : public GeometricFit<ChernovLesort<A>, A> {
    friend class GeometricFit<ChernovLesort<A>, A>;
    typedef GeometricFit<ChernovLesort<A>, A> Base;

    public:
        ChernovLesort<A>() : Base() {}
        ChernovLesort<A>(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        ChernovLesort<A>& compute(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> initalGuess) {
            const std::size_t n = data.rows() * data.cols();
            const std::size_t rows = data.rows();

            double Xshift = 0, Yshift = 0;
            double anew = initalGuess.getA() + Xshift;
            std::cout << "anew: " << anew << std::endl;
            double bnew = initalGuess.getB() + Yshift;
            std::cout << "bnew: " << bnew << std::endl;

            double Anew = 1/(2 * initalGuess.getRadius());
            std::cout << "Anew: " << Anew << std::endl;
            double aabb = anew*anew + bnew*bnew;
            std::cout << "aabb: " << aabb << std::endl;
            double Fnew = (aabb - initalGuess.getRadius()* initalGuess.getRadius()) * Anew;
            std::cout << "Fnew: " << Fnew << std::endl;
            double Tnew = std::acos(-anew/std::sqrt(aabb));
            std::cout << "Tnew: " << Tnew << std::endl;

            if (bnew > 0) {
                Tnew = 2* M_PI - Tnew;
            }

            double VarNew = computeSampleVariance(data, initalGuess);

            const int IterMAX = 50;
            double lambda = 0.01, finish = 0;
            for (int i = 0; i < IterMAX; ++i) {
                double Aold = Anew, Fold = Fnew, Told = Tnew, VarOld = VarNew;

                double H = std::sqrt(1 + 4 * Aold * Fold);
                double aold = -H * std::cos(Told)/(Aold + Aold) - Xshift;
                double bold = -H * std::sin(Told)/(Aold + Aold) - Yshift;
                double Rold = 1/ std::abs(Aold + Aold);

                // computing moments
                double DD = 1 + 4*Aold*Fold;
                double D = std::sqrt(DD);
                double CT = std::cos(Told);
                double ST = std::sin(Told);

                double H11 = 0, H12 = 0, H13 = 0, H22 = 0, H23 = 0, H33 = 0, F1 = 0, F2 = 0, F3 = 0;

                Eigen::VectorX<double> Xi = data.col(0).array() + Xshift;
                Eigen::VectorX<double> Yi = data.col(1).array() + Yshift;
                Eigen::VectorX<double> Zi = Xi.cwiseProduct(Xi).array() + Yi.cwiseProduct(Yi).array();

                Eigen::VectorX<double> Ui = (Xi.array()) * CT + (Yi.array() * ST);
                Eigen::VectorX<double> Vi = (-Xi.array()) * ST + (Yi.array() * CT);

                Eigen::VectorX<double> ADF = D*Ui.array() + (Fold + Aold*Zi(0));
                Eigen::VectorX<double> SQ = (4*Aold*ADF.array() + 1).cwiseSqrt();
                Eigen::VectorX<double> DEN = SQ.array() + 1;
                Eigen::VectorX<double> Gi = 2* ADF.array() / DEN.array();
                Eigen::VectorX<double> FACT = 2/DEN.array() * (1 - Aold*Gi.array()/SQ.array());
                //std::cout << FACT << std::endl;
                std::cout << Gi << std::endl;
                exit(1);

            }
            return *this;
        }

        double computeSampleVariance(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> guess) {
            const std::size_t rows = data.rows();

            Eigen::MatrixXd Dxy = data.rowwise() - guess.getCenter();

            Eigen::MatrixXd D = Dxy.cwiseProduct(Dxy);
            Eigen::VectorXd DSum = (D.rowwise().sum()).cwiseSqrt().array() - guess.getRadius();
            Eigen::VectorXd dtd = DSum.transpose() * DSum;
            return dtd(0) / (rows - 3);
        }
};
}
#endif /* CHERNOV_LESORT_HPP */
