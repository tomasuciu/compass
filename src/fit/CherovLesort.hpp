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
            const std::size_t n = data.rows();

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

            std::cout << Tnew << std::endl;

            if (bnew > 0) {
                Tnew = 2* M_PI - Tnew;
            }

            double VarNew = computeSampleVariance(data, initalGuess);

            return *this;
        }

        double computeSampleVariance(const Eigen::Ref<const DataMatrixD>& data, const Circle<double> guess) {
            const std::size_t n = data.rows();

            Eigen::MatrixXd Dxy = data.rowwise() - guess.getCenter();

            std::cout << Dxy << std::endl;
            //auto D = (Dxy.cwiseProduct(Dxy)).rowwise().sum().cwiseSqrt().array() - guess.getRadius();
            //std::cout << D << std::endl;

            //auto variance = (D.transpose() * D).array() / (n - 3);
            //std::cout << D.transpose() * D << std::endl;
            return 0;
        }
};
}
#endif /* CHERNOV_LESORT_HPP */
