#ifndef GANDER_GOLUB_STREBEL_HPP
#define GANDER_GOLUB_STREBEL_HPP
#include "../fit.hpp"

namespace compass {
class GanderGolubStrebel : public AlgebraicFit<GanderGolubStrebel> {
    friend class AlgebraicFit<GanderGolubStrebel>;
    typedef AlgebraicFit<GanderGolubStrebel> Base;

    public:
        GanderGolubStrebel() : Base() {}
        GanderGolubStrebel(const Eigen::Ref<const DataMatrixD>& data) : Base(data) {}

    protected:
        GanderGolubStrebel& compute (const Eigen::Ref<const DataMatrixD>& data) {
            // TODO implement!
            return *this;
        }
};
}
#endif /*GANDER_GOLUB_STREBEL_HPP */
