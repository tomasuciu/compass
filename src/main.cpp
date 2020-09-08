#include <Python.h>
#include "util.hpp"
#include "fit.hpp"

#include "fit/Spath.hpp"
#include "fit/Kasa.hpp"
#include "fit/Taubin.hpp"
#include "fit/LevenbergMarquardt.hpp"
#include "fit/Pratt.hpp"
#include "fit/Hyper.hpp"

#include <Eigen/Dense>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace Eigen;
int main() {

    //auto mat = compass::Generate<double>::SimulateArc(100, 0, 0, 1, 0, M_PI/2, 2);

    compass::DataMatrix benchmark(2, 6);

    benchmark << 1.0, 2.0, 5.0, 7.0, 9.0, 3.0,
              7.0, 6.0, 8.0, 7.0, 5.0, 7.0;

    std::cout << benchmark << "\n" << std::endl;

    auto A = benchmark.row(0);
    auto B = benchmark.row(1);

    std::vector<double> X(A.data(), A.data() + A.size());
    std::vector<double> Y(B.data(), B.data() + B.size());

    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> benchmarkTransposed(6, 2);
    benchmarkTransposed << benchmark.transpose();


    //compass::HyperSimple h(benchmarkTransposed);
    //compass::HyperSVD H(benchmarkTransposed);

    compass::Spath<compass::Kasa> s;
    //compass::Spath<compass::LevenbergMarquardtFull> l;

}


