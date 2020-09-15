//#include <Python.h>
//#include "util.hpp"
//#include "fit.hpp"

/*#include "fit/Riemann.hpp"
#include "fit/Inversion.hpp"
#include "fit/Landau.hpp"
#include "fit/Spath.hpp"
#include "fit/Kasa.hpp"
#include "fit/Taubin.hpp"
#include "fit/LevenbergMarquardt.hpp"
#include "fit/Pratt.hpp"
#include "fit/KukushMarkovskyHuffel.hpp"
#include "fit/Hyper.hpp"
#include "fit/ChernovLesort.hpp"
#include "fit/Nievergelt.hpp"
#include <Eigen/Dense>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;*/
#include <iostream>

//using namespace Eigen;
int main() {

    //auto mat = compass::Generate<double>::SimulateArc(100, 0, 0, 1, 0, M_PI/2, 2);

    //compass::DataMatrix benchmark(2, 6);

    /*benchmark << 1.0, 2.0, 5.0, 7.0, 9.0, 3.0,
              7.0, 6.0, 8.0, 7.0, 5.0, 7.0;

    std::cout << benchmark << "\n" << std::endl;

    auto A = benchmark.row(0);
    auto B = benchmark.row(1);

    std::vector<double> X(A.data(), A.data() + A.size());
    std::vector<double> Y(B.data(), B.data() + B.size());

    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> benchmarkTransposed(6, 2);
    benchmarkTransposed << benchmark.transpose();

    Eigen::Vector2d pole;
    pole << 0, 0; */

    //compass::LevenbergMarquardtReduced<compass::TaubinSVD> L;
    //L.fit(benchmarkTransposed);

    //compass::ChernovLesort<compass::TaubinSVD> C;
    //C.fit(benchmarkTransposed);

    /*compass::TaubinNewton TN;
    TN.fit(benchmarkTransposed);

    compass::Kasa K;
    K.fit(benchmarkTransposed);
    std::cout << K.getCircle() << std::endl;

    compass::KasaConsistent KC(benchmarkTransposed);
    std::cout << KC.getCircle() << std::endl;

    //compass::ChernovLesort<compass::TaubinSVD> C_T;
    //C_T.fit(benchmarkTransposed);

    compass::KukushMarkovskyHuffel KMH;
    KMH.fit(benchmarkTransposed); */

    //compass::HyperSimple h(benchmarkTransposed);
    //compass::HyperSVD H(benchmarkTransposed);
    //std::cout << compass::Kasa().fit(benchmarkTransposed).getCircle() << std::endl;
    //compass::Kasa k(benchmarkTransposed);
    //std::cout << k.getCircle() << std::endl;

    //Eigen::LDLT(benchmarkTransposed);

    //std::cout << compass::Kasa k(benchmarkTransposed) << '\n';
    //compass::Spath<compass::Kasa> s(benchmarkTransposed);
    //compass::Spath<compass::Kasa>(benchmarkTransposed).fit();
    //compass::Spath<compass::LevenbergMarquardtFull> l;

}


