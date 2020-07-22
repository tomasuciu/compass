#pragma once
#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <Eigen/Dense>

namespace compass {
//TODO: Figure out generic type harmony between util and fit
using DataMatrix = Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::RowMajor>;
using DataMatrix3 = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;

using DesignMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using ExtendedDesignMatrix = Eigen::Matrix<double, Eigen::Dynamic, 4>;

/*template <typename T>
[[nodiscard]] static DesignMatrix createDesignMatrix(Eigen::Ref<Eigen::MatrixX<double>> target) {
    const int rows = target.rows();
    return Eigen::Matrix<T, rows, target.cols() + 1> (Eigen::Matrix<T, rows, target.cols() + 1>()
            << target, Eigen::MatrixX<T>::Ones(target.rows())).finished();

}*/

template <typename T>
[[nodiscard]] static ExtendedDesignMatrix createExtendedDesignMatrix(const DataMatrix& data) {
    // TODO: Utilize the center(DataMatrix&) function instead
    Eigen::Vector2<T> mean {data.row(0).mean(), data.row(1).mean()};
    Eigen::Matrix2X<T> centered = data.colwise() - mean;

    Eigen::Matrix2X<T> Z = centered.cwiseProduct(centered);
    Eigen::VectorX<T> Mz = Z.colwise().sum();

    ExtendedDesignMatrix dMat(data.cols(), 4);
    dMat << Mz, centered.transpose(), Eigen::Vector<T, 6>::Ones();
    return dMat;
}

static void clamp(Eigen::Ref<Eigen::MatrixX<double>> target, float threshold = 1e-15) {
    target = (threshold < target.array().abs()).select(target, 0.0f);
}

template <typename T>
[[nodiscard]] static Eigen::RowVector2<T> center(Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>> data) {
    Eigen::RowVector2<T> mean {data.col(0).mean(), data.col(1).mean()};
    data = data.rowwise() - mean;
    return mean;
}

[[nodiscard]] static Eigen::Matrix<double, 4, 4> computeMatrixM(const DataMatrix& data){
    ExtendedDesignMatrix designMat = createExtendedDesignMatrix<double>(data);
    Eigen::Matrix<double, 4, 4> M = (designMat.transpose() * designMat).array() / data.cols();
    //instead of this, write clamp as a functor
    clamp(M);

    return M;
}

//TODO: implement support for different norms - hence p parameter
template <typename T>
[[nodiscard]] static T computeNorm (Eigen::VectorX<double> vec, int p=2) {
    return vec.cwiseAbs2().sum();
}


}

#endif
