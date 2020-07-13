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

static void clamp(Eigen::Ref<Eigen::MatrixX<double>> target, float threshold = -1e-15) {
    target = (threshold < target.array().abs()).select(target, 0.0f);
}

[[nodiscard]] static Eigen::Matrix<double, 4, 4> computeMatrixM(const DataMatrix& data){
    ExtendedDesignMatrix designMat = createExtendedDesignMatrix<double>(data);
    Eigen::Matrix<double, 4, 4> M = (designMat.transpose() * designMat).array() / data.cols();
    //instead of this, write clamp as a functor
    clamp(M);

    return M;
}

template <typename T>
static Eigen::Vector2<T> center(DataMatrix& data) {
    Eigen::Vector2<T> mean {data.row(0).mean(), data.row(1).mean()};
    data.colwise() - mean;
    return mean;
}

}

#endif
