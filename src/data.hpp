#ifndef DATA_HPP
#define DATA_HPP

#include <type_traits>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

namespace compass {

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>

struct Generate {

    inline int fastrand(unsigned int g_seed) {
        g_seed = (214013*g_seed+2531011);
        return (g_seed>>16)&0x7FFF;
    }

    // Polar form of Box-Muller transformation to generate a pair of iid random variables.
    static std::pair<T, T> RandomNormalPair() {
        T rand1, rand2, wrand;

        do {
            rand1 = (2 * static_cast<T>(rand()))/RAND_MAX - 1;
            rand2 = (2 * static_cast<T>(rand()))/RAND_MAX - 1;
            wrand = rand1*rand1 + rand2*rand2;

        } while (wrand >= 1);

        wrand = sqrt((-2*log(wrand))/wrand);
        return std::make_pair<T, T>(rand1*wrand, rand2*wrand);
    }

    /**
     * Generate data points equally spaced along a circular arc with Gaussian noise.
     * @param n         the number of points
     * @param a, b      coordinates of the center
     * @param radius    the radius of the circle
     * @param theta1    first endpoint of the arc
     * @param theta2    second endpoint of the arc
     * @param sigma     standard deviation of residuals
     */
    [[nodiscard]] static Eigen::Ref<Eigen::Matrix2X<T>> SimulateArc(const std::size_t n,
            T a, T b, T R, T theta1, T theta2, T sigma) {
        Eigen::Matrix2X<T> data = Eigen::Matrix2X<T>::Zero(2, n);

        auto dataColWise = data.colwise();

        std::transform(dataColWise.begin(), dataColWise.end(), dataColWise.begin(), [&](const auto &col) {
            static size_t counter = 0;

            T theta = theta1 + (theta2 - theta1) * counter++/(n-1);
            auto [dx, dy] = RandomNormalPair();

            Eigen::Vector2<T> updatedCol;
            auto newX = a + R*cos(theta) + sigma*dx;
            auto newY = b + R*cos(theta) + sigma*dy;
            updatedCol << newX, newY;
            return updatedCol;
        });

        return data;
    }

    // Overload of SimulateArc, taking endpoints as Cartesian coordinates
    [[nodiscard]] static Eigen::Ref<Eigen::Matrix2X<T>> SimulateArc(const std::size_t n,
            T a, T b, T R, T x1, T y1, T x2, T y2, T sigma) {
        auto theta1 = CartesianToPolar(x1, y1).second;
        auto theta2 = CartesianToPolar(x2, y2).second;
        return SimulateArc(n, a, b, R, theta1, theta2, sigma);
    }

    [[nodiscard]] static std::pair<T, T> CartesianToPolar(T x, T y) {
        auto r = std::sqrt(std::pow(x, 2) + std::pow(y, 2));
        auto theta = std::atan(y/x);
        return std::make_pair(r, theta);
    }

    /**
     * Computes the terminal endpoint of an arc on a circle, given the radius and arc angle.
     * @param a, b      coordinates of the center
     * @param x, y      coordinates of the provided endpoint
     * @param radius    the radius of the circle
     * @param theta     the arc angle, wrt the center
     */
    [[nodiscard]] static std::pair<T, T> ComputeArcEndpoint(T a, T b, T x, T y, T radius, T theta) {
        T alpha = std::atan((y - b)/(x - a));
        T x1 = a + radius * std::cos(alpha + theta);
        T y1 = b + radius * std::sin(alpha + theta);
        return std::make_pair(x1, y1);
    }


    // Simulate data points with uniform distribution in the square |x| < window, |y| < window
    [[nodiscard]] static Eigen::Ref<Eigen::Matrix2X<T>> SimulateRandom(const std::size_t n, T window) {
        Eigen::Matrix2X<T> data = Eigen::Matrix2X<T>::Zero(2, n);

        auto dataColWise = data.colwise();

        std::transform(dataColWise.begin(), dataColWise.end(), dataColWise.begin(), [&](const auto &col) {
            Eigen::Vector2<T> updatedCol;
            auto newX = window * (2* static_cast<T>(rand())/RAND_MAX - 1);
            auto newY = window * (2* static_cast<T>(rand())/RAND_MAX - 1);
            updatedCol << newX, newY;
            return updatedCol;
        });
        return data;
    }
};


}

#endif
