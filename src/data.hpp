#include <type_traits>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

namespace circle {

/*template<typename I>
std::enable_if_t<std::is_integral_v<I>>*/

template<typename T,
typename = std::enable_if_t<std::is_floating_point_v<T>>>

class Data {
public:

    Data();
    Data(int N);

    Data(int N, std::vector<T> X, std::vector<T> Y): X(X), Y(Y) {}
    ~Data();

    void means();
    void center();
    void scale();
    void print();

private:
    int n;
    std::vector<T> X, Y;
    T meanX, meanY;
};

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

    // Simulate data points equally spaced along a circular arc with Gaussian noise.
    void SimulateArc(Data<T>& data, T a, T b, T R, T theta1, T theta2, T sigma) {
        int N = data.n;
        T theta;

        for (int i=0; i < N; ++i) {
            theta = theta1 + (theta2 - theta1)* i/(N-1);

            // isotropic Gaussian noise
            auto [dx, dy] = RandomNormalPair();

            data.X.at(i) = a + R*cos(theta) + sigma*dx;
            data.Y.at(i) = b + R*sin(theta) + sigma*dy;
        }
    }

    // Simulate data points with uniform distribution in the square |x| < window, |y| < window
    void SimulateRandom(Data<T>& data, T window) {
        for (int i = 0; i < data.n; ++i) {
            data.X.at(i) = window * (2* static_cast<T>(rand())/RAND_MAX - 1);
            data.Y.at(i) = window * (2* static_cast<T>(rand())/RAND_MAX - 1);
        }
    }

};


}

