/*

PICCANTE
The hottest HDR imaging library!
http://vcg.isti.cnr.it/piccante

Copyright (C) 2014
Visual Computing Laboratory - ISTI CNR
http://vcg.isti.cnr.it
First author: Francesco Banterle

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

*/

#ifndef PIC_UTIL_POLYNOMIAL_HPP
#define PIC_UTIL_POLYNOMIAL_HPP


#ifndef PIC_DISABLE_EIGEN

#include <vector>

#include "externals/Eigen/QR"

namespace pic {

/**
 * @brief polynomialFit
 * @param x
 * @param y
 * @param n is the degree of the polynomial
 * @return
 */
std::vector<float> polynomialFit(std::vector<float> &x, std::vector<float> &y, unsigned int n)
{
    std::vector<float> poly;

    if(n == 0) {
        return poly;
    }

    if(x.size() != y.size()) {
        return poly;
    }

    unsigned int np1 = n + 1;

    unsigned int s = x.size();
    Eigen::MatrixXf A(s, np1);
    Eigen::VectorXf b(s);

    for(unsigned int i = 0; i < s; i++) {
        b(i) = y[i];
        A(i, n) = 1.0f;
    }

    for(int j = (n - 1); j >= 0; j--) {
        for(unsigned int i = 0; i < s; i++) {
            A(i, j) = x[i] * A(i, j + 1);
        }
    }

    Eigen::VectorXf _x = A.colPivHouseholderQr().solve(b);

    for(int i = n; i >= 0; i--) {
        poly.push_back(_x(i));
    }

    return poly;
}

/**
 * @brief polyval
 * @param poly
 * @return
 */
float polynomialVal(std::vector< float > & poly, float x)
{
    float val = 0.f;
    float M = 1.f;
    for (const float &c : poly) {
        val += c * M;
        M *= x;
    }
    return val;
}

} // end namespace pic

#endif

#endif //PIC_UTIL_POLYNOMIAL_HPP
