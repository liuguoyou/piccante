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

#ifndef PIC_ALGORITHMS_WEIGHT_FUNCTION_HPP
#define PIC_ALGORITHMS_WEIGHT_FUNCTION_HPP

namespace pic {

/**
 * @brief The CRF_WEIGHT enum
 */
enum CRF_WEIGHT {CW_ALL, CW_HAT, CW_DEB97, CW_DEB97p01, CW_ROBERTSON, CW_GAUSS};

/**
 * @brief WeightFunction computes weight functions for x in [0,1].
 * @param x is an input value in [0, 1].
 * @param type is the type of the function.
 * @return It returns a weight for x.
 */
inline float WeightFunction(float x, CRF_WEIGHT type)
{
    switch(type) {

    case CW_GAUSS: {
        float sigma = 0.5f;
        float mu = 0.5f;
        float sigma_sq_2 = 2.0f * (sigma * sigma);
        float x_mu = (x - mu);
        return expf(-4.0f * (x_mu * x_mu) / sigma_sq_2);
    }
    break;

    case CW_ROBERTSON: {
        float sigma = 0.5f;
        float mu = 0.5f;
        float mu_sq = mu * mu;
        float sigma_sq_2 = 2.0f * (sigma * sigma);

        float x_mu = (x - mu);

        float y =  expf(-4.0f * (x_mu * x_mu) / sigma_sq_2);

        float shift_val = expf(-4.0f * mu_sq / sigma_sq_2);
        //float scale_val = expf(0.0f); --> 1.0

        y = (y - shift_val) / (1.0f - shift_val);
        return  CLAMPi(y, 0.0f, 1.0f);
    }
    break;

    case CW_HAT: {
        float val = (2.0f * x - 1.0f);
        float val_squared = val * val;
        float val_quartic = val_squared * val_squared;
        return (1.0f - val_quartic * val_quartic * val_quartic);
    }

    case CW_DEB97: {
        float Zmin = 0.0f;
        float Zmax = 1.0f;
        float tr = (Zmin + Zmax) / 2.0f;

        if(x <= tr) {
            return x - Zmin;
        } else {
            return Zmax - x;
        }
    }
    break;

    case CW_DEB97p01: {
        float Zmin = 0.01f;
        float Zmax = 0.99f;
        float tr = (Zmin + Zmax) / 2.0f;

        if(x <= tr) {
            return x - Zmin;
        } else {
            return Zmax - x;
        }
    }

    }

    return 1.0f;
}

} // end namespace pic

#endif /* PIC_ALGORITHMS_WEIGHT_FUNCTION_HPP */

