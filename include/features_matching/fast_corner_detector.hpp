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

#ifndef PIC_FEATURES_MATCHING_FAST_CORNER_DETECTOR_HPP
#define PIC_FEATURES_MATCHING_FAST_CORNER_DETECTOR_HPP

#include "util/vec.hpp"

#include "image.hpp"
#include "filtering/filter_luminance.hpp"

#include "features_matching/general_corner_detector.hpp"

#ifndef PIC_DISABLE_EIGEN
#include "externals/Eigen/Dense"
#endif

namespace pic {

#ifndef PIC_DISABLE_EIGEN

class FastCornerDetector: public GeneralCornerDetector
{
protected:
    Image *lum_flt;
    bool      bComputeThreshold;

    float     sigma, threshold;
    int       radius;

public:
    FastCornerDetector(float sigma = 1.0f, int radius = 1, float threshold = 0.001f) : GeneralCornerDetector()
    {
        bComputeThreshold = true;
        Update(sigma, radius, threshold);
    }

    ~FastCornerDetector()
    {
        if(bLum) {
            delete lum;
        }
    }

    void Update(float sigma = 1.0f, int radius = 1, float threshold = 0.001f)
    {
        if(sigma > 0.0f) {
            this->sigma = sigma;
        } else {
            this->sigma = 1.0f;
        }

        if(radius > 0) {
            this->radius = radius;
        } else {
            this->radius = 1;
        }

        if(threshold > 0.0f) {
            this->threshold = threshold;
        } else {
            this->threshold = 0.001f;
        }
    }

    void Compute(Image *img, std::vector< Eigen::Vector3f > *corners)
    {
        if(img == NULL) {
            return;
        }

        if(img->channels == 1) {
            bLum = false;
            lum = img;
        } else {
            bLum = true;
            lum = FilterLuminance::Execute(img, lum, LT_CIE_LUMINANCE);
        }

        if(corners == NULL) {
            corners = new std::vector< Eigen::Vector3f >;
        }

        corners->clear();

        //Filtering the image
        FilterGaussian2D flt(sigma);
        lum_flt = flt.ProcessP(Single(lum), NULL);

        int x[] = {0, 1, 2, 3, 3,  3,  2,  1,  0, -1, -2, -3, -3, -3, -2, -1};
        int y[] = {3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1,  0,  1,  2,  3};

        int width  = lum_flt->width;
        int height = lum_flt->height;

        bool *corners_map = new bool[width * height];

        memset(corners_map, 0, sizeof(bool) * width *  height);

        Image V(1,width, height, 1);
        V.SetZero();

        for(int i=3; i<(height - 3); i++) {
            for(int j=3; j<(width - 3); j++) {
                int ind = i * width + j;

                float p = lum_flt->data[ind];

                bool bDark[16];
                bool bBright[16];
                float v[16];

                float sum = p;

                for(int k=0; k<16; k++){
                    float value = (*lum_flt)(j + x[k], i + y[k])[0];
                    v[k] = value;
                    sum += value;
                }

                //computing the threshold
                float thr;
                if(bComputeThreshold) {
                    thr = 0.2f * sum / 16.0f;
                    thr = (thr > 1e-9f ) ? thr : threshold;

                } else {
                    thr = threshold;
                }

                //testing
                float p_thr_dark   = p - thr;
                float p_thr_bright = p + thr;

                for(int k=0; k<16; k++){
                    bDark[k]   = v[k] <= p_thr_dark;
                    bBright[k] = v[k] >= p_thr_bright;
                }

                //first test: 0, 4, 8, 12
                int cDark   = bDark[0]   + bDark[4]   + bDark[8]   + bDark[12];
                int cBright = bBright[0] + bBright[4] + bBright[8] + bBright[12];

                if((cDark < 3) && (cBright < 3) ){
                    //corners_map[ind] = false;
                    continue;
                }

                //second test: check for 12 continuous true values in bDark or cBright
                int counter_dark   = 0;
                int counter_bright = 0;

                //corners_map[ind] = false;
                for(int k=0; k<16; k++){
                    if(bDark[k]){
                        counter_dark++;
                    } else {
                        counter_dark = 0;
                    }

                    if(bBright[k]){
                        counter_bright++;
                    } else {
                        counter_dark = 0;
                    }

                    if((counter_bright > 11) || (counter_dark > 11)) {
                        corners_map[ind] = true;

                        //Computing V function
                        float V_dark   = 0.0f;
                        float V_bright = 0.0f;
                        for(int k=0; k<16; k++){
                            if(bDark[k]) {
                                V_bright += fabsf(v[k] - p) - thr;
                            }

                            if(bBright[k]) {
                                V_dark += fabsf(p - v[k]) - thr;
                            }
                        }
                        V.data[ind] = MAX(V_bright, V_dark);

                        break;
                    }
                }
            }
        }

        //Maximal supression
        int side = radius * 2 + 1;
        int *indices = new int[side * side];

        for(int i=3; i<(height - 3); i++) {
            for(int j=3; j<(width - 3); j++) {
                int ind = i * width + j;

                if(!corners_map[ind]){
                    continue;
                }

                indices[0] = ind;
                int counter = 1;

                for(int k = -radius; k <= radius; k++) {
                    int yy = CLAMP(i + k, height);

                    for(int l = -radius; l <= radius; l++) {

                        if((l == 0)&&(k == 0)) {
                            continue;
                        }

                        int xx = CLAMP(j + l, width);

                        ind = yy * width + xx;

                        if(corners_map[ind]){
                            indices[counter] = ind;
                            counter++;
                        }

                    }
                }

                //are other corners near-by?
                if(counter > 1) {
                    //finding the maximum value
                    float V_value = V.data[indices[0]];
                    int index = 0;

                    for(int k=1; k<counter; k++){
                        if(V.data[indices[k]] > V_value) {
                            V_value = V.data[indices[k]];
                            index = k;
                        }
                    }

                    if(index == 0){
                        corners->push_back(Eigen::Vector3f (float(j), float(i), 1.0f) );
                    }
                } else {
                    corners->push_back(Eigen::Vector3f (float(j), float(i), 1.0f) );
                }

            }
        }
    }
\
};

#endif

} // end namespace pic

#endif /* PIC_FEATURES_MATCHING_FAST_CORNER_DETECTOR_HPP */

