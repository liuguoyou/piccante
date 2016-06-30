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

#ifndef PIC_ALGORITHMS_CAMERA_RESPONSE_FUNCTION_HPP
#define PIC_ALGORITHMS_CAMERA_RESPONSE_FUNCTION_HPP

#include "image.hpp"
#include "point_samplers/sampler_random.hpp"
#include "histogram.hpp"
#include "filtering/filter_mean.hpp"

#include "algorithms/sub_sample_stack.hpp"
#include "algorithms/weight_function.hpp"

#ifndef PIC_DISABLE_EIGEN
    #include "externals/Eigen/SVD"
#endif

namespace pic {

enum IMG_LIN {IL_LIN, IL_2_2, IL_LUT_8_BIT, IL_POLYNOMIAL};

/**
 * @brief The CameraResponseFunction class
 */
class CameraResponseFunction
{
protected:

    /**
    * \brief gsolve computes the inverse CRF of a camera.
    */
    float *gsolve(unsigned char *samples, float *log_exposure, float lambda,
                  int nSamples, int nExposure)
    {
		#ifndef PIC_DISABLE_EIGEN

        int n = 256;
        int rows = nSamples * nExposure + n + 1;
        int cols = n + nSamples;

        #ifdef PIC_DEBUG
            printf("Matrix size: (%d, %d)\n", rows, cols);
        #endif

        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(rows, cols);
        Eigen::VectorXf b = Eigen::VectorXf::Zero(rows);

        int k = 0;

        for(int i = 0; i < nSamples; i++) {
            for(int j = 0; j < nExposure; j++) {
                int tmp = samples[i * nExposure + j];

                float w_ij = w[tmp];

                A.coeffRef(k, tmp)   =  w_ij;
                A.coeffRef(k, n + i) = -w_ij;
                
                b[k]                 =  w_ij * log_exposure[j];

                k++;
            }
        }

        A.coeffRef(k, 128) = 1.0f;
        k++;

        //Smoothness term
        for(int i = 0; i < (n - 2); i++) {
            float w_l = lambda * w[i + 1];
            A.coeffRef(k, i)     =         w_l;
            A.coeffRef(k, i + 1) = -2.0f * w_l;
            A.coeffRef(k, i + 2) =         w_l;
            k++;
        }

        //Solving the linear system
        Eigen::JacobiSVD< Eigen::MatrixXf > svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::VectorXf x = svd.solve(b);

        float *ret = new float[n];

        for(int i = 0; i < n; i++) {
            ret[i] = expf(x[i]);
        }
		#else
            float *ret = NULL;
        #endif

        return ret;
    }

    /**
     * @brief Destroy frees memory.
     */
    void Destroy()
    {
        for(unsigned int i=0; i<icrf.size(); i++) {
            if(icrf[i] != NULL) {
                delete[] icrf[i];
            }
        }
    }

    IMG_LIN                 type_linearization;
    float                   w[256];

public:

    std::vector<float *>    icrf;
    
    /**
     * @brief CameraResponseFunction
     */
    CameraResponseFunction()
    {
        type_linearization = IL_LIN;
    }

    ~CameraResponseFunction()
    {
        //Destroy();
    }

    /**
     * @brief Remove removes a camera resposnse function to a value.
     * @param x is an intensity value in [0,1].
     * @param channel
     * @return It returns x in the linear domain.
     */
    inline float Remove(float x, int channel)
    {
        switch(type_linearization) {
            case IL_LIN: {
                return x;
            }
            break;

            case IL_LUT_8_BIT: {
                int index =  CLAMP(int(std::round(x * 255.0f)), 256);
                return icrf.at(channel)[index];
            }
            break;

            case IL_2_2: {
                return powf(x, 2.2f);
            }
            break;

        }

        return x;
    }

    /**
     * @brief Apply
     * @param x a value in [0, 1]
     * @param channel
     * @return
     */
    inline float Apply(float x, int channel)
    {
        switch(type_linearization) {
            case IL_LIN: {
                return x;
            }
            break;

            case IL_LUT_8_BIT: {
                float *ptr = std::lower_bound(&icrf[channel][0], &icrf[channel][255], x);
                int offset = CLAMPi((int)(ptr - icrf[channel]), 0, 255);

                return float(offset) / 255.0f;
            }
            break;

            case IL_2_2: {
                return powf(x, 1.0f / 2.2f);
            }
            break;

        }

        return x;
    }

    /**
     * @brief FromRAWJPEG computes the CRF by exploiting the couple RAW/JPEG from cameras.
     * @param img_raw is a RAW image.
     * @param img_jpg is a JPEG compressed image.
     * @param filteringSize
     */
    void FromRAWJPEG(Image *img_raw, Image *img_jpg, int filteringSize = 11)
    {
        if((img_raw == NULL) || (img_jpg == NULL))
            return;

        if(!img_raw->SimilarType(img_jpg))
            return;
        
        icrf.clear();

        int width    = img_raw->width;
        int height   = img_raw->height;
        int channels = img_raw->channels;

        int crf_size = 256 * 256 * channels;
        unsigned int *crf = new unsigned int[crf_size];

        for(int i=0;i<crf_size;i++) {
            crf[i] = 0;
        }
               
        for(int i=0; i<height; i++) {
            for(int j=0; j<width; j++) {

                float *data_raw = (*img_raw)(j, i);
                float *data_jpg = (*img_jpg)(j, i);               

                for(int k=0;k<channels;k++) {
                    int i_raw = CLAMPi(int(255.0f * data_raw[k]), 0, 255);
                    int i_jpg = CLAMPi(int(255.0f * data_jpg[k]), 0, 255);

                    int addr = (i_raw * 256 + i_jpg ) * channels;

                    crf[addr + k ]++;
                }
            }
        }
       
        //computing the result
        std::vector< int > coords;

        for(int k=0;k<channels;k++) {

            float *ret_c = new float[256];

            for(int j=0;j<256;j++) {
                coords.clear();

                for(int i=0;i<256;i++) {

                    int addr = (i * 256 + j ) * channels + k;

                    if(crf[addr] > 0) {
                        coords.push_back(i);                        
                    }

                }

                if(!coords.empty()) {//getting the median value
                    std::sort (coords.begin(), coords.end());  
                    ret_c[j] = float(coords[coords.size() >> 1]) / 255.0f;
                }
            }
            
            if(filteringSize > 0) {
                Image toBeFiltered(1, 256, 1, 1, ret_c);

                Image *filtered = FilterMean::Execute(&toBeFiltered, NULL, filteringSize);
                
                icrf.push_back(filtered->data);

            } else {
                icrf.push_back(ret_c);
            }
        }
    }

    /**
     * @brief setCRFtoGamma2_2
     */
    void setCRFtoGamma2_2()
    {
        type_linearization = IL_2_2;
    }

    /**
     * @brief setCRFtoLinear
     */
    void setCRFtoLinear()
    {
        type_linearization = IL_LIN;
    }

    /**
     * @brief DebevecMalik computes the CRF of a camera using multiple exposures value following Debevec and Malik
    1997's method.
     * @param stack
     * @param exposure
     * @param type
     * @param nSamples
     * @param lambda
     */
    void DebevecMalik(ImageVec stack, CRF_WEIGHT type = CW_DEB97, int nSamples = 256, float lambda = 20.0f)
    {
        if(stack.empty()) {
            return;
        }

        if(nSamples < 1) {
            nSamples = 256;
        }

        icrf.clear();

        this->type_linearization = IL_LUT_8_BIT;

        //Subsampling the image stack
        unsigned char *samples = SubSampleStack::Grossberg(stack, nSamples);
        
        //Computing CRF using Debevec and Malik
        int channels = stack[0]->channels;

        //pre-computing the weight function
        for(int i = 0; i < 256; i++) {
            w[i] = WeightFunction(float(i) / 255.0f, type);
        }

        unsigned int nExposure = stack.size();

        //log domain exposure time        
        float *log_exposure = new float[nExposure];

        for(unsigned int i = 0; i < nExposure; i++) {
            log_exposure[i] = logf(stack[i]->exposure);
        }

        int stride = nSamples * nExposure;

        #ifdef PIC_DEBUG
            printf("nSamples: %d\n", nSamples);
        #endif

        for(int i = 0; i < channels; i++) {
            float *icrf_channel = gsolve(&samples[i * stride], log_exposure, lambda, nSamples,
                                        nExposure);

            /*//Wrapping into an Image for normalization
            Image img(1, 256, 1, 1, icrf_channel);

            float *max_val = img.getMaxVal(NULL, NULL);
            if(max_val[0] > 0.0f) {
                img /= max_val[0];
            }*/

            icrf.push_back(icrf_channel);
        }
        
        delete[] log_exposure;
        delete[] samples;
    }

    /**
     * @brief MitsunagaNayar
     * @param stack
     * @param polynomial_degree
     * @param nSamples
     */
    void MitsunagaNayar(ImageVec stack, int polynomial_degree = 3, int nSamples = 100)
    {
    }
};

} // end namespace pic

#endif /* PIC_ALGORITHMS_CAMERA_RESPONSE_FUNCTION_HPP */

