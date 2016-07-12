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
    #include "externals/Eigen/LU"
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
    float *gsolve(int *samples, float *log_exposure, float lambda,
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
     * @brief MitsunagaNayarClassic computes the inverse CRF of a camera as a polynomial function.
     * @param samples           Sample array of size nSamples x #exposures.
     * @param nSamples          Number of samples, for each exposure.
     * @param exposures         Array of exposure timings (size: #exposures = 'Q' as in the Mitsunaga & Nayar paper).
     * @param coefficients      The output coefficients ('c' in the paper) resulting from the computation.
     * @param R                 The output estimated exposure ratios, i.e. R[q] = 'R_{q,q+1}' as in the paper.
     * @param eps               Threshold for stopping the approximation process.
     * @param max_iterations    Maximum number of iterations.
     * @return The error as in the paper.
     */
    float MitsunagaNayarClassic(int *samples, const std::size_t nSamples, const std::vector<float> &exposures,
                                std::vector<float> &coefficients, std::vector<float> &R,
                                const float eps, const std::size_t max_iterations)
    {
        float eval, val;
        const std::size_t Q = exposures.size();
        const std::size_t N = coefficients.size() - 1;

        const float Mmax = 1.f;

        for (float &_c : coefficients)
            _c = 0.f;
        R.assign(Q < 2 ? 0 : Q - 1, 1.f);
        for (int q = 0; q < R.size(); ++q)
            R[q] = exposures[q] / exposures[q+1];
        if (!samples || Q < 2 || coefficients.size() < 2)
            return std::numeric_limits<float>::infinity();

        //Precompute test with exponentials
        std::vector<Eigen::VectorXf> test(256, Eigen::VectorXf::Zero(N+1));
        for (std::size_t i = 0; i < 256; ++i) {
            test[i][0] = 1.f;
            for (std::size_t n = 1; n <= N; ++n)
                test[i][n] = (float)(i / 255.f) * test[i][n-1];
        }

        //Check valid samples
        std::vector<std::vector<float>> g(nSamples, std::vector<float>(Q, 0.f));
        std::size_t P = 0;
        for (std::size_t p = 0; p < nSamples; ++p) {
            bool valid = true;
            for (std::size_t q = 0; q < Q; ++q)
                if (samples[p * Q + q] < 0 || samples[p * Q + q] > 255) {
                    valid = false;
                    break;
                }
            if (valid) {
                for (std::size_t q = 0; q < Q; ++q)
                    g[P][q] = samples[p * Q + q] / 255.f;
                ++P;
            }
        }
        g.resize(P);

        if (g.empty())
            return std::numeric_limits<float>::infinity();

        //Precompute M with exponentials
        std::vector<std::vector<std::vector<float>>> M(P,
                                                       std::vector<std::vector<float>>(Q,
                                                                                       std::vector<float>(N+1, 0.f)));
        for (std::size_t p = 0; p < P; ++p)
            for (std::size_t q = 0; q < Q; ++q) {
                M[p][q][0] = 1.f;
                for (std::size_t n = 1; n <= N; ++n)
                    M[p][q][n] = g[p][q] * M[p][q][n-1];
            }
        g.clear();

        std::vector<std::vector<std::vector<float>>> d(P,
                                                       std::vector<std::vector<float>>(Q-1,
                                                                                       std::vector<float>(N+1, 1.f)));
        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N, N);
        Eigen::VectorXf x, b = Eigen::VectorXf::Zero(N);
        Eigen::VectorXf c(N+1), prev_c = Eigen::VectorXf::Zero(N+1);

        std::size_t iter = 0;

        do {
            //Compute d
            for (std::size_t p = 0; p < P; ++p)
                for (std::size_t q = 0; q < Q-1; ++q)
                    for (std::size_t n = 0; n <= N; ++n)
                        d[p][q][n] = M[p][q][n] - R[q] * M[p][q+1][n];

            //Build the matrix A of the linear system
            A.setZero(N, N);
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j)
                    for (std::size_t p = 0; p < P; ++p)
                        for (std::size_t q = 0; q < Q - 1; ++q)
                            A(i, j) += d[p][q][i] * (d[p][q][j] - d[p][q][N]);

            //Build the vector of knowns b
            b.setZero(N);
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t p = 0; p < P; ++p)
                    for (std::size_t q = 0; q < Q - 1; ++q)
                        b(i) -= Mmax * d[p][q][i] * d[p][q][N];

            //Solve the linear system
            x = A.partialPivLu().solve(b);
            c << x, Mmax - x.sum();

            //Evaluate approximation increment
            eval = std::numeric_limits<float>::lowest();
            for (const Eigen::VectorXf &_M : test) {
                val = std::abs((c - prev_c).dot(_M));
                if (val > eval)
                    eval = val;
            }

            prev_c = c;

            ++iter;
        } while (eval > eps && iter < max_iterations);

        for (std::size_t n = 0; n <= N; ++n)
            coefficients[n] = c[n];

        //Evaluate error
        eval = 0.f;
        for (std::size_t q = 0; q < Q-1; ++q)
            for (std::size_t p = 0; p < P; ++p) {
                val = 0.f;
                for (std::size_t n = 0; n <= N; ++n)
                    val += coefficients[n] * (M[p][q][n] - R[q] * M[p][q+1][n]);
                eval += val * val;
            }

        return eval;
    }

    /**
     * @brief MitsunagaNayarFull computes the inverse CRF of a camera as a polynomial function, using all exposure ratios.
     * @param samples           Sample array of size nSamples x #exposures.
     * @param nSamples          Number of samples, for each exposure.
     * @param exposures         Array of exposure timings (size: #exposures = 'Q' as in the Mitsunaga & Nayar paper).
     * @param coefficients      The output coefficients ('c' in the paper) resulting from the computation.
     * @param R                 The output estimated exposure ratios, i.e. R[q1][q2] = 'R_{q1,q2}' as in the book.
     * @param eps               Threshold for stopping the approximation process.
     * @param max_iterations    Maximum number of iterations.
     * @return The error as in the paper.
     */
    float MitsunagaNayarFull(int *samples, const std::size_t nSamples, const std::vector<float> &exposures,
                                std::vector<float> &coefficients, std::vector<std::vector<float>> &R,
                                const float eps, const std::size_t max_iterations)
    {
        float eval, val;
        const std::size_t Q = exposures.size();
        const std::size_t N = coefficients.size() - 1;

        const float Mmax = 1.f;

        for (float &_c : coefficients)
            _c = 0.f;
        R.assign(Q < 2 ? 0 : Q - 1, std::vector<float>(Q < 2 ? 0 : Q - 1, 1.f));
        for (int q1 = 0; q1 < R.size(); ++q1)
            for (int q2 = 0; q2 < R[q1].size(); ++q2) {
                if (q2 == q1)
                    R[q1][q2] = 1.f;
                else
                    R[q1][q2] = exposures[q1] / exposures[q2];
            }
        if (!samples || Q < 2 || coefficients.size() < 2)
            return std::numeric_limits<float>::infinity();

        //Precompute test with exponentials
        std::vector<Eigen::VectorXf> test(256, Eigen::VectorXf::Zero(N+1));
        for (std::size_t i = 0; i < 256; ++i) {
            test[i][0] = 1.f;
            for (std::size_t n = 1; n <= N; ++n)
                test[i][n] = (float)(i / 255.f) * test[i][n-1];
        }

        //Check valid samples
        std::vector<std::vector<float>> g(nSamples, std::vector<float>(Q, 0.f));
        std::size_t P = 0;
        for (std::size_t p = 0; p < nSamples; ++p) {
            bool valid = true;
            for (std::size_t q = 0; q < Q; ++q)
                if (samples[p * Q + q] < 0) {
                    valid = false;
                    break;
                }
            if (valid) {
                for (std::size_t q = 0; q < Q; ++q)
                    g[P][q] = samples[p * Q + q] / 255.f;
                ++P;
            }
        }
        g.resize(P);

        if (g.empty())
            return std::numeric_limits<float>::infinity();

        //Precompute M with exponentials
        std::vector<std::vector<std::vector<float>>> M(P,
                                                       std::vector<std::vector<float>>(Q,
                                                                                       std::vector<float>(N+1, 0.f)));
        for (std::size_t p = 0; p < P; ++p)
            for (std::size_t q = 0; q < Q; ++q) {
                M[p][q][0] = 1.f;
                for (std::size_t n = 1; n <= N; ++n)
                    M[p][q][n] = g[p][q] * M[p][q][n-1];
            }

        g.clear();

        std::vector<std::vector<std::vector<std::vector<float>>>> d(P,
                                                                    std::vector<std::vector<std::vector<float>>>(Q-1,
                                                                                std::vector<std::vector<float>>(Q-1,
                                                                                            std::vector<float>(N+1, 1.f))));
        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N, N);
        Eigen::VectorXf x, b = Eigen::VectorXf::Zero(N);
        Eigen::VectorXf c(N+1), prev_c = Eigen::VectorXf::Zero(N+1);

        std::size_t iter = 0;

        do {
            //Compute d
            for (std::size_t p = 0; p < P; ++p)
                for (std::size_t q1 = 0; q1 < Q-1; ++q1)
                    for (std::size_t q2 = 0; q2 < Q-1; ++q2) {
                        d[p][q1][q2].assign(N+1, 0.f);
                        if (q2 != q1)
                            for (std::size_t n = 0; n <= N; ++n)
                                d[p][q1][q2][n] = M[p][q1][n] - R[q1][q2] * M[p][q2][n];
                    }

            //Build the matrix A of the linear system
            A.setZero(N, N);
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j)
                    for (std::size_t p = 0; p < P; ++p)
                        for (std::size_t q1 = 0; q1 < Q - 1; ++q1)
                            for (std::size_t q2 = 0; q2 < Q - 1; ++q2)
                                if (q2 != q1)
                                    A(i, j) += d[p][q1][q2][i] * (d[p][q1][q2][j] - d[p][q1][q2][N]);

            //Build the vector of knowns b
            b.setZero(N);
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t p = 0; p < P; ++p)
                    for (std::size_t q1 = 0; q1 < Q - 1; ++q1)
                        for (std::size_t q2 = 0; q2 < Q - 1; ++q2)
                            if (q2 != q1)
                                b(i) -= Mmax * d[p][q1][q2][i] * d[p][q1][q2][N];

            //Solve the linear system
            x = A.partialPivLu().solve(b);
            c << x, Mmax - x.sum();

            //Evaluate approximation increment
            eval = std::numeric_limits<float>::lowest();
            for (const Eigen::VectorXf &_M : test) {
                val = std::abs((c - prev_c).dot(_M));
                if (val > eval)
                    eval = val;
            }

            prev_c = c;

            ++iter;
        } while (eval > eps && iter < max_iterations);

        for (std::size_t n = 0; n <= N; ++n)
            coefficients[n] = c[n];

        //Evaluate error
        eval = 0.f;
        for (std::size_t q1 = 0; q1 < Q-1; ++q1)
            for (std::size_t q2 = 0; q2 < Q-1; ++q2)
                if (q2 != q1)
                    for (std::size_t p = 0; p < P; ++p) {
                        val = 0.f;
                        for (std::size_t n = 0; n <= N; ++n)
                            val += coefficients[n] * (M[p][q1][n] - R[q1][q2] * M[p][q2][n]);
                        eval += val * val;
                    }

        return eval;
    }

    /**
     * @brief Destroy frees memory.
     */
    void Destroy()
    {
        stackOut.Destroy();
        for(unsigned int i = 0; i < icrf.size(); i++) {
            if(icrf[i] != NULL) {
                delete[] icrf[i];
            }
        }

        icrf.clear();
        poly.clear();
    }

    /**
     * @brief stackCheck
     * @param stack
     * @return
     */
    bool stackCheck(ImageVec &stack)
    {
        if(stack.size() < 2) {
            return false;
        }

        for (size_t i=1; i<stack.size(); i++) {
            if (!stack[0]->SimilarType(stack[i])) {
                return false;
            }
        }

        return true;
    }

    SubSampleStack          stackOut;
    IMG_LIN                 type_linearization;
    float                   w[256];

public:

    std::vector<float *>    icrf;

    std::vector<std::vector<float>>  poly;
    
    /**
     * @brief CameraResponseFunction
     */
    CameraResponseFunction()
    {
        type_linearization = IL_LIN;
    }

    ~CameraResponseFunction()
    {
        Destroy();
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

            case IL_POLYNOMIAL: {
                float val = 0.f;
                float M = 1.f;
                for (const float &c : poly[channel]) {
                    val += c * M;
                    M *= x;
                }
                return val;
            }
            break;

            default:
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
                constexpr float inv_gamma = 1.0f / 2.2f;
                return powf(x, inv_gamma);
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
        Destroy();

        if(!stackCheck(stack)) {
            return;
        }

        if(nSamples < 1) {
            nSamples = 256;
        }

        this->type_linearization = IL_LUT_8_BIT;

        //Subsampling the image stack
        stackOut.Compute(stack, nSamples, false);

        int *samples = stackOut.get();
        nSamples = stackOut.getNSamples();
        
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

        #ifdef PIC_DEBUG
            printf("nSamples: %d\n", nSamples);
        #endif

        int stride = nSamples * nExposure;
        for(int i = 0; i < channels; i++) {
            float *icrf_channel = gsolve(&samples[i * stride], log_exposure, lambda, nSamples,
                                        nExposure);

            icrf.push_back(icrf_channel);
        }
        
        delete[] log_exposure;
    }

    /**
     * @brief MitsunagaNayar computes the inverse CRF of a camera as a polynomial function.
     * @param stack Array of images with associated exposure. Note that this array will be sorted with increasing exposure.
     * @param polynomial_degree Degree of the polynomial. If negative, the best degree will be selected in [1, -polynomial_degree] for each channel.
     * @param nSamples Number of samples to extract from each image.
     * @param full true for computing all exposure ratios (as in book "High Dynamic Range Imaging", second edition, Reinhard et al.), false as in
     *          the original paper (only among successive exposures).
     * @param eps Threshold on the difference among successive approximations for stopping the computation.
     * @param max_iterations Stop the computation after this number of iterations.
     * @return true if successfully computed, false otherwise.
     */
    bool MitsunagaNayar(ImageVec &stack, int polynomial_degree = -3, int nSamples = 256, const bool full = false,
                        const float eps = 0.0001f, const std::size_t max_iterations = 100)
    {
        Destroy();

        if(!stackCheck(stack)) {
            return false;
        }

        if(nSamples < 1) {
            nSamples = 256;
        }

        type_linearization = IL_POLYNOMIAL;

        //Sort the array by exposure
        std::sort(stack.begin(), stack.end(), [](const Image *l, const Image *r)->bool{
            if (!l || !r) {
                return false;
            }
            return l->exposure < r->exposure;
        });

        //Subsampling the image stack
        stackOut.Compute(stack, nSamples, true);
        int *samples = stackOut.get();
        nSamples = stackOut.getNSamples();

        if (nSamples < 1) {
            return false;
        }

        //Computing CRF using Mitsunaga and Nayar
        int channels = stack[0]->channels;

        std::size_t nExposures = stack.size();

        std::vector<float> exposures(nExposures, 0.f);
        for (std::size_t t = 0; t < nExposures; ++t) {
            exposures[t] = stack[t]->exposure;
        }

        int stride = nSamples * nExposures;

        float error = std::numeric_limits<float>::infinity();
        std::vector<float> R(nExposures - 1);
        std::vector<std::vector<float>> RR(nExposures - 1, std::vector<float>(nExposures - 1));

        poly.resize(channels);

        if (polynomial_degree > 0) {
            for (int i = 0; i < channels; ++i) {
                error = 0.f;
                for (int i = 0; i < channels; ++i) {
                    poly[i].assign(polynomial_degree + 1, 0.f);
                    if (full) {
                        error += MitsunagaNayarFull(&samples[i * stride], nSamples, exposures, poly[i], RR, eps, max_iterations);
                    } else {
                        error += MitsunagaNayarClassic(&samples[i * stride], nSamples, exposures, poly[i], R, eps, max_iterations);
                    }
                }
            }
        } else if (polynomial_degree < 0) {
            error = std::numeric_limits<float>::infinity();
            std::vector<std::vector<float>> tmpCoefficients(channels);
            for (int degree = 1; degree <= -polynomial_degree && degree <= 6; ++degree) {
                float tmpError = 0.f;
                for (int i = 0; i < channels; ++i) {
                    tmpCoefficients[i].resize(degree + 1);
                    if (full) {
                        tmpError += MitsunagaNayarFull(&samples[i * stride], nSamples, exposures, tmpCoefficients[i], RR, eps, max_iterations);
                    } else {
                        tmpError += MitsunagaNayarClassic(&samples[i * stride], nSamples, exposures, tmpCoefficients[i], R, eps, max_iterations);
                    }
                }
                if (tmpError < error) {
                    error = tmpError;
                    poly = std::move(tmpCoefficients);
                    std::cout << degree << std::endl;
                }
            }
        }

        return error < std::numeric_limits<float>::infinity();
    }

    /**
     * @brief Robertson computes the CRF of a camera using all multiple exposures value Robertson et al
       1999's method (Dynamic range improvement through multiple exposures).
     * @param stack
     * @param maxIterations
     */
    void Robertson(ImageVec &stack, const size_t maxIterations = 50)
    {
        Destroy();

        if(!stackCheck(stack)) {
            return;
        }

        this->type_linearization = IL_LUT_8_BIT;

        const int channels   = stack[0]->channels;
        const int pixelcount = stack[0]->nPixels();

        // precompute robertson weighting function
        for (size_t i=0; i<256; i++) {
            this->w[i] = pic::WeightFunction(i / 255.0, pic::CW_ROBERTSON);
        }

        // avoid saturation
        int minM = 0;
        int maxM = 255;
        for (int m=0; m<256; m++) {
            if (this->w[m] > 0) {
                minM = m;
                break;
            }
        }

        for (int m=255; m>=0; m--) {
            if (this->w[m] > 0) {
                maxM = m;
                break;
            }
        }

        // avoid ghosting (for each exposure get the index for the immediately higher and lower exposure)
        int lower [stack.size()];
        int higher[stack.size()];

        for (size_t i=0; i<stack.size(); i++) {
            lower[i]  = -1;
            higher[i] = -1;
            float t = stack[i]->exposure;
            float tHigh = stack[0]->exposure;
            float tLow  = tHigh;

            for (size_t j=0; j<stack.size(); j++) {
                if (i != j) {
                    float tj = stack[j]->exposure;

                    if (tj > t && tj < tHigh) {
                        tHigh = tj;
                        higher[i] = j;
                    }
                    if (tj < t && tj > tLow) {
                        tLow = tj;
                        lower[i] = j;
                    }
                }
            }

            if (lower[i]  == -1) {
                lower[i]  = i;
            }

            if (higher[i] == -1) {
                higher[i] = i;
            }
        }

        // create initial inv response function
        {
            float * lin = new float[256];
            for (int i=0; i<256; i++) {
                lin[i] = float(2.0 * i / 255.0);
            }
            this->icrf.push_back(lin);

            for (int i=1; i<channels; i++) {
                float * col = new float[256];
                BufferAssign(col, lin, 256);
                this->icrf.push_back(col);
            }
        }

        // create quantized stack
        std::vector<unsigned char *> qstack;
        for (Image * slice : stack) {
            assert(slice->frames == 1);
            unsigned char * q = pic::ConvertHDR2LDR(slice->data, NULL, slice->size(), pic::LT_NOR);
            qstack.push_back(q);
        }

        // iterative gauss-seidel
        for (int ch=0; ch<channels; ch++) {
            float * fun = this->icrf[ch];
            float funPrev[256];
            BufferAssign(funPrev, fun, 256);

            std::vector<float> x(pixelcount);

            float prevDelta = 0.0f;
            for (size_t iter=0; iter<maxIterations; iter++) {
                // Normalize inv crf to midpoint
                {
                    // find min max
                    size_t minIdx, maxIdx;
                    for (minIdx = 0   ; minIdx < 255 && fun[minIdx]==0 ; minIdx++);
                    for (maxIdx = 255 ; maxIdx > 0   && fun[maxIdx]==0 ; maxIdx--);

                    size_t midIdx = minIdx+(maxIdx-minIdx)/2;
                    float  mid = fun[midIdx];

                    if (mid == 0.0f) {
                        // find first non-zero middle response
                        while (midIdx < maxIdx && fun[midIdx] == 0.0f) {
                            midIdx++;
                        }
                        mid = fun[midIdx];
                    }

                    if (mid != 0.0f) {
                        BufferDiv(fun, 256, mid);
                    }
                }

                // Update x
                for (int i=0; i<pixelcount; i++) {
                    float sum     = 0.0f;
                    float divisor = 0.0f;

                    float maxt = -1.0f;
                    float mint = FLT_MAX;

                    int ind = i * channels + ch;

                    for (size_t s=0; s<qstack.size(); s++) {
                        unsigned char * qslice = qstack[s];
                        const float     t      = stack[s]->exposure;

                        int m = qslice[ind];

                        // compute max/min time for under/over exposed pixels
                        if (m > maxM) {
                            mint = std::min(mint, t);
                        }

                        if (m < minM) {
                            maxt = std::max(maxt, t);
                        }

                        // to avoid ghosting
                        int mLow  = qstack[lower [s]][ind];
                        int mHigh = qstack[higher[s]][ind];
                        if (mLow > m || mHigh < m) {
                            continue;
                        }

                        const float wm = this->w[m];

                        sum     += wm * t * fun[m];
                        divisor += wm * t * t;
                    }

                    if (divisor == 0.0f) {
                        // avoid saturation
                        if (maxt > -1.0f) {
                            x[i] = fun[minM] / maxt;
                        }

                        if (mint < FLT_MAX) {
                            x[i] = fun[maxM] / mint;
                        }
                    } else if (divisor < 1e-4f) {
                        x[i] = -1.0f;
                    } else {
                        x[i] = sum / divisor;
                    }
                }

                // Update inv crf
                {
                    size_t cardEm[256] = { 0 };
                    float  sum[256]    = { 0.0f };
                    float minSatTime = FLT_MAX;
                    for (size_t s=0; s<qstack.size(); s++) {
                        unsigned char * qslice = qstack[s];
                        const float     t      = stack[s]->exposure;

                        for (int i=0; i<pixelcount; i++) {
                            if (x[i] < 0.0f) {
                                continue;
                            }

                            const int m = int(qslice[i*channels+ch]);

                            if (m == 255) {
                                if (t < minSatTime) {
                                    minSatTime = t;
                                    sum[m] = t * x[i];
                                    cardEm[m] = 1;
                                }
                                else if (t == minSatTime)
                                {
                                    sum[m] = std::min(sum[m], t * x[i]);
                                }
                            } else {
                                sum[m] += t * x[i];
                                cardEm[m]++;
                            }
                        }
                    }

                    // compute average and fill undefined values with previous one
                    float prev = 0.0f;
                    for (int m=0; m<256; m++) {
                        if (cardEm[m] != 0) {
                            fun[m] = prev = sum[m] / cardEm[m];
                        } else {
                            fun[m] = prev;
                        }
                    }
                }

                // check residuals
                {
                    static const float MaxDelta = 1e-7f;

                    float delta = 0.0f;
                    int count   = 0;
                    for (int m=0; m<256; m++) {
                        if( fun[m] != 0.0f ) {
                            float diff = fun[m] - funPrev[m];
                            delta += diff * diff;
                            funPrev[m] = fun[m];
                            count++;
                        }
                    }
                    delta /= count;

                    if (delta < MaxDelta) {
                        break;
                    }

                    prevDelta = delta;
                }
            }
        }
        // estimation complete!

        // normalize response function keeping relative scale between colors
        float maxV = -1.0f;
        for (int ch=0; ch<channels; ch++) {
            int ind;
            maxV = std::max(pic::Array<float>::getMax(this->icrf[ch], 256, ind), maxV);
        }

        for (int ch=0; ch<channels; ch++) {
            BufferDiv(this->icrf[ch], 256, maxV);
            this->icrf[ch][255] = 1.0f;
        }

        // clean quantized stack
        for (unsigned char * qslice : qstack) {
            delete qslice;
        }
    }
};

} // end namespace pic

#endif /* PIC_ALGORITHMS_CAMERA_RESPONSE_FUNCTION_HPP */

