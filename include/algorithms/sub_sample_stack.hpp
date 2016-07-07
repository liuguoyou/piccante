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

#ifndef PIC_ALGORITHMS_SUB_SAMPLE_STACK_HPP
#define PIC_ALGORITHMS_SUB_SAMPLE_STACK_HPP

#include "util/math.hpp"

#include "image.hpp"
#include "point_samplers/sampler_random.hpp"
#include "histogram.hpp"

namespace pic {

/**
 * @brief The SubSampleStack class
 */
class SubSampleStack
{
public:
    
    /**
     * @brief SubSampleStack
     */
    SubSampleStack(ImageVec stack, int &nSamples)
    {
        
    }

    /**
    * \brief This function creates a low resolution version of the stack using Grossberg and Nayar sampling.
    * \param stack is a stack of Image* at different exposures
    * \param nSamples output number of samples
    * \return samples an array of unsigned char values which is the low resolution stack
    */
    static unsigned char *Grossberg(ImageVec stack, int &nSamples)
    {
        if(stack.size() < 1) {
            return NULL;
        }

        if(nSamples < 1) {
            nSamples = 256;
        }

        int channels  = stack[0]->channels;
        unsigned int exposures = stack.size();

        Histogram *h = new Histogram[exposures * channels];

        int c = 0;

        #ifdef PIC_DEBUG
            printf("Computing histograms...");
        #endif

        for(int j = 0; j < channels; j++) {
            for(unsigned int i = 0; i < exposures; i++) {
                h[c].Calculate(stack[i], VS_LDR, 256, j);
                h[c].cumulativef(true);
                c++;
            }
        }

        #ifdef PIC_DEBUG
            printf("Ok\n");
        #endif

        unsigned char *samples = new unsigned char[nSamples * channels * exposures];

        #ifdef PIC_DEBUG
            printf("Sampling...");
        #endif

        c = 0;

        float div = float(nSamples - 1);

        for(int k = 0; k < channels; k++) {
            for(int i = 0; i < nSamples; i++) {

                float u = float(i) / div;

                for(unsigned int j = 0; j < exposures; j++) {

                    int ind = k * exposures + j;

                    float *bin_c = h[ind].getCumulativef();

                    float *ptr = std::upper_bound(&bin_c[0], &bin_c[255], u);

                    samples[c] = CLAMPi((int)(ptr - bin_c), 0, 255);
                    c++;
                }
            }
        }

        #ifdef PIC_DEBUG
            printf("Ok\n");
        #endif

        delete[] h;

        return samples;
    }

    /**
     * @brief Spatial creates a low resolution version of the stack.
    * \param stack is a stack of Image* at different exposures
    * \param nSamples output number of samples
    * \return samples an array of unsigned char values which is the low resolution stack
     * @return
     */
    static unsigned char *Spatial(ImageVec stack, int &nSamples, SAMPLER_TYPE sub_type = ST_MONTECARLO_S)
    {
        if(stack.size() < 1) {
            return NULL;
        }

        if(nSamples < 1) {
            nSamples = 256;
        }

        int width    = stack[0]->width;
        int height   = stack[0]->height;
        int channels = stack[0]->channels;

        Vec<2, int> vec(width, height);

        RandomSampler<2> *sampler = new RandomSampler<2>(sub_type, vec, nSamples, 1, 0);

        #ifdef PIC_DEBUG
            int oldNSamples = nSamples;
        #endif

        nSamples = sampler->getSamplesPerLevel(0);

        #ifdef PIC_DEBUG
            printf("--subSample samples: %d \t \t old samples: %d\n", nSamples, oldNSamples);
        #endif

        unsigned char *samples = new unsigned char[nSamples * channels * stack.size()];

        int c = 0;

        for(int k = 0; k < channels; k++) {
            for(int i = 0; i < nSamples; i++) {

                int x, y;
                sampler->getSampleAt(0, i, x, y);

                for(unsigned int j = 0; j < stack.size(); j++) {
                    float fetched = (*stack[j])(x, y)[k];
                    float tmp = lround(fetched * 255.0f);
                    samples[c] = CLAMPi(int(tmp), 0, 255);
                    c++;
                }
            }
        }

        delete sampler;

        return samples;
    }
};

} // end namespace pic

#endif /* PIC_ALGORITHMS_SUB_SAMPLE_STACK_HPP */

