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

#ifndef PIC_BRUTE_FORCE_FEATURE_MATCHER_BINARY_HPP
#define PIC_BRUTE_FORCE_FEATURE_MATCHER_BINARY_HPP

#include <vector>

namespace pic{

/**
 * @brief The BruteForceFeatureMatcherBinary class
 */
class BruteForceFeatureMatcherBinary
{
protected:

    std::vector<unsigned int *> *descs;
    unsigned int n;

public:

    /**
     * @brief BruteForceFeatureMatcherBinary
     * @param descs
     * @param n
     */
    BruteForceFeatureMatcherBinary(std::vector<unsigned int *> *descs, unsigned int n)
    {
        this->n = n;
        this->descs = descs;
    }

    /**
     * @brief getMatch
     * @param desc0
     * @param matched_j
     * @param dist_1
     * @return
     */
    bool getMatch(unsigned int *desc, int &matched_j, unsigned int &dist_1)
    {
        unsigned int dist_2 = 0;

        dist_1 = 0;

        matched_j = -1;

        for(unsigned int j = 0; j < descs->size(); j++) {
            unsigned int dist = BRIEFDescriptor::match(desc, descs->at(j), n);

            if(dist > 0) {
                if(dist > dist_1) {
                    dist_2 = dist_1;
                    dist_1 = dist;
                    matched_j = j;
                } else {
                    if(dist > dist_2) {
                        dist_2 = dist;
                    }
                }
            }
        }
        return ((dist_1 * 100 > dist_2 * 105) && matched_j != -1);
    }
};

}

#endif // PIC_BRUTE_FORCE_FEATURE_MATCHER_BINARY_HPP
