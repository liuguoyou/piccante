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

#ifndef PIC_FEATURES_MATCHING_BINARY_FEATURE_MATCHER
#define PIC_FEATURES_MATCHING_BINARY_FEATURE_MATCHER

#include <vector>

namespace pic{

/**
 * @brief The BinaryFeatureMatcher class
 */
class BinaryFeatureMatcher
{
protected:
    std::vector<unsigned int *> *descs;
    unsigned int desc_size;

public:

    /**
     * @brief BinaryFeatureMatcher
     * @param descs
     * @param n
     */
    BinaryFeatureMatcher(std::vector<unsigned int *> *descs, unsigned int desc_size)
    {
        this->desc_size = desc_size;
        this->descs = descs;
    }

    /**
     * @brief getMatch
     * @param desc0
     * @param matched_j
     * @param dist_1
     * @return
     */
    virtual bool getMatch(unsigned int *desc, int &matched_j, unsigned int &dist_1)
    {
        return false;
    }

    void getAllMatches(std::vector<unsigned int *> *descs0, std::vector< Eigen::Vector3i > &matches)
    {
        matches.clear();

        for(unsigned int i = 0; i< descs0->size(); i++) {
            int matched_j;
            unsigned int dist_1;

            if(getMatch(descs0->at(i), matched_j, dist_1)) {
                matches.push_back(Eigen::Vector3i(i, matched_j, dist_1));
            }
        }
    }
};

}

#endif // PIC_FEATURES_MATCHING_BINARY_FEATURE_MATCHER
