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

#ifndef PIC_FEATURES_MATCHING_BINARY_FEATURE_LSH_MATCHER_HPP
#define PIC_FEATURES_MATCHING_BINARY_FEATURE_LSH_MATCHER_HPP

#include <vector>

#include "features_matching/hash_table_lsh.hpp"

namespace pic {

/**
 * @brief The LSH class
 */
class BinaryFeatureLSHMatcher
{
protected:
    std::vector< HashTableLSH* > tables;

public:

    /**
     * @brief LSH
     */
    BinaryFeatureLSHMatcher(std::vector< unsigned int *> *descs, unsigned int desc_size, unsigned int nTables = 8, unsigned int hash_size = 8)
    {
        std::mt19937 m_rnd(1);

        printf("C");

        for(unsigned int i=0; i < nTables; i++) {
            printf("B");

            unsigned int n = desc_size * sizeof(unsigned int) * 8;
            unsigned int *g_f = getHash(m_rnd, n, hash_size);
            HashTableLSH *tmp = new HashTableLSH(hash_size, g_f, descs, desc_size);
            tables.push_back(tmp);
            printf("A");
        }
    }

    /**
     * @brief getHash
     * @param dim
     * @param hash_size
     * @param seed
     * @return
     */
    static unsigned int *getHash(std::mt19937 &m, unsigned int dim, unsigned int hash_size = 0)
    {
        if(hash_size == 0) {
            hash_size = 8;
        }

        unsigned int *out = new unsigned int[hash_size];

        std::set<unsigned int> tmp;

        int c = 0;
        while(tmp.size() < dim) {
            unsigned int val = m() % dim;
            auto result = tmp.insert(val);

            if(result.second) {
                out[c] = val;
                c++;
            }
        }

        return out;
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
        dist_1 = 0;
        matched_j = -1;

        for(unsigned int i=0; i<tables.size(); i++) {
            tables[i]->getNearest(desc, matched_j, dist_1);
        }

        return true;
    }
};

} // end namespace pic

#endif /* PIC_FEATURES_MATCHING_BINARY_FEATURE_LSH_MATCHER_HPP */

