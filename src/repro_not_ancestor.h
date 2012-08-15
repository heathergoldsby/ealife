//
//  repro_not_ancestor.h
//  ealife
//
//  Created by Heather Goldsby on 8/15/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_REPRO_NOT_ANCESTOR_H_
#define _EALIFE_REPRO_NOT_ANCESTOR_H_

#include <ea/artificial_life.h>
#include <ea/meta_data.h>


namespace ea {
    /*! Generates a self-replicating ancestor that performs not.
     */
    struct repro_not_ancestor {
        template <typename EA>
        typename EA::population_entry_type operator()(EA& ea) {
            typedef typename EA::representation_type representation_type;
            typename EA::individual_type ind;
            ind.name() = next<INDIVIDUAL_COUNT>(ea);
            
            representation_type& repr=ind.repr();
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), 3);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            // not
            repr[24] = ea.isa()["input"]; // input
            repr[25] = ea.isa()["input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            repr[99] =  ea.isa()["repro"]; // repro
            
            ind.hw().initialize();
            
            return make_population_entry(ind, ea);
        }
    };
    
}

#endif
