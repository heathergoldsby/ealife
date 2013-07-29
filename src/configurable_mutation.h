//
//  germ_mutation.h
//  ealife
//
//  Created by Heather Goldsby on 10/4/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_CONFIGURABLE_MUTATION_H_
#define _EALIFE_CONFIGURABLE_MUTATION_H_


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/isa.h>
#include <ea/digital_evolution/spatial.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/meta_population.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>

using namespace ealib;


/*! Mutation - Per-site mutation at a configurable rate
 */
struct configurable_per_site {            
    typedef mutation::uniform_integer mutation_type;
    
    configurable_per_site(double prob) : _mp(prob) {
    }
    
    //! Iterate through all elements in the given representation, possibly mutating them.
    template <typename EA>
    void operator()(typename EA::individual_type& ind, EA& ea) {
        typename EA::representation_type& repr=ind.repr();
        for(typename EA::representation_type::iterator i=repr.begin(); i!=repr.end(); ++i){
            if(ea.rng().p(_mp)) {
                _mt(repr, i, ea);
            }
        }
    }
    
    mutation_type _mt;
    double _mp; //! Mutation probability
};


#endif
