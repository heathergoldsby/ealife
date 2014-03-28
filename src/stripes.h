//
//  stripes.h
//  ealife
//
//  Created by Heather Goldsby on 11/22/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//


#ifndef _EALIFE_STRIPES_H_
#define _EALIFE_STRIPES_H_
#include "selfrep_not_ancestor.h"
#include "resource_consumption.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/isa.h>
#include <ea/digital_evolution/spatial.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/meta_population.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/mutation.h>
#include <ea/recombination.h>

using namespace ealib;


LIBEA_MD_DECL(STRIPE_FIT, "ea.stripes.fit", int); // count the number of organisms that have the right color stripe
LIBEA_MD_DECL(ANCESTOR, "ea.stripes.ancestor", int);



/*! Compete to evolve stripes -- even number rows nand; odd number rows not
 */
template <typename EA>
struct vert_stripes : periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
    vert_stripes(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("vert_stripes.dat") {
        _df.add_field("update")
        .add_field("mean_fitness")
        .add_field("max_fitness");
    }
    
    virtual ~vert_stripes() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::max> > fit;
        
        
        // calculate "fitness":
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            int tmp_fit = 0;
            for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                
                std::string lt = get<LAST_TASK>(**j,"");
                if (((i->env().handle2ptr((*j)->location())->y) % 2) == 0) {

                    if (lt == "NAND") {
                        ++tmp_fit;
                    }
                } else {
                    if(lt == "NOT") {
                        ++tmp_fit;
                    }
                }
            }
            fit(tmp_fit);
            put<STRIPE_FIT>(tmp_fit,*i);
        }
        
        _df.write(ea.current_update())
        .write(mean(fit))
        .write(max(fit))
        .endl();
        
        std::size_t n=get<META_POPULATION_SIZE>(ea);
        typename EA::population_type offspring;
        recombine_n(ea.population(), offspring,
                    selection::tournament < access::meta_data<STRIPE_FIT> > (n, ea.population(), ea),
                    recombination::propagule_without_replacement(),
                    n, ea);
        

        //mutate
        //mutate(offspring->begin(), offspring->end(), ea);
        
        
        
        // swap populations
        std::swap(ea.population(), offspring);
        
        
    }
    
    datafile _df;
};

#endif
