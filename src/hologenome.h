//
//  hologenome.h
//  ealife
//
//  Created by Heather Goldsby on 7/8/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_HOLOGENOME_H_
#define _EALIFE_HOLOGENOME_H_
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




/*! Competition to perform NOT as many times as possible.
 */
template <typename EA>
struct not_strength : periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
    not_strength(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("not_strength.dat") {
        _df.add_field("update")
        .add_field("mean_fitness")
        .add_field("max_fitness");
    }
    
    virtual ~not_strength() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::max> > fit;
        
        // calculate "fitness":
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {
            fit(get<TASK_NOT>(**i,0));
        }
        
        _df.write(ea.current_update())
        .write(mean(fit))
        .write(max(fit))
        .endl();
             
        std::size_t n=get<META_POPULATION_SIZE>(ea);
        typename EA::population_type offspring;
        recombine_n(ea.population(), offspring,
                    selection::tournament < access::meta_data<TASK_NOT> > (n, ea.population(), ea),
                    recombination::propagule_without_replacement(),
                    n, ea);
        
        // mutate? if desired?
        
        // swap populations
        std::swap(ea.population(), offspring);
  

    }
    
    datafile _df;
};

/*! Competition to maximize the product of nots and nands.
 */

LIBEA_MD_DECL(NOT_NAND_PROD, "ea.hologenome.not_nand_prod", double);


template <typename EA>
struct not_nand_balance : periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
    not_nand_balance(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("not_nand_balance.dat") {
        _df.add_field("update")
        .add_field("mean_fitness")
        .add_field("max_fitness");
    }
    
    virtual ~not_nand_balance() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::max> > fit;
        
        // calculate "fitness":
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {
            double nn = get<TASK_NOT>(**i,0) * get<TASK_NAND>(**i,0);
            put<NOT_NAND_PROD>(nn, **i);
            fit(nn);
        }
        
        _df.write(ea.current_update())
        .write(mean(fit))
        .write(max(fit))
        .endl();
        
        std::size_t n=get<META_POPULATION_SIZE>(ea);
        typename EA::population_type offspring;
        recombine_n(ea.population(), offspring,
                    selection::tournament < access::meta_data<NOT_NAND_PROD> > (n, ea.population(), ea),
                    recombination::propagule_without_replacement(),
                    n, ea);
        
        // mutate? if desired?
        
        // swap populations
        std::swap(ea.population(), offspring);
        
        
    }
    
    datafile _df;
};


#endif
