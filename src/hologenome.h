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


LIBEA_MD_DECL(GENDER, "ea.hologenome.gender", bool); // false = male; true = female
LIBEA_MD_DECL(GENDER_FIT, "ea.hologenome.gender_fit", bool); // a groups gender fitness. 

/*! Execute the next instruction if the ea is male.
 */
DIGEVO_INSTRUCTION_DECL(if_male) {
    if(get<GENDER>(ea,false) == false) {
        hw.advanceHead(Hardware::IP);
    }
}

//! Get the the gender of the ea.
DIGEVO_INSTRUCTION_DECL(get_gender){
    hw.setRegValue(hw.modifyRegister(), get<GENDER>(ea, 0));
}

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


/*! Females perform nand as many times as possible; Males perform not as many times as possible. 
 */
template <typename EA>
struct gender_compete: periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
    gender_compete(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("gender_compete.dat") {
        _df.add_field("update")
        .add_field("mean_fitness")
        .add_field("max_fitness")
        .add_field("num_females")
        .add_field("num_males")
        .add_field("mean_fitness_females")
        .add_field("max_fitness_females")
        .add_field("mean_fitness_males")
        .add_field("max_fitness_males");
    }
    
    virtual ~gender_compete() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::max> > fit;
        accumulator_set<double, stats<tag::mean, tag::max> > female_fit;
        accumulator_set<double, stats<tag::mean, tag::max> > male_fit;

        int num_females = 0;
        int num_males = 0;
        
        // compute gender fitness. for females = nand/not; for males not/nand
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {
            double score = 0.0;
            // males
            if (get<GENDER>(**i,0) == 0) {
                double nand_count = get<TASK_NAND>(**i,0);
                double not_count = get<TASK_NOT>(**i,1);
                if (nand_count && not_count) {
                    score = get<TASK_NOT>(**i,0) / get<TASK_NAND>(**i,1);
                }
                ++num_males;
                male_fit(score);
            } else { // females
                double nand_count = get<TASK_NAND>(**i,0);
                double not_count = get<TASK_NOT>(**i,1);
                if (nand_count && not_count) {
                    score = get<TASK_NAND>(**i,0) / get<TASK_NOT>(**i,1);
                }
                ++num_females;
                female_fit(score);
            }
            fit(score);
            put<GENDER_FIT>(score,**i);
        }
        
        
        _df.write(ea.current_update())
        .write(mean(fit))
        .write(max(fit))
        .write(num_females)
        .write(num_males);
        
        if (num_females) {
            _df.write(mean(female_fit))
            .write(max(female_fit));
        } else {
            _df.write(0)
            .write(0);
        }
        
        if (num_males) {
            _df.write(mean(male_fit))
            .write(max(male_fit))
            .endl();
        } else {
            _df.write(0)
            .write(0)
            .endl();
        }
        
        
        std::size_t n=get<META_POPULATION_SIZE>(ea);
        typename EA::population_type offspring;
        recombine_n(ea.population(), offspring,
                    selection::tournament < access::meta_data<GENDER_FIT> > (n, ea.population(), ea),
                    recombination::propagule_without_replacement(),
                    n, ea);
        
        // set gender
        for(typename EA::population_type::iterator i=offspring.begin(); i!=offspring.end(); ++i) {
            bool g = ea.rng().bit();
            put<GENDER>(g,**i);
        }
        
        
        // swap populations
        std::swap(ea.population(), offspring);
        
        
    }
    
    datafile _df;
};




#endif
