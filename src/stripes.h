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
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/discrete_spatial_environment.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/mutation.h>
#include <ea/recombination.h>

using namespace ealib;


LIBEA_MD_DECL(STRIPE_FIT, "ea.stripes.fit", int); // count the number of organisms that have the right color stripe
LIBEA_MD_DECL(ANCESTOR, "ea.stripes.ancestor", int);

LIBEA_MD_DECL(NUM_PROPAGULE_CELL, "ea.stripes.num_propagule_cell", int);




/*! Compete to evolve stripes -- even number rows nand; odd number rows not
 */
template <typename EA>
struct permute_stripes : periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
    permute_stripes(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("permute_stripes.dat") {
        _df.add_field("update")
        .add_field("mean_fitness")
        .add_field("max_fitness")
        .add_field("mean_one_fitness")
        .add_field("max_one_fitness")
        .add_field("mean_two_fitness")
        .add_field("max_two_fitness")
        .add_field("mean_three_fitness")
        .add_field("max_three_fitness")
        .add_field("mean_four_fitness")
        .add_field("max_four_fitness")
        .add_field("mean_five_fitness")
        .add_field("max_five_fitness")
        .add_field("mean_six_fitness")
        .add_field("max_six_fitness")
        .add_field("mean_num_not")
        .add_field("max_num_not")
        .add_field("mean_num_nand")
        .add_field("max_num_nand")
        .add_field("mean_num_org")
        .add_field("max_num_org");

    }
    
    virtual ~permute_stripes() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::max> > fit;
        accumulator_set<double, stats<tag::mean, tag::max> > one_fit;
        accumulator_set<double, stats<tag::mean, tag::max> > two_fit;
        accumulator_set<double, stats<tag::mean, tag::max> > three_fit;
        accumulator_set<double, stats<tag::mean, tag::max> > four_fit;
        accumulator_set<double, stats<tag::mean, tag::max> > five_fit;
        accumulator_set<double, stats<tag::mean, tag::max> > six_fit;
        accumulator_set<double, stats<tag::mean, tag::max> > num_not;
        accumulator_set<double, stats<tag::mean, tag::max> > num_nand;
        accumulator_set<double, stats<tag::mean, tag::max> > num_org;

        
        
        
        // calculate "fitness":
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            // vert stripes
            double one_fit_not = 0;
            double one_fit_nand = 0;
            double two_fit_not = 0;
            double two_fit_nand = 0;
            // horizontal stripes
            double three_fit_not = 0;
            double three_fit_nand = 0;
            double four_fit_not = 0;
            double four_fit_nand = 0;
            // diagonal stripes
            double five_fit_not = 0;
            double five_fit_nand = 0;
            double six_fit_not = 0;
            double six_fit_nand = 0;
            
            double tmp_num_not = 0;
            double tmp_num_nand = 0;
            
            int tmp_num_org = 0;


            for(typename EA::individual_type::ea_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                ++tmp_num_org;
                std::string lt = get<LAST_TASK>(**j,"");
                if (((i->ea().env().location((**j).position())->y) % 2) == 0) {
                    if (lt == "nand") {
                        ++one_fit_nand;
                        ++tmp_num_nand;
                    }
                    if (lt == "not") {
                        ++two_fit_not;
                        ++tmp_num_not;
                    }
                } else {
                    if(lt == "not") {
                        ++one_fit_not;
                        ++tmp_num_not;
                    }
                    if (lt == "nand") {
                        ++two_fit_nand;
                        ++tmp_num_nand;
                    }
                }
                
                if (((i->ea().env().location((**j).position())->x) % 2) == 0) {
                    if (lt == "nand") {
                        ++three_fit_nand;
                    }
                    if (lt == "not") {
                        ++four_fit_not;
                    }
                } else {
                    if(lt == "not") {
                        ++three_fit_not;
                    }
                    if (lt == "nand") {
                        ++four_fit_nand;
                    }
                }
                
                
                if (((i->ea().env().location((**j).position())->x) % 2) ==
                    ((i->ea().env().location((**j).position())->y) % 2)) {
                    
                    if(lt == "not") {
                        ++five_fit_not;
                    }
                    if (lt == "nand") {
                        ++six_fit_nand;
                    }
                } else {
                    if(lt == "nand") {
                        ++five_fit_nand;
                    }
                    if (lt == "not") {
                        ++six_fit_not;
                    }
                }

            }
            double tmp_one_fit = (one_fit_not + 1)  * (one_fit_nand + 1);
            double tmp_two_fit = (two_fit_not + 1)  * (two_fit_nand + 1);
            double tmp_three_fit = (three_fit_not + 1)  * (three_fit_nand + 1);
            double tmp_four_fit = (four_fit_not + 1)  * (four_fit_nand + 1);
            double tmp_five_fit = (five_fit_not + 1)  * (five_fit_nand + 1);
            double tmp_six_fit = (six_fit_not + 1)  * (six_fit_nand + 1);
            double tmp_fit = std::max(tmp_one_fit, tmp_two_fit);
            tmp_fit = std::max(tmp_fit, tmp_three_fit);
            tmp_fit = std::max(tmp_fit, tmp_four_fit);
            tmp_fit = std::max(tmp_fit, tmp_five_fit);
            tmp_fit = std::max(tmp_fit, tmp_six_fit);
            
            
            fit(tmp_fit);
            one_fit(tmp_one_fit);
            two_fit(tmp_two_fit);
            three_fit(tmp_three_fit);
            four_fit(tmp_four_fit);
            five_fit(tmp_five_fit);
            six_fit(tmp_six_fit);
            
            num_org(tmp_num_org);
            num_nand(tmp_num_nand);
            num_not(tmp_num_not);


            put<STRIPE_FIT>(tmp_fit,*i);
        }
        
        

        
        _df.write(ea.current_update())
        .write(mean(fit))
        .write(max(fit))
        .write(mean(one_fit))
        .write(max(one_fit))
        .write(mean(two_fit))
        .write(max(two_fit))
        .write(mean(three_fit))
        .write(max(three_fit))
        .write(mean(four_fit))
        .write(max(four_fit))
        .write(mean(five_fit))
        .write(max(five_fit))
        .write(mean(six_fit))
        .write(max(six_fit))
        .write(mean(num_not))
        .write(max(num_not))
        .write(mean(num_nand))
        .write(max(num_nand))
        .write(mean(num_org))
        .write(max(num_org))
        .endl();
        
        std::size_t n=get<META_POPULATION_SIZE>(ea);
        typename EA::population_type offspring; // container of (pointers to) subpopulations
        recombine_n(ea.population(), offspring,
                    selection::tournament < access::meta_data<STRIPE_FIT> > (n, ea.population(), ea),
                    recombination::propagule_without_replacement(),
                    n, ea);
        

        
        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));

        
        // Mutate and fill each offspring group.
        for(typename EA::population_type::iterator i=offspring.begin(); i!=offspring.end(); ++i) {
            assert((*i)->ea().population().size() == 1);
            
            // mutate it:
            mutate(**((*i)->ea().population().begin()),m,(*i)->ea());
            typename EA::individual_type::ea_type::individual_type g = (**((*i)->ea().population().begin()));
            
            // and fill up the offspring population with copies of the germ:
            for (int k=1; k<get<NUM_PROPAGULE_CELL>(ea); ++k) {
                typename EA::individual_type::ea_type::individual_ptr_type o = (*i)->ea().copy_individual(g);
                (*i)->insert((*i)->end(), o);
            }
        }
        

        

        //mutate
        //mutate(offspring.begin(), offspring.end(), ea);
        

        
        
        // swap populations
        std::swap(ea.population(), offspring);
        
        
    }
    
    datafile _df;
};

#endif
