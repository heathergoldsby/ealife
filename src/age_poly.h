//
//  ts.h
//  ealife
//
//  Created by Heather Goldsby on 8/23/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_AGE_POLY_H_
#define _EALIFE_AGE_POLY_H_



#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/algorithm/string.hpp>

#include "resource_consumption.h"
#include "configurable_mutation.h"


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/isa.h>
#include <ea/digital_evolution/discrete_spatial_environment.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>

using namespace ealib;
using namespace boost::accumulators;


LIBEA_MD_DECL(LAST_TASK, "ea.ape.last_task", std::string);
LIBEA_MD_DECL(GERM_MUTATION_PER_SITE_P, "ea.ape.germ_mutation_per_site_p", double);
LIBEA_MD_DECL(TASK_LETHALITY_PROB, "ea.ape.task_lethality_prob", double);
LIBEA_MD_DECL(EACH_TASK_THRESH, "ea.ape.each_task_thresh", double);

LIBEA_MD_DECL(NOT_LETHALITY_PROB, "ea.ape.not_lethality_prob", double);
LIBEA_MD_DECL(NAND_LETHALITY_PROB, "ea.ape.nand_lethality_prob", double);
LIBEA_MD_DECL(AND_LETHALITY_PROB, "ea.ape.and_lethality_prob", double);
LIBEA_MD_DECL(ORNOT_LETHALITY_PROB, "ea.ape.ornot_lethality_prob", double);
LIBEA_MD_DECL(OR_LETHALITY_PROB, "ea.ape.or_lethality_prob", double);
LIBEA_MD_DECL(ANDNOT_LETHALITY_PROB, "ea.ape.andnot_lethality_prob", double);
LIBEA_MD_DECL(NOR_LETHALITY_PROB, "ea.ape.nor_lethality_prob", double);
LIBEA_MD_DECL(XOR_LETHALITY_PROB, "ea.ape.xor_lethality_prob", double);
LIBEA_MD_DECL(EQUALS_LETHALITY_PROB, "ea.ape.equals_lethality_prob", double);

LIBEA_MD_DECL(NOT_AGE, "ea.ape.not_age", double);
LIBEA_MD_DECL(NAND_AGE, "ea.ape.nand_age", double);
LIBEA_MD_DECL(AND_AGE, "ea.ape.and_age", double);
LIBEA_MD_DECL(ORNOT_AGE, "ea.ape.ornot_age", double);
LIBEA_MD_DECL(OR_AGE, "ea.ape.or_age", double);
LIBEA_MD_DECL(ANDNOT_AGE, "ea.ape.andnot_age", double);
LIBEA_MD_DECL(NOR_AGE, "ea.ape.nor_age", double);
LIBEA_MD_DECL(XOR_AGE, "ea.ape.xor_age", double);
LIBEA_MD_DECL(EQUALS_AGE, "ea.ape.equals_age", double);

LIBEA_MD_DECL(RES_INITIAL_AMOUNT, "ea.ape.res_initial_amount", double);
LIBEA_MD_DECL(RES_INFLOW_AMOUNT, "ea.ape.res_inflow_amount", double);
LIBEA_MD_DECL(RES_OUTFLOW_FRACTION, "ea.ape.res_outflow_fraction", double);
LIBEA_MD_DECL(RES_FRACTION_CONSUMED, "ea.ape.res_fraction_consumed", double);
LIBEA_MD_DECL(ANCESTOR, "ea.ape.ancestor", int);







/*! This event triggers a task to have a lethal consequence.
 */

template <typename EA>
struct task_lethality : reaction_event<EA> {
    
    task_lethality(EA& ea) : reaction_event<EA>(ea) {
    }
    
    virtual ~task_lethality() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer
                            double r, // amount of resource consumed
                            EA& ea) {
        put<LAST_TASK>(task->name(), ind);
        
        // Grab this task's lethality load
        double lethality_prob = get<TASK_LETHALITY_PROB>(*task);
        // Check if this org dies as a consequence of performing the task...
        if (ea.rng().p(lethality_prob)) {
            ind.alive() = false;
        }
        
    }
};

/*! Tracks the first age at which an organism performed a task.
 */

template <typename EA>
struct task_first_age : reaction_event<EA> {
    task_first_age(EA& ea) : reaction_event<EA>(ea) {
    }
    
    virtual ~task_first_age() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer
                            double r,
                            EA& ea) {
        std::string t = task->name();
        if (t == "not") {
            if (!exists<NOT_AGE>(ind)) {
                put<NOT_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "nand") {
            if (!exists<NAND_AGE>(ind)) {
                put<NAND_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "and") {
            if (!exists<AND_AGE>(ind)) {
                put<AND_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "ornot") {
            if (!exists<ORNOT_AGE>(ind)) {
                put<ORNOT_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "or") {
            if (!exists<OR_AGE>(ind)) {
                put<OR_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "andnot") {
            if (!exists<ANDNOT_AGE>(ind)) {
                put<ANDNOT_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "nor") {
            if (!exists<NOR_AGE>(ind)) {
                put<NOR_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "xor") {
            if (!exists<XOR_AGE>(ind)) {
                put<XOR_AGE>(ind.hw().age(), ind) ;
            }
        } else if (t == "equals") {
            if (!exists<EQUALS_AGE>(ind)) {
                put<EQUALS_AGE>(ind.hw().age(), ind) ;
            }
        }
    }
};


/*! Prints information about the mean age at which tasks are first performed at the metapop level.
 */


template <typename EA>
struct task_first_age_tracking : end_of_update_event<EA> {
    task_first_age_tracking(EA& ea) : end_of_update_event<EA>(ea), _df("tasks_first_age.dat") {
        _df.add_field("update")
        .add_field("sub_pop_size")
        .add_field("pop_size")
        .add_field("not_age")
        .add_field("nand_age")
        .add_field("and_age")
        .add_field("ornot_age")
        .add_field("or_age")
        .add_field("andnot_age")
        .add_field("nor_age")
        .add_field("xor_age")
        .add_field("equals_age")
        ;
    }
    
    //! Destructor.
    virtual ~task_first_age_tracking() {
    }
    
    //! Track resources!
    virtual void operator()(EA& ea) {
        if ((ea.current_update() % 100) == 0) {
            double t_not = 0.0;
            double t_not_count = 0.0;
            double t_nand = 0.0;
            double t_nand_count = 0.0;
            double t_and = 0.0;
            double t_and_count = 0.0;
            double t_ornot = 0.0;
            double t_ornot_count = 0.0;
            double t_or = 0.0;
            double t_or_count = 0.0;
            double t_andnot = 0.0;
            double t_andnot_count = 0.0;
            double t_nor = 0.0;
            double t_nor_count = 0.0;
            double t_xor = 0.0;
            double t_xor_count = 0.0;
            double t_equals = 0.0;
            double t_equals_count = 0.0;
            
            int sub_pop_size = 0;
            int pop_size = 0;
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                ++sub_pop_size;
                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j){
                    typename EA::individual_type::individual_type& ind=**j;
                    
                    if (ind.alive()) {
                    
                        ++pop_size;
                    
                        if (exists<NOT_AGE>(ind)) {
                            t_not += get<NOT_AGE>(ind);
                            ++t_not_count;
                        }
                        if (exists<NAND_AGE>(ind)) {
                            t_nand += get<NAND_AGE>(ind);
                            ++t_nand_count;
                        }
                        if (exists<AND_AGE>(ind)) {
                            t_and += get<AND_AGE>(ind);
                            ++t_and_count;
                        }
                        if (exists<ORNOT_AGE>(ind)) {
                            t_ornot += get<ORNOT_AGE>(ind);
                            ++t_ornot_count;
                        }
                        if (exists<OR_AGE>(ind)) {
                            t_or += get<OR_AGE>(ind);
                            ++t_or_count;
                        }
                        if (exists<ANDNOT_AGE>(ind)) {
                            t_andnot += get<ANDNOT_AGE>(ind);
                            ++t_andnot_count;
                        }
                        if (exists<NOR_AGE>(ind)) {
                            t_nor += get<NOR_AGE>(ind);
                            ++t_nor_count;
                        }
                        if (exists<XOR_AGE>(ind)) {
                            t_xor += get<XOR_AGE>(ind);
                            ++t_xor_count;
                        }
                        if (exists<EQUALS_AGE>(ind)) {
                            t_equals += get<EQUALS_AGE>(ind);
                            ++t_equals_count;
                        }
                    }
                    
                }
            }
            
            if (t_not_count) { t_not /= t_not_count; }
            if (t_nand_count) { t_nand /= t_nand_count; }
            if (t_and_count) { t_and /= t_and_count; }
            if (t_ornot_count) { t_ornot /= t_ornot_count; }
            if (t_or_count) { t_or /= t_or_count; }
            if (t_andnot_count) { t_andnot /= t_andnot_count; }
            if (t_nor_count) { t_nor /= t_nor_count; }
            if (t_xor_count) { t_xor /= t_xor_count; }
            if (t_equals_count) { t_equals /= t_equals_count; }

            
            _df.write(ea.current_update())
            .write(sub_pop_size)
            .write(pop_size)
            .write(t_not)
            .write(t_nand)
            .write(t_and)
            .write(t_ornot)
            .write(t_or)
            .write(t_andnot)
            .write(t_nor)
            .write(t_xor)
            .write(t_equals)
            .endl();
        }
        
    }
    datafile _df;
    
};


//! Performs group replication using germ lines.
template <typename EA>
struct ape_two_task_replication : end_of_update_event<EA> {
    //! Constructor.
    ape_two_task_replication(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ape_two_task_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (i->population().size() < 2) continue;
            
            double t_not = get<TASK_NOT>(*i, 0.0);
            double t_nand = get<TASK_NAND>(*i, 0.0);
            
            // change this snippet...
            // each task must have been done more than EACH_TASK_THRESH
            if ((t_not > get<EACH_TASK_THRESH>(*i)) &&
                (t_nand > get<EACH_TASK_THRESH>(*i))){
                
                
                typename EA::individual_type::individual_type prop = (*i).founder();
                prop.repr().resize((*i).founder().hw().original_size());
                
                //typename EA::individual_type::individual_type j = **(i->population().begin());
                //typename EA::individual_type::individual_type prop = j;
                //prop.repr().resize(j.hw().original_size());
                prop.hw().initialize();
                
                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p = ea.make_individual();
                
                // mutate it:
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
                mutate(prop,m,*p);
                
                // and fill up the offspring population with copies of the germ:
                typename EA::individual_type::individual_ptr_type o=p->make_individual(prop.repr());
                p->append(o);
                offspring.push_back(p);
                
                // reset resource units
                i->env().reset_resources();
                put<GROUP_RESOURCE_UNITS>(0,*i);
                put<TASK_NOT>(0.0,*i);
                put<TASK_NAND>(0.0,*i);
                
                
                // i == parent individual;
                typename EA::population_type parent_pop, offspring_pop;
                parent_pop.push_back(*i.base());
                offspring_pop.push_back(p);
                inherits(parent_pop, offspring_pop, ea);
            }
        }
        
        
        // select surviving parent groups
        if (offspring.size() > 0) {
            int n = get<META_POPULATION_SIZE>(ea) - offspring.size();
            
            typename EA::population_type survivors;
            select_n<selection::random>(ea.population(), survivors, n, ea);
            
            // add the offspring to the list of survivors:
            survivors.insert(survivors.end(), offspring.begin(), offspring.end());
            
            // and swap 'em in for the current population:
            std::swap(ea.population(), survivors);
        }
        
    }
    
};



//! Performs group replication using germ lines and three tasks.
template <typename EA>
struct ape_three_task_replication : end_of_update_event<EA> {
    //! Constructor.
    ape_three_task_replication(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ape_three_task_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (i->population().size() < 2) continue;
            
            double t_not = get<TASK_NOT>(*i, 0.0);
            double t_nand = get<TASK_NAND>(*i, 0.0);
            double t_ornot = get<TASK_ORNOT>(*i, 0.0);
            
            // change this snippet...
            // each task must have been done more than EACH_TASK_THRESH
            if ((t_not > get<EACH_TASK_THRESH>(*i)) &&
                (t_nand > get<EACH_TASK_THRESH>(*i)) &&
                (t_ornot> get<EACH_TASK_THRESH>(*i))){
                
                
                typename EA::individual_type::individual_type prop = (*i).founder();
                prop.repr().resize((*i).founder().hw().original_size());
                
                //typename EA::individual_type::individual_type j = **(i->population().begin());
                //typename EA::individual_type::individual_type prop = j;
                //prop.repr().resize(j.hw().original_size());
                prop.hw().initialize();
                
                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p = ea.make_individual();
                
                // mutate it:
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
                mutate(prop,m,*p);
                
                // and fill up the offspring population with copies of the germ:
                typename EA::individual_type::individual_ptr_type o=p->make_individual(prop.repr());
                p->append(o);
                offspring.push_back(p);
                
                // reset resource units
                i->env().reset_resources();
                put<GROUP_RESOURCE_UNITS>(0,*i);
                put<TASK_NOT>(0.0,*i);
                put<TASK_NAND>(0.0,*i);
                put<TASK_ORNOT>(0.0,*i);
                
                
                // i == parent individual;
                typename EA::population_type parent_pop, offspring_pop;
                parent_pop.push_back(*i.base());
                offspring_pop.push_back(p);
                inherits(parent_pop, offspring_pop, ea);
            }
        }
        
        
        // select surviving parent groups
        if (offspring.size() > 0) {
            int n = get<META_POPULATION_SIZE>(ea) - offspring.size();
            
            typename EA::population_type survivors;
            select_n<selection::random>(ea.population(), survivors, n, ea);
            
            // add the offspring to the list of survivors:
            survivors.insert(survivors.end(), offspring.begin(), offspring.end());
            
            // and swap 'em in for the current population:
            std::swap(ea.population(), survivors);
        }
        
    }
    
};

//! Performs group replication using germ lines.
template <typename EA>
struct ape_lr_replication : end_of_update_event<EA> {
    //! Constructor.
    ape_lr_replication(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ape_lr_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (i->population().size() < 2) continue;
            
            if (exists<GROUP_RESOURCE_UNITS>(*i) &&
                (get<GROUP_RESOURCE_UNITS>(*i) > get<GROUP_REP_THRESHOLD>(*i))){
                
                // grab a copy of the founder!
                
                typename EA::individual_type::individual_type prop = (*i).founder();
                prop.repr().resize((*i).founder().hw().original_size());
                
                prop.hw().initialize();
                
                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p = ea.make_individual();
                
                // mutate it:
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
                mutate(prop,m,*p);
                
                // and fill up the offspring population with copies of the germ:
                typename EA::individual_type::individual_ptr_type o=p->make_individual(prop.repr());
                p->append(o);
                offspring.push_back(p);
                
                // reset resource units
                i->env().reset_resources();
                put<GROUP_RESOURCE_UNITS>(0,*i);
                
                // i == parent individual;
                typename EA::population_type parent_pop, offspring_pop;
                parent_pop.push_back(*i.base());
                offspring_pop.push_back(p);
                inherits(parent_pop, offspring_pop, ea);
            }
        }
        
        
        // select surviving parent groups
        if (offspring.size() > 0) {
            int n = get<META_POPULATION_SIZE>(ea) - offspring.size();
            
            typename EA::population_type survivors;
            select_n<selection::random>(ea.population(), survivors, n, ea);
            
            // add the offspring to the list of survivors:
            survivors.insert(survivors.end(), offspring.begin(), offspring.end());
            
            // and swap 'em in for the current population:
            std::swap(ea.population(), survivors);
        }
        
    }
    
    
    
};



#endif
