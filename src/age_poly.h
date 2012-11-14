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

#include "multi_birth_selfrep_not_nand_ancestor.h"
#include "resource_consumption.h"
#include "configurable_mutation.h"


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/isa.h>
#include <ea/digital_evolution/spatial.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/meta_population.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>

using namespace ea;
using namespace boost::accumulators;


LIBEA_MD_DECL(LAST_TASK, "ea.ape.last_task", std::string);
LIBEA_MD_DECL(GERM_MUTATION_PER_SITE_P, "ea.ape.germ_mutation_per_site_p", double);
LIBEA_MD_DECL(TASK_LETHALITY_PROB, "ea.ape.task_lethality_prob", double);
LIBEA_MD_DECL(EACH_TASK_THRESH, "ea.ape.each_task_thresh", double);
LIBEA_MD_DECL(NOT_LETHALITY_PROB, "ea.ape.not_lethality_prob", double);
LIBEA_MD_DECL(NAND_LETHALITY_PROB, "ea.ape.nand_lethality_prob", double);
LIBEA_MD_DECL(ORNOT_LETHALITY_PROB, "ea.ape.ornot_lethality_prob", double);
LIBEA_MD_DECL(NOT_AGE, "ea.ape.not_age", double);
LIBEA_MD_DECL(NAND_AGE, "ea.ape.nand_age", double);
LIBEA_MD_DECL(ORNOT_AGE, "ea.ape.ornot_age", double);

/*! Divide this organism's memory between parent and offspring.
 
 Instructions from the beginning of the organism's memory to the current
 position of the read head are preserved for the parent, while instructions
 between the read head and the write head are split off to form the
 offspring's genome; the offspring is then ``born.''
 */
DIGEVO_INSTRUCTION_DECL(h_divide_soft_parent_reset) {
    if(hw.age() >= (0.8 * hw.original_size())) {
        typename Hardware::representation_type& r=hw.repr();
        
        // Check to see if the offspring would be a good length.
        int divide_pos = hw.getHeadLocation(Hardware::RH);
        int extra_lines = r.size() - hw.getHeadLocation(Hardware::WH);
        
        int child_size = r.size() - divide_pos - extra_lines;
        int parent_size = r.size() - child_size - extra_lines;
        double ratio = 2.0;
        
        if ((child_size < (hw.original_size()/ratio)) ||
            (child_size > (hw.original_size()*ratio)) ||
            (parent_size < (hw.original_size()/ratio)) ||
            (parent_size > (hw.original_size()*ratio))){
            // fail!
            return;
        }
        
        
        typename Hardware::representation_type::iterator f=r.begin(),l=r.begin();
        std::advance(f, hw.getHeadLocation(Hardware::RH));
        std::advance(l, hw.getHeadLocation(Hardware::WH));
        typename Hardware::representation_type offr(f, l);
        
        
        r.resize(parent_size);
        replicate(p, offr, ea);
        hw.replicated_soft_reset();
        
    }
}





/*! This event triggers a task to have a lethal consequence.
 */

template <typename EA>
struct task_lethality : task_performed_event<EA> {
    
    task_lethality(EA& ea) : task_performed_event<EA>(ea) {
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
struct task_first_age : task_performed_event<EA> {
    task_first_age(EA& ea) : task_performed_event<EA>(ea) {
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
        } else if (t == "ornot") {
            if (!exists<ORNOT_AGE>(ind)) {
                put<ORNOT_AGE>(ind.hw().age(), ind) ;
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
        .add_field("not age")
        .add_field("nand age")
        .add_field("ornot age")
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
            double t_ornot = 0.0;
            double t_ornot_count = 0.0;
            
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j){
                    typename EA::individual_type::individual_type& ind=**j;
                    
                    if (exists<NOT_AGE>(ind)) {
                        t_not += get<NOT_AGE>(ind);
                        ++t_not_count;
                    }
                    if (exists<NAND_AGE>(ind)) {
                        t_nand += get<NAND_AGE>(ind);
                        ++t_nand_count;
                    }
                    if (exists<ORNOT_AGE>(ind)) {
                        t_ornot += get<ORNOT_AGE>(ind);
                        ++t_ornot_count;
                    }
                    
                }
            }
            
            if (t_not_count) { t_not /= t_not_count; }
            if (t_nand_count) { t_nand /= t_nand_count; }
            if (t_ornot_count) { t_ornot /= t_ornot_count; }
            
            _df.write(ea.current_update())
            .write(t_not)
            .write(t_nand)
            .write(t_ornot)
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

#endif