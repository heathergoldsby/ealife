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
#include <ea/datafiles/generation_priority.h>
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
                
                // grab a copy of the founder!
                typename EA::individual_type::individual_type prop = (*i).founder();
                
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
