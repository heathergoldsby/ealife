//
//  resource_phenotype.h
//  ealife
//
//  Created by Heather Goldsby on 10/4/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_RESOURCE_PHENOTYPE_H_
#define _EALIFE_RESOURCE_PHENOTYPE_H_

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



/*! Spots: Controls how many resources an organism recieves on the basis of how many neighbors
    have done a task.
 
    Also tracks resources at the subpop level.
 
 */

template <typename EA>
struct spots : reaction_event<EA> {
    spots(EA& ea) : reaction_event <EA>(ea) {
    }
    
    virtual ~spots() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer
                            double r,
                            EA& ea) {
        
        // rather than using r, we are just assigning a value for task performance on the
        // basis of how many neighbors have performed the task (not)
        
        // check neigbhors
        double neighbor_task_count = 0.0;
        double potential_neighbors = 0.0;
        
        // spin around the circle (find even locations), check if they performed NOT
        for (int i=0; i<=6; i+=2) {
            typedef typename EA::environment_type::iterator env_it;
            env_it n = ea.env().direction_neighbor(ind, i, ea);
            ++potential_neighbors;
            if(n->occupied() && get<TASK_NOT>(*(n->inhabitant()),0.0)) {
                    neighbor_task_count++;
            }
        }
      
        
        
        double res = (potential_neighbors - neighbor_task_count + 1.0)/potential_neighbors;
        
        get<SAVED_RESOURCES>(ind, 0.0) += res;
        std::string t = task->name();
        if (t == "not") { get<TASK_NOT>(ea,0.0) += 1.0; get<TASK_NOT>(ind,0.0) += 1.0; }
        
    }
};


#endif
