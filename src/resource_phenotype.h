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

LIBEA_MD_DECL(TASK_NOT_REWARD, "ea.not_reward", double);


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
        
        if (get<TASK_NOT_REWARD>(ind, 0.0) >= 50) { return; }

        
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
        get<TASK_NOT_REWARD>(ind, 0.0) += res;
        get<TASK_NOT_REWARD>(ea, 0.0) += res;
        get<TASK_NOT>(ea, 0.0) += 1.0;
        get<TASK_NOT>(ind, 0.0) += 1.0;
        
    }
};

template <typename EA>
struct stripes : reaction_event<EA> {
    stripes(EA& ea) : reaction_event <EA>(ea) {
    }
    
    virtual ~stripes() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer
                            double r,
                            EA& ea) {
        
        
        if (get<TASK_NOT_REWARD>(ind, 0.0) >= 50) { return; } 

        // rather than using r, we are just assigning a value for task performance on the
        // basis of how many neighbors have performed the task (not)
        
         
        /*!
         Unit circle:
         3  |  2  |  1
         4  |  Or.|  0
         5  |  6  |  7
         */
                
        double res = 0.25;
        typedef typename EA::environment_type::iterator env_it;
        
        // north
        env_it n = ea.env().direction_neighbor(ind, 2, ea);
        if (n->occupied() && get<TASK_NOT>(*(n->inhabitant()),0.0)) { res += 0.5; }
        
        // south
        env_it n2 = ea.env().direction_neighbor(ind, 6, ea);
        if (n2->occupied() && get<TASK_NOT>(*(n2->inhabitant()),0.0)) { res += 0.5; }
        
        // spin around the remainder of the neighbors...
        for (int i=0; i<=6; i+=2) {
            if (i==2 || i ==6 || res<0.25) { continue; }

            env_it n3 = ea.env().direction_neighbor(ind, i, ea);
            if(n3->occupied() && get<TASK_NOT>(*(n3->inhabitant()),0.0)) {
                res -= 0.25;
            }
        }
        
        get<SAVED_RESOURCES>(ind, 0.0) += res;
        get<TASK_NOT_REWARD>(ind, 0.0) += res;
        get<TASK_NOT_REWARD>(ea, 0.0) += res;
        get<TASK_NOT>(ea, 0.0) += 1.0;
        get<TASK_NOT>(ind, 0.0) += 1.0;
        
    }
};


/*! Prints information about the mean amount of reward organisms got for their tasks. 
 (higher means that 'spots' are appearing more frequently.
 */


template <typename EA>
struct reward_tracking : end_of_update_event<EA> {
    reward_tracking(EA& ea) : end_of_update_event<EA>(ea), _df("ps.dat") {
        _df.add_field("update")
        .add_field("sub_pop_size")
        .add_field("pop_size")
        .add_field("mean_reward");
        
    }
    
    //! Destructor.
    virtual ~reward_tracking() {
    }
    
    //! Track how many task-switches are being performed!
    virtual void operator()(EA& ea) {
        if ((ea.current_update() % 100) == 0) {
            double org = 0;
            
            int sub_pop_size = 0;
            double rew = 0.0;
            double not_count = 0.0;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                ++sub_pop_size;
                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j){
                    rew += get<TASK_NOT_REWARD>(*i,0.0);
                    not_count += get<TASK_NOT>(*i,0.0);
                    typename EA::individual_type::individual_type& ind=**j;
                    if (ind.alive()) {
                        ++org;
                    }
                }
            }
            if (not_count) {
                rew /= not_count; 
            }
            _df.write(ea.current_update())
            .write(sub_pop_size)
            .write(org)
            .write(rew)
            .endl();
        }
        
    }
    datafile _df;
    
};




#endif
