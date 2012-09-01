/* main.cpp
 * 
 * This file is part of EALife.
 * 
 * Copyright 2012 David B. Knoester, Heather J. Goldsby.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "gls.h"



/*! Triggers a task having a mutagenic effect on a *GROUP*.
 This is a control, since it eliminates the ability of organism's to
 create pseudo-soma. 
 */

template <typename EA>
struct task_mutagenesis_group : task_performed_event<EA> {
    
    task_mutagenesis_group(EA& ea) : task_performed_event<EA>(ea) {
    }
    
    virtual ~task_mutagenesis_group() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer   
                            double r, // amount of resource consumed
                            EA& ea) {
        
        double mult = get<TASK_MUTATION_MULT>(*task);
        double prob = get<TASK_MUTATION_PER_SITE_P>(ea) * mult;
        if (prob > 0) {
            configurable_per_site m(prob); 
            // grab a random individual from the population... (sacrificial lamb)
            typename EA::individual_type& ind_mut = **ea.rng().choice(ea.population().begin(), ea.population().end());
            mutate(ind_mut,m,ea);
            get<WORKLOAD>(ind_mut,0.0) += 1.0;
        }
    }
};


/*! Triggers a task having a mutagenic effect on a given category (germ/soma)
 within a group. This is a control, since it eliminates the ability of 
 organism's to create pseudo-soma. 
 */

template <typename EA>
struct task_mutagenesis_category : task_performed_event<EA> {
    
    task_mutagenesis_category(EA& ea) : task_performed_event<EA>(ea) {
    }
    
    virtual ~task_mutagenesis_category() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer   
                            double r, // amount of resource consumed
                            EA& ea) {
        
        double mult = get<TASK_MUTATION_MULT>(*task);
        double prob = get<TASK_MUTATION_PER_SITE_P>(ea) * mult;
        std::random_shuffle(ea.population().begin(), ea.population().end(), ea.rng());
        bool gs = get<GERM_STATUS>(ind, true);

        if (prob > 0) {
            configurable_per_site m(prob); 
            
            for(typename EA::population_type::iterator j=ea.population().begin(); j!=ea.population().end(); ++j) {
                typename EA::individual_type& org=**j;
                if (get<GERM_STATUS>(org, true) == gs) {            
                    // grab a random individual of the same category from the population... (sacrificial lamb)
                    mutate(org,m,ea);
                    get<WORKLOAD>(org,0.0) += 1.0;
                    return;
                }
            }
        }
    }
};




//! Configuration object for an EA.
template <typename EA>
struct gls_configuration : public abstract_configuration<EA> {
    
    typedef typename EA::tasklib_type::task_ptr_type task_ptr_type;
    typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
    
    
    //! Called as the final step of EA construction.
    void construct(EA& ea) {
        using namespace ea::instructions;
        append_isa<nop_a>(0,ea); // 0
        append_isa<nop_b>(0,ea);
        append_isa<nop_c>(0,ea);
        append_isa<nop_x>(ea);
        append_isa<mov_head>(ea);
        append_isa<if_label>(ea); //5
        append_isa<h_search>(ea);
        append_isa<nand>(ea);
        append_isa<push>(ea);
        append_isa<pop>(ea);
        append_isa<swap>(ea);//10
        append_isa<inc>(ea);
        append_isa<dec>(ea);
        append_isa<tx_msg>(ea); 
        append_isa<rx_msg>(ea); //15
        append_isa<bc_msg>(ea);
        append_isa<rotate>(ea);
        append_isa<rotate_cw>(ea);
        append_isa<rotate_ccw>(ea);
        append_isa<if_less>(ea); //20
        append_isa<h_alloc>(ea);             
        append_isa<h_copy>(ea);
        append_isa<h_divide>(ea);
        append_isa<input>(ea);
        append_isa<output>(ea);//25
        append_isa<become_soma>(ea);
        append_isa<if_germ>(ea);
        append_isa<if_soma>(ea);
        append_isa<donate_res_to_group>(ea);
        append_isa<get_xy>(ea);
        
        // Add tasks
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", ea);
        task_ptr_type task_and = make_task<tasks::task_and,catalysts::additive<0> >("and", ea);
        task_ptr_type task_ornot = make_task<tasks::task_ornot,catalysts::additive<0> >("ornot", ea);
        task_ptr_type task_or = make_task<tasks::task_or,catalysts::additive<0> >("or", ea);
        task_ptr_type task_andnot = make_task<tasks::task_andnot,catalysts::additive<0> >("andnot", ea);
        task_ptr_type task_nor = make_task<tasks::task_nor,catalysts::additive<0> >("nor", ea);
        task_ptr_type task_xor = make_task<tasks::task_xor,catalysts::additive<0> >("xor", ea);
        task_ptr_type task_equals = make_task<tasks::task_equals,catalysts::additive<0> >("equals", ea);
        
        put<TASK_MUTATION_MULT>(0, *task_not);
        put<TASK_MUTATION_MULT>(1, *task_nand);
        put<TASK_MUTATION_MULT>(1, *task_and);
        put<TASK_MUTATION_MULT>(1, *task_ornot);
        put<TASK_MUTATION_MULT>(1, *task_or);
        put<TASK_MUTATION_MULT>(1, *task_andnot);
        put<TASK_MUTATION_MULT>(1, *task_nor);
        put<TASK_MUTATION_MULT>(1, *task_xor);
        put<TASK_MUTATION_MULT>(1, *task_equals);
        
        resource_ptr_type resA = make_resource("resA", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resB = make_resource("resB", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resC = make_resource("resC", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resD = make_resource("resD", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resE = make_resource("resE", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resF = make_resource("resF", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resG = make_resource("resG", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resH = make_resource("resH", 100.0, 1.0, 0.01, 0.05, ea);
        resource_ptr_type resI = make_resource("resI", 100.0, 1.0, 0.01, 0.05, ea);
        
        task_not->consumes(resA);
        task_nand->consumes(resB);
        task_and->consumes(resC);
        task_ornot->consumes(resD);
        task_or->consumes(resE);
        task_andnot->consumes(resF);
        task_nor->consumes(resG);
        task_xor->consumes(resH);
        task_equals->consumes(resI);
        
        add_event<task_mutagenesis_category>(this,ea);
        add_event<gs_inherit_event>(this,ea);
        add_event<task_resource_consumption>(this,ea);
        
    }
    
    //! Called to generate the initial EA population.
    void initial_population(EA& ea) {
        generate_ancestors(selfrep_not_ancestor(), 1, ea);
    }
};


/*! Artificial life simulation definition.
 */
typedef digital_evolution<
gls_configuration, spatial, empty_neighbor, round_robin
> ea_type;

template <typename EA>
struct mp_configuration : public abstract_configuration<EA> {
};


//! Meta-population definition.
typedef meta_population<ea_type, mp_configuration> mea_type;




/*! 
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    
    
    virtual void gather_options() {
        add_option<SPATIAL_X>(this);
        add_option<SPATIAL_Y>(this);
        add_option<META_POPULATION_SIZE>(this);
        add_option<POPULATION_SIZE>(this);
        add_option<INITIAL_POPULATION_SIZE>(this);
        add_option<REPRESENTATION_SIZE>(this);
        add_option<SCHEDULER_TIME_SLICE>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_INSERTION_P>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<MUTATION_UNIFORM_INT_MIN>(this);
        add_option<MUTATION_UNIFORM_INT_MAX>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);        
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        
        // gls specific options
        add_option<TASK_MUTATION_PER_SITE_P>(this);
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        add_option<GROUP_REP_THRESHOLD>(this);
        
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        
        add_event<gls_replication>(this,ea);
        add_event<task_performed_tracking>(this,ea);

        
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
