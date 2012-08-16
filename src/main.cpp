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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "selfrep_not_ancestor.h"
#include "repro_not_ancestor.h"
#include <ea/artificial_life.h>
#include <ea/artificial_life/hardware.h>
#include <ea/artificial_life/isa.h>
#include <ea/artificial_life/spatial.h>
#include <ea/datafiles/reactions.h>
#include <ea/datafiles/generation_priority.h>
#include <ea/cmdline_interface.h>
#include <ea/meta_population.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>

using namespace ea;
using namespace boost::accumulators;


/*
 Track results - tasks?
 Fixed inputs?
 
 */

LIBEA_MD_DECL(GERM_STATUS, "ea.gls.germ_status", bool);
LIBEA_MD_DECL(GROUP_RESOURCE_UNITS, "ea.gls.group_resource_units", double);
LIBEA_MD_DECL(SAVED_RESOURCES, "ea.gls.organism_saved_resources", double);
LIBEA_MD_DECL(TASK_MUTATION_PER_SITE_P, "ea.gls.task_mutation_per_site_p", double);
LIBEA_MD_DECL(GERM_MUTATION_PER_SITE_P, "ea.gls.germ_mutation_per_site_p", double);
LIBEA_MD_DECL(GROUP_REP_THRESHOLD, "ea.gls.group_rep_threshold", double);


// Germ instructions!

/*! Mark an organism as soma.
 */

DIGEVO_INSTRUCTION_DECL(become_soma) {
    put<GERM_STATUS>(false,*p);
}


/*! Execute the next instruction if the organism is marked as germ.
 */

DIGEVO_INSTRUCTION_DECL(if_germ) {
    if(!get<GERM_STATUS>(*p)) {
        hw.advanceHead(Hardware::IP);
    }
}


/*! Execute the next instruction if the organism is marked as soma.
 */
DIGEVO_INSTRUCTION_DECL(if_soma){
    if(get<GERM_STATUS>(*p)) {
        hw.advanceHead(Hardware::IP);
    }
}

/*! Donate an organism's resources to the group. 
 */

DIGEVO_INSTRUCTION_DECL(donate_res_to_group){
    if(exists<SAVED_RESOURCES>(ind(p,ea))) {
        double group_res = 0.0;
        if (exists<GROUP_RESOURCE_UNITS>(ea)) {
            group_res = get<GROUP_RESOURCE_UNITS>(ea); 
        }
        group_res += get<SAVED_RESOURCES>(*p); 
        put<GROUP_RESOURCE_UNITS>(group_res,ea);
        put<SAVED_RESOURCES>(0,*p);
    }
}



/*! Mutation - Per-site mutation at a gls rate.
 */
struct per_site_gls {            
    typedef mutation::uniform_integer mutation_type;
    
    //! Iterate through all elements in the given representation, possibly mutating them.
    template <typename Representation, typename EA>
    void operator()(Representation& repr, EA& ea) {
        const double per_site_p=get<TASK_MUTATION_PER_SITE_P>(ea);
        for(typename Representation::iterator i=repr.begin(); i!=repr.end(); ++i){
            if(ea.rng().p(per_site_p)) {
                _mt(repr, i, ea);
            }
        }
    }
    
    mutation_type _mt;
};

struct per_site_germ {            
    typedef mutation::uniform_integer mutation_type;
    
    //! Iterate through all elements in the given representation, possibly mutating them.
    template <typename Representation, typename EA>
    void operator()(Representation& repr, EA& ea) {
        const double per_site_p=get<GERM_MUTATION_PER_SITE_P>(ea);
        for(typename Representation::iterator i=repr.begin(); i!=repr.end(); ++i){
            if(ea.rng().p(per_site_p)) {
                _mt(repr, i, ea);
            }
        }
    }
    
    mutation_type _mt;
};

// Events!


/*! An organism inherits its parent's germ/soma status. If it is undefined, 
 then it is set to germ.
 */
template <typename EA>
struct gs_inherit_event : inheritance_event<EA> {
    
    //! Constructor.
    gs_inherit_event(EA& ea) : inheritance_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~gs_inherit_event() {
    }
    
    /*! Called for every inheritance event. We are using the germ/soma status 
     of the first parent
     */
    virtual void operator()(typename EA::population_type& parents,
                            typename EA::individual_type& offspring,
                            EA& ea) {
        if(!exists<GERM_STATUS>(ind(parents.begin(),ea))) {
            put<GERM_STATUS>(true,offspring);
        } else {
            bool parent_status = get<GERM_STATUS>(ind(parents.begin(),ea)); 
            put<GERM_STATUS>(parent_status,offspring);
        }
        
    }
};

/*! Triggers a task having a mutagenic effect on a organism.
 */

template <typename EA>
struct task_mutagenesis : task_performed_event<EA> {
    
    task_mutagenesis(EA& ea) : task_performed_event<EA>(ea) {
    }
    
    virtual ~task_mutagenesis() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            double r, // amount of resource consumed
                            const std::string& task_name, // task name
                            EA& ea) {
        // Check if task name should be mutated. 
        if (task_name != "not") {
            per_site_gls m; 
            mutate(ind,m,ea);
        }
    }
};

/*! Tracks an organism's resources.
 */

template <typename EA>
struct task_resource_consumption : task_performed_event<EA> {
    task_resource_consumption(EA& ea) : task_performed_event<EA>(ea) {
    }
    
    virtual ~task_resource_consumption() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            double r, // amount of resource consumed
                            const std::string& task_name, // task name
                            EA& ea) {
        double res_amount = r;
        if(exists<SAVED_RESOURCES>(ind)) {
            res_amount += get<SAVED_RESOURCES>(ind); 
        }
        put<SAVED_RESOURCES>(res_amount, ind);
    
    }
};


//! Performs group replication using germ lines.
template <typename EA>
struct gls_replication : end_of_update_event<EA> {
    //! Constructor.
    gls_replication(EA& ea) : end_of_update_event<EA>(ea), _df("gls.dat") {
        _df.add_field("update")
        .add_field("mean_germ_num")
        .add_field("mean_pop_num")
        .add_field("mean_germ_percent")
        .add_field("replication_count");
        
        num_rep = 0;
    }
    
    
    //! Destructor.
    virtual ~gls_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            if (exists<GROUP_RESOURCE_UNITS>(*i) && 
                (get<GROUP_RESOURCE_UNITS>(*i) > get<GROUP_REP_THRESHOLD>(*i))){
                
                
                // grab a copy of the first individual: 
                typename EA::individual_type::individual_type germ;
                bool germ_present = false;
                
                // If so, setup a new replicate pop.
                // Find a germ...
                std::random_shuffle(i->population().begin(), i->population().end(), ea.rng());
                
                int germ_count = 0;
                int pop_count = 0;
                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                    typename EA::individual_type::individual_type& org=**j;
                    if (get<GERM_STATUS>(org)) {
                        ++germ_count;
                        if (!germ_present){
                            germ = org;
                            // Makes sure that we keep the size of the organism and discard its in-memory offspring
                            germ.repr().resize(org.hw().original_size());
                            germ.hw().initialize();
                            germ_present = true;
                        }
                    } 
                    pop_count++;
                    
                }
                
                if (!germ_present) continue;
                
                pop_num.push_back(pop_count);
                germ_num.push_back(germ_count);
                germ_percent.push_back(germ_count/((double) i->population().size())*100.0); 
                ++num_rep;
                
                if (germ_num.size() > 100) {
                    germ_num.pop_front();
                    germ_percent.pop_front();
                    pop_num.pop_front();
                }
                
                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p(new typename EA::individual_type());
                
                p->md() = ea.md();
                p->rng().reset(ea.rng()(std::numeric_limits<int>::max()));
                p->initialize();
                
                // mutate it:
                per_site_germ m; 
                mutate(germ,m,*p);
                
                // and fill up the offspring population with copies of the germ:
                
                typename EA::individual_type::individual_ptr_type o=make_population_entry(germ,*p);
                p->population().push_back(o);
                p->env().insert(o);
                
                offspring.push_back(p);
                i->env().reset_resources();
                put<GROUP_RESOURCE_UNITS>(0,*i);
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
        
        //        assert(ea.population().size() == 10); 
        
        if ((ea.current_update() % 100) == 0) {
            if (germ_num.size() > 0) {
                _df.write(ea.current_update())
                .write(std::accumulate(germ_num.begin(), germ_num.end(), 0.0)/germ_num.size())
                .write(std::accumulate(pop_num.begin(), pop_num.end(), 0.0)/pop_num.size())
                .write(std::accumulate(germ_percent.begin(), germ_percent.end(), 0.0)/germ_percent.size())
                .write(num_rep)
                .endl();
                num_rep = 0;
            } else {
                _df.write(ea.current_update())
                .write(0)
                .write(0)
                .write(0)
                .write(0)
                .endl();     
            }
        }
    }
    datafile _df;    
    std::deque<double> germ_num; 
    std::deque<double> germ_percent;
    std::deque<double> pop_num;

    int num_rep;
    
    
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
        append_isa<repro>(ea);
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
        
        add_event<task_mutagenesis>(this,ea);
        add_event<gs_inherit_event>(this,ea);
        add_event<task_resource_consumption>(this,ea);
        
    }
    
    //! Called to generate the initial EA population.
    void initial_population(EA& ea) {
        alife_population<selfrep_not_ancestor> init;
        init(ea);
    }
};


/*! Artificial life simulation definition.
 */
typedef artificial_life<
gls_configuration, spatial, empty_neighbor, round_robin
> al_type;

namespace ea {
    template < > int al_type::alife_allocated = 0;
    template < > int al_type::alife_deallocated = 0;
    template < > int al_type::individual_type::org_allocated = 0;
    template < > int al_type::individual_type::org_deallocated = 0;
}
    
typedef meta_population<al_type> mea_type;




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
        /* for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {                    
         
         add_event<datafiles::generation_priority>(this,*i);
         add_event<datafiles::record_reactions_event>(this,ea);

         add_event<task_mutagenesis>(this,*i);
         add_event<gs_inherit_event>(this,*i);
         add_event<task_resource_consumption>(this,*i);
         }*/
        add_event<gls_replication>(this,ea);

        
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
