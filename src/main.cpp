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
#include <boost/graph/adjacency_list.hpp>

#include <ea/artificial_life/artificial_life.h>
#include <ea/artificial_life/digital_evolution/hardware.h>
#include <ea/artificial_life/digital_evolution/isa.h>
#include <ea/artificial_life/spatial.h>
#include <ea/artificial_life/datafiles/reactions.h>
#include <ea/artificial_life/datafiles/generation_priority.h>
#include <ea/cmdline_interface.h>
#include <ea/meta_population.h>
#include <ea/selection/random.h>
#include <ea/mutation.h>

using namespace ea;


/*
 Move ancestor to selfrep_not_ancestor
 Double check group replication
 Track results
 */


/*! Generates a self-replicating ancestor that performs not.
 */
struct selfrep_not_ancestor {
    template <typename EA>
    typename EA::population_entry_type operator()(EA& ea) {
        typedef typename EA::representation_type representation_type;
        typename EA::individual_type ind;
        ind.name() = next<INDIVIDUAL_COUNT>(ea);
        
        representation_type& repr=ind.repr();
        repr.resize(get<REPRESENTATION_SIZE>(ea));
        std::fill(repr.begin(), repr.end(), 3);
        
        // Must use representation size of 100.
        assert(repr.size() == 100);
        
        repr[0] = 21; // h_alloc
        repr[1] = 2; // nopc
        repr[2] = 0; // nopa
        repr[3] = 6; // hsearch
        repr[4] = 2; // nopc
        repr[5] = 4; // movhead
        
        // not
        repr[24] = 24; // input
        repr[25] = 24; // input
        repr[26] = 8; // push
        repr[27] = 2; // nopc
        repr[28] = 9; // pop
        repr[29] = 7; // nand
        repr[30] = 25; //output
        repr[31] = 29; // donate_res_to_group
        
        repr[91] = 6; // hsearch
        repr[92] = 22; // hcopy
        repr[93] = 2; // nopc
        repr[94] = 0; // nopa
        repr[95] = 5; // iflabel
        repr[96] = 23; // hdivide
        repr[97] = 4; // movhead
        repr[98] = 0; // nopa
        repr[99] = 1; // nopb
        
        ind.hw().initialize();
        
        return make_population_entry(ind, ea);
    }
};



/*! Artificial life simulation definition.
 */
typedef artificial_life<
hardware, isa, spatial, empty_neighbor, round_robin, 
mutation::per_site<mutation::uniform_integer>, task_library, organism, population, 
alife_population<selfrep_not_ancestor>
> al_type;


typedef meta_population<al_type> mea_type;

LIBEA_MD_DECL(GERM_STATUS, "ea.gls.germ_status", bool);
LIBEA_MD_DECL(GROUP_RESOURCE_UNITS, "ea.gls.group_resource_units", double);
LIBEA_MD_DECL(SAVED_RESOURCES, "ea.gls.organism_saved_resources", double);
LIBEA_MD_DECL(TASK_MUTATION_PER_SITE_P, "ea.gls.task_mutation_per_site_p", double);
LIBEA_MD_DECL(GERM_MUTATION_PER_SITE_P, "ea.gls.germ_mutation_per_site_p", double);
LIBEA_MD_DECL(GROUP_REP_THRESHOLD, "ea.gls.group_rep_threshold", double);



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
struct germline_replication : end_of_update_event<EA> {
    //! Constructor.
    germline_replication(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~germline_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            if (exists<GROUP_RESOURCE_UNITS>(*i) && 
                (get<GROUP_RESOURCE_UNITS>(*i) > get<GROUP_REP_THRESHOLD>(*i))){
                
                typename EA::individual_ptr_type p(new typename EA::individual_type());
                
                // grab a copy of the first individual: 
                typename EA::individual_type::individual_type germ;
                bool germ_present = false;
                
                // If so, setup a new replicate pop.
                // Find a germ...
                std::random_shuffle(i->population().begin(), i->population().end(), ea.rng());
                
                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                    typename EA::individual_type::individual_type& org=**j;
                    if (get<GERM_STATUS>(org)) {
                        germ = **j;
                        germ_present = true;
                        break;
                    }
                }
                
                if (!germ_present) continue;
                
                // setup the population (really, an ea):
                p->md() = ea.md();
                p->rng().reset(ea.rng()(std::numeric_limits<int>::max()));
                p->initialize();
                
                // mutate it:
                per_site_germ m; 
                mutate(germ,m,*i);
                
                // and fill up the offspring population with copies of the germ:
                
                typename EA::individual_type::individual_ptr_type o=make_population_entry(germ,*p);
                p->population().push_back(o);
                p->env().insert(o);
                
                offspring.push_back(p);
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

// Germ instructions!

/*! Mark an organism as soma.
 */

DIGEVO_INSTRUCTION_DECL(become_soma) {
    put<GERM_STATUS>(false,*p);
    return 1;
}


/*! Execute the next instruction if the organism is marked as germ.
 */

DIGEVO_INSTRUCTION_DECL(if_germ) {
    if(!get<GERM_STATUS>(*p)) {
        hw.advanceHead(Hardware::IP);
    }
    return 1;
}


/*! Execute the next instruction if the organism is marked as soma.
 */
DIGEVO_INSTRUCTION_DECL(if_soma){
    if(get<GERM_STATUS>(*p)) {
        hw.advanceHead(Hardware::IP);
    }
    return 1;
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
    return 1;
}

template <typename EA>
void add_instructions(EA& ea) {
    using namespace ea::instructions;
    append_isa<nop_a>(ea); // 0
    append_isa<nop_b>(ea);
    append_isa<nop_c>(ea);
    append_isa<nop_x>(ea);
    append_isa<mov_head>(ea);
    append_isa<if_label>(ea); //5
    append_isa<h_search>(ea);
    append_isa<nand>(ea);
    append_isa<push>(ea);
    append_isa<pop>(ea);
    append_isa<swap>(ea);//10
    append_isa<latch_ldata>(ea);
    append_isa<inc>(ea);
//    append_isa<repro>(ea); //13
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
    
}



/*! 
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    typedef typename EA::individual_type::tasklib_type::task_ptr_type task_ptr_type;
    typedef typename EA::individual_type::environment_type::resource_ptr_type resource_ptr_type;
    
    virtual void postinitialization(EA& ea) {
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {            
            add_instructions(**i);
        }
    }
    
    virtual void configure(EA& ea) {
        
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {                    
            task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", **i);
            task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", **i);
            task_ptr_type task_and = make_task<tasks::task_and,catalysts::additive<0> >("and", **i);
            task_ptr_type task_ornot = make_task<tasks::task_ornot,catalysts::additive<0> >("ornot", **i);
            task_ptr_type task_or = make_task<tasks::task_or,catalysts::additive<0> >("or", **i);
            task_ptr_type task_andnot = make_task<tasks::task_andnot,catalysts::additive<0> >("andnot", **i);
            task_ptr_type task_nor = make_task<tasks::task_nor,catalysts::additive<0> >("nor", **i);
            task_ptr_type task_xor = make_task<tasks::task_xor,catalysts::additive<0> >("xor", **i);
            task_ptr_type task_equals = make_task<tasks::task_equals,catalysts::additive<0> >("equals", **i);
            
            resource_ptr_type resA = make_resource("resA", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resB = make_resource("resB", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resC = make_resource("resC", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resD = make_resource("resD", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resE = make_resource("resE", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resF = make_resource("resF", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resG = make_resource("resG", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resH = make_resource("resH", 100.0, 1.0, 0.01, 0.05, **i);
            resource_ptr_type resI = make_resource("resI", 100.0, 1.0, 0.01, 0.05, **i);
            
            task_not->consumes(resA);
            task_nand->consumes(resB);
            task_and->consumes(resC);
            task_ornot->consumes(resD);
            task_or->consumes(resE);
            task_andnot->consumes(resF);
            task_nor->consumes(resG);
            task_xor->consumes(resH);
            task_equals->consumes(resI);
        }
        
    }
    
    /*
    typedef typename EA::tasklib_type::task_ptr_type task_ptr_type;
    typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
    
    virtual void postinitialization(EA& ea) {
        add_instructions(ea);
        
    }
    
    virtual void configure(EA& ea) {
        
        //for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {                    
            task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<2> >("not", ea);
            task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<2> >("nand", ea);
            task_ptr_type task_and = make_task<tasks::task_and,catalysts::additive<2> >("and", ea);
            task_ptr_type task_ornot = make_task<tasks::task_ornot,catalysts::additive<2> >("ornot", ea);
            task_ptr_type task_or = make_task<tasks::task_or,catalysts::additive<2> >("or", ea);
            task_ptr_type task_andnot = make_task<tasks::task_andnot,catalysts::additive<2> >("andnot", ea);
            task_ptr_type task_nor = make_task<tasks::task_nor,catalysts::additive<2> >("nor", ea);
            task_ptr_type task_xor = make_task<tasks::task_xor,catalysts::additive<2> >("xor", ea);
            task_ptr_type task_equals = make_task<tasks::task_equals,catalysts::additive<2> >("equals", ea);
            
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
        //}
        
    }
     */

    
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

            add_event<datafiles::record_reactions_event>(this,*i);
            add_event<datafiles::generation_priority>(this,*i);
            add_event<task_mutagenesis>(this,*i);
            add_event<gs_inherit_event>(this,*i);
            add_event<task_resource_consumption>(this,*i);
        }*/
        add_event<germline_replication>(this,ea);

    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
