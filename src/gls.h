//
//  gls.h
//  ealife
//
//  Created by Heather Goldsby on 8/23/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_GLS_H_
#define _EALIFE_GLS_H_


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/algorithm/string.hpp>

#include "selfrep_not_ancestor.h"
#include "repro_not_ancestor.h"
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
LIBEA_MD_DECL(WORKLOAD, "ea.gls.workload", double);
LIBEA_MD_DECL(TASK_MUTATION_MULT, "ea.gls.task_mutation_mult", double);

LIBEA_MD_DECL(TASK_NOT, "ea.not", double);
LIBEA_MD_DECL(TASK_NAND, "ea.nand", double);
LIBEA_MD_DECL(TASK_AND, "ea.and", double);
LIBEA_MD_DECL(TASK_ORNOT, "ea.ornot", double);
LIBEA_MD_DECL(TASK_OR, "ea.or", double);
LIBEA_MD_DECL(TASK_ANDNOT, "ea.andnot", double);
LIBEA_MD_DECL(TASK_NOR, "ea.nor", double);
LIBEA_MD_DECL(TASK_XOR, "ea.xor", double);
LIBEA_MD_DECL(TASK_EQUALS, "ea.equals", double);


// Germ instructions!

/*! Mark an organism as soma.
 */

DIGEVO_INSTRUCTION_DECL(become_soma) {
    put<GERM_STATUS>(false,*p);
}


/*! Execute the next instruction if the organism is marked as germ.
 */

DIGEVO_INSTRUCTION_DECL(if_germ) {
    if(!get<GERM_STATUS>(*p,true)) {
        hw.advanceHead(Hardware::IP);
    }
}


/*! Execute the next instruction if the organism is marked as soma.
 */
DIGEVO_INSTRUCTION_DECL(if_soma){
    if(get<GERM_STATUS>(*p, true)) {
        hw.advanceHead(Hardware::IP);
    }
}

/*! Donate an organism's resources to the group. 
 */

DIGEVO_INSTRUCTION_DECL(donate_res_to_group){
    
    get<GROUP_RESOURCE_UNITS>(ea, 0.0) += get<SAVED_RESOURCES>(*p,0.0);
    put<SAVED_RESOURCES>(0,*p);
}


/*! Mutation - Per-site mutation at a configurable rate
 */
struct configurable_per_site {            
    typedef mutation::uniform_integer mutation_type;
    
    configurable_per_site(double prob) : _mp(prob) {
    }
    
    //! Iterate through all elements in the given representation, possibly mutating them.
    template <typename Representation, typename EA>
    void operator()(Representation& repr, EA& ea) {
        for(typename Representation::iterator i=repr.begin(); i!=repr.end(); ++i){
            if(ea.rng().p(_mp)) {
                _mt(repr, i, ea);
            }
        }
    }
    
    mutation_type _mt;
    double _mp; //! Mutation probability
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
        
        get<GERM_STATUS>(offspring, true) = get<GERM_STATUS>(ind(parents.begin(),ea), true);
    }
};


/*! Triggers a task having a mutagenic effect on a organism.
 One mutagenic rate for all tasks other than NOT.
 */

template <typename EA>
struct task_mutagenesis : task_performed_event<EA> {
    
    task_mutagenesis(EA& ea) : task_performed_event<EA>(ea) {
    }
    
    virtual ~task_mutagenesis() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer
                            double r, // amount of resource consumed
                            EA& ea) {

        double mult = get<TASK_MUTATION_MULT>(*task);
        double prob = get<TASK_MUTATION_PER_SITE_P>(ea) * mult;
        if (prob > 0) {
            configurable_per_site m(prob); 
            mutate(ind,m,ea);
            get<WORKLOAD>(ind,0.0) += mult;
        }
    }
};



/*! Triggers a task having a mutagenic effect on a *GROUP*.
 This is a control, since it eliminates the ability of organism's to
 create pseudo-soma. There is one mutagenic rate for all tasks other 
 than NOT.
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

/*! Tracks an organism's resources and tasks. 
 */

template <typename EA>
struct task_resource_consumption : task_performed_event<EA> {
    task_resource_consumption(EA& ea) : task_performed_event<EA>(ea) {
    }
    
    virtual ~task_resource_consumption() { }
    virtual void operator()(typename EA::individual_type& ind, // individual
                            typename EA::tasklib_type::task_ptr_type task, // task pointer   
                            double r,
                            EA& ea) {
        get<SAVED_RESOURCES>(ind, 0.0) += r;
        std::string t = task->name();
        if (t == "not") { get<TASK_NOT>(ea,0.0) += 1.0; }
        else if (t == "nand") {get<TASK_NAND>(ea,0.0) += 1.0;}
        else if (t == "and") {get<TASK_AND>(ea,0.0) += 1.0;}
        else if (t == "ornot") {get<TASK_ORNOT>(ea,0.0) += 1.0;}
        else if (t == "or") {get<TASK_OR>(ea,0.0) += 1.0;}
        else if (t == "andnot") {get<TASK_ANDNOT>(ea,0.0) += 1.0;}
        else if (t == "nor") {get<TASK_NOR>(ea,0.0) += 1.0;}
        else if (t == "xor") {get<TASK_XOR>(ea,0.0) += 1.0;}
        else if (t == "equals") {get<TASK_EQUALS>(ea,0.0) += 1.0;}
        
    }
};

/*! Prints information about the aggregate task performance of the group.
 */


template <typename EA>
struct task_performed_tracking : end_of_update_event<EA> {
    task_performed_tracking(EA& ea) : end_of_update_event<EA>(ea), _df("tasks.dat") { 
        _df.add_field("update")
        .add_field("not")
        .add_field("nand")
        .add_field("and")
        .add_field("ornot")
        .add_field("or")
        .add_field("andnot")
        .add_field("nor")
        .add_field("xor")
        .add_field("equals");
    }
    
    //! Destructor.
    virtual ~task_performed_tracking() {
    }
    
    //! Track resources!
    virtual void operator()(EA& ea) {
        if ((ea.current_update() % 100) == 0) {
            int t_not = 0;
            int t_nand = 0;
            int t_and = 0;
            int t_ornot = 0; 
            int t_or = 0;
            int t_andnot = 0;
            int t_nor = 0;
            int t_xor = 0;
            int t_equals = 0;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                t_not += get<TASK_NOT>(*i, 0.0);
                t_nand += get<TASK_NAND>(*i, 0.0);
                t_and += get<TASK_AND>(*i, 0.0);
                t_ornot += get<TASK_ORNOT>(*i, 0.0);
                t_or += get<TASK_OR>(*i, 0.0);
                t_andnot += get<TASK_ANDNOT>(*i, 0.0);
                t_nor += get<TASK_NOR>(*i, 0.0);
                t_xor += get<TASK_XOR>(*i, 0.0);
                t_equals += get<TASK_EQUALS>(*i, 0.0);
            }
            
            _df.write(ea.current_update())
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
struct gls_replication : end_of_update_event<EA> {
    //! Constructor.
    gls_replication(EA& ea) : end_of_update_event<EA>(ea), _df("gls.dat") {
        _df.add_field("update")
        .add_field("mean_germ_num")
        .add_field("mean_pop_num")
        .add_field("mean_germ_percent")
        .add_field("mean_germ_workload")
        .add_field("mean_germ_workload_var")
        .add_field("mean_soma_workload")
        .add_field("mean_soma_workload_var")
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
                accumulator_set<double, stats<tag::mean, tag::variance> > germ_workload_acc; 
                accumulator_set<double, stats<tag::mean, tag::variance> > soma_workload_acc; 
                
                
                
                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                    typename EA::individual_type::individual_type& org=**j;
                    if (get<GERM_STATUS>(org, true)) {
                        ++germ_count;
                        germ_workload_acc(get<WORKLOAD>(org, 0.0));
                        if (!germ_present){
                            germ = org;
                            // Makes sure that we keep the size of the organism and discard its in-memory offspring
                            germ.repr().resize(org.hw().original_size());
                            germ.hw().initialize();
                            germ_present = true;
                        }
                    } else {
                        soma_workload_acc(get<WORKLOAD>(org, 0.0));
                    }
                    ++pop_count;
                }
                
                if (!germ_present) continue;
                
                pop_num.push_back(pop_count);
                germ_num.push_back(germ_count);
                germ_percent.push_back(germ_count/((double) i->population().size())*100.0); 
                germ_workload.push_back(mean(germ_workload_acc));
                germ_workload_var.push_back(variance(germ_workload_acc));
                
                if (germ_count != pop_count) {
                    soma_workload.push_back(mean(soma_workload_acc));
                    soma_workload_var.push_back(variance(soma_workload_acc));
                } else {
                    soma_workload.push_back(0);
                    soma_workload_var.push_back(0);
                }
                
                ++num_rep;
                
                
                if (germ_num.size() > 100) {
                    germ_num.pop_front();
                    germ_percent.pop_front();
                    pop_num.pop_front();
                    germ_workload.pop_front();
                    germ_workload_var.pop_front();
                    soma_workload.pop_front();
                    soma_workload_var.pop_front();
                }
                
                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p(new typename EA::individual_type());
                
                p->md() = ea.md();
                p->rng().reset(ea.rng()(std::numeric_limits<int>::max()));
                p->initialize();
                
                // mutate it:
                configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea)); 
                mutate(germ,m,*p);
                
                // and fill up the offspring population with copies of the germ:
                typename EA::individual_type::individual_ptr_type o=p->make_individual(germ.repr());
                p->append(o);
                offspring.push_back(p);
                
                // reset resource units
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
                .write(std::accumulate(germ_workload.begin(), germ_workload.end(), 0.0)/germ_workload.size())
                .write(std::accumulate(germ_workload_var.begin(), germ_workload_var.end(), 0.0)/germ_workload.size())
                .write(std::accumulate(soma_workload.begin(), soma_workload.end(), 0.0)/soma_workload.size())
                .write(std::accumulate(soma_workload_var.begin(), soma_workload_var.end(), 0.0)/soma_workload.size())
                .write(num_rep)
                .endl();
                num_rep = 0;
            } else {
                _df.write(ea.current_update())
                .write(0)
                .write(0)
                .write(0)
                .write(0)
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
    std::deque<double> germ_workload;
    std::deque<double> germ_workload_var;
    std::deque<double> soma_workload;
    std::deque<double> soma_workload_var;
    int num_rep;
    
    
};


#endif
