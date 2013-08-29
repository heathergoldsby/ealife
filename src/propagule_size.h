//
//  propagule_size.h
//  ealife
//
//  Created by Heather Goldsby on 5/22/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#ifndef ealife_propagule_size_h
#define ealife_propagule_size_h



#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/algorithm/string.hpp>

#include "selfrep_not_ancestor.h"
#include "repro_not_ancestor.h"
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

using namespace ealib;
using namespace boost::accumulators;


//! the number of organisms used to start a group
LIBEA_MD_DECL(PROP_SIZE, "ea.ps.propagule_size", double);
//! the composition of the propagule: 0: clonal, 1: randomly selected
LIBEA_MD_DECL(PROP_COMPOSITION, "ea.ps.propagule_composition", double);
//! the germ status of an org. true = germ, false = soma
LIBEA_MD_DECL(GERM_STATUS, "ea.gls.germ_status", bool);
//! the base number of units needed for propagule replication
LIBEA_MD_DECL(PROP_BASE_REP_UNITS, "ea.ps.prop_base_rep", double);
//! the extra number of units needed per cell for propagule replication
LIBEA_MD_DECL(PROP_CELL_REP_UNITS, "ea.ps.prop_cell_rep", double);
//! the actual number of orgs in the propagule
LIBEA_MD_DECL(ACTUAL_PROP_SIZE, "ea.ps.actual_prop_size", double);



//! epigenetic info
LIBEA_MD_DECL(EPIGENETIC_INFO, "ea.ts.epigenetic_info", double);

//! Set organism's epigenetic info
DIGEVO_INSTRUCTION_DECL(set_epigenetic_info){
    put<EPIGENETIC_INFO>(hw.getRegValue(hw.modifyRegister()), *p);
}

//! Get organism's epigenetic info
DIGEVO_INSTRUCTION_DECL(get_epigenetic_info){
    if(exists<EPIGENETIC_INFO>(*p)) {
        hw.setRegValue(hw.modifyRegister(), get<EPIGENETIC_INFO>(*p));
    }
}

//! Increment the propagule size suggested by the organism.
DIGEVO_INSTRUCTION_DECL(inc_propagule_size){
    get<PROPAGULE_SIZE>(*p,0)++;
}

//! Decrement the propagule size suggested by the organism.
DIGEVO_INSTRUCTION_DECL(dec_propagule_size){
    int prop_size = get<PROPAGULE_SIZE>(*p, 1);
    if (prop_size > 2) {
        get<PROPAGULE_SIZE>(*p)--;
    }
}

//! Get the propagule size suggested by the organism.
DIGEVO_INSTRUCTION_DECL(get_propagule_size){
    hw.setRegValue(hw.modifyRegister(), get<PROPAGULE_SIZE>(*p, 1));
}


//! Mark an organism as soma.
DIGEVO_INSTRUCTION_DECL(become_soma) {
    put<GERM_STATUS>(false,*p);
}


//! Execute the next instruction if the organism is marked as germ.
DIGEVO_INSTRUCTION_DECL(if_germ) {
    if(!get<GERM_STATUS>(*p,true)) {
        hw.advanceHead(Hardware::IP);
    }
}


//! Execute the next instruction if the organism is marked as soma.
DIGEVO_INSTRUCTION_DECL(if_soma){
    if(get<GERM_STATUS>(*p, true)) {
        hw.advanceHead(Hardware::IP);
    }
}



//! Performs group replication.
template <typename EA>
struct ps_size_propagule : end_of_update_event<EA> {
    //! Constructor.
    ps_size_propagule(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ps_size_propagule() {
    }
    
    //! Perform replication among populations.
    virtual void operator()(EA& ea) {
        // For each multicell:
        // (1) figure out its propagule size (mean of prop sizes of individuals)
        // (2) figure out how many resources it currently has
        // (3) can it replicate?
        
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            // Do not replicate if the 'founding org' is sterile.
            // Note may need to adjust this to include starting colony size.
            if (i->population().size() == get<ACTUAL_PROP_SIZE>(*i, 1.0)) continue;
            
            double group_res = get<GROUP_RESOURCE_UNITS>(*i,0.0);
            double prop_base = get<PROP_BASE_REP_UNITS>(*i,0.0);
            // If there aren't enough resources for the base replication rate, skip
            // assessing propagule size.
            if (group_res < get<PROP_BASE_REP_UNITS>(*i)) {
                continue;
            }
            
            double desired_prop_size = 0.0;
            int num_germ = 0;
            int num_org = 0;
            for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                typename EA::individual_type::individual_type& o=**j;
                
                ++num_org;
                desired_prop_size += get<PROPAGULE_SIZE>(o, 1);
                if (get<GERM_STATUS>(o,true)) {
                    ++num_germ;
                }
            }
            
            if (num_germ == 0) { break; }
            desired_prop_size = floor(desired_prop_size/num_org);
            if(desired_prop_size > num_germ) {
                desired_prop_size = num_germ;
            }
            if (desired_prop_size < 1) { desired_prop_size = 1; }
            
            
            
            // Can this multicell replicate
            double res_required = get<PROP_BASE_REP_UNITS>(*i) + (desired_prop_size * get<PROP_CELL_REP_UNITS>(*i));
            if (group_res < res_required) {
                continue; // multicell cannot replicate. continue
            }
            
            get<NUM_GROUP_REPLICATIONS>(ea,0) ++;
            // setup the population (really, an ea):
            typename EA::individual_ptr_type p = ea.make_individual();
            
            std::random_shuffle(i->population().begin(), i->population().end(), ea.rng());
            
            int p_size = 0;
            
            typename EA::individual_type::individual_type org;
            
            for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                typename EA::individual_type::individual_type& prop_org=**j;
                
                if (get<GERM_STATUS>(prop_org,true)) {
                    org = prop_org;
                    org.repr().resize(org.hw().original_size());
                    org.hw().initialize();
                
                    // mutate it:
                    configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
                    mutate(org,m,*p);
                    typename EA::individual_type::individual_ptr_type o=p->make_individual(org.repr());
                    if(exists<EPIGENETIC_INFO>(org)) {
                        put<EPIGENETIC_INFO>(get<EPIGENETIC_INFO>(org),*o);
                    }
                
                    p->append(o);
                    ++p_size;
                }
                
                if (p_size >= desired_prop_size) {
                    put<ACTUAL_PROP_SIZE>(p_size, *p);
                    break;
                }
            }
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



//! Performs group replication.
template <typename EA>
struct ts_replication_propagule : end_of_update_event<EA> {
    //! Constructor.
    ts_replication_propagule(EA& ea) : end_of_update_event<EA>(ea) {
    }
    
    
    //! Destructor.
    virtual ~ts_replication_propagule() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
        // See if any subpops have exceeded the resource threshold
        typename EA::population_type offspring;
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            // Do not replicate if the 'founding org' is sterile.
            if (i->population().size() < (get<PROP_SIZE>(*i) + 1)) continue;
            
            double x = get<GROUP_RESOURCE_UNITS>(*i,0.0);
            
            if (exists<GROUP_RESOURCE_UNITS>(*i) &&
                (get<GROUP_RESOURCE_UNITS>(*i) > get<GROUP_REP_THRESHOLD>(*i))){
                
                // the number of parents selected is the propagule size or 1, if
                // the propagule's composition is clonal.
                int num_parents = get<PROP_SIZE>(*i);
                if (get<PROP_COMPOSITION>(*i) == 0) {
                    num_parents = 1;
                }
                
                get<NUM_GROUP_REPLICATIONS>(ea,0) ++;

                
                // setup the population (really, an ea):
                typename EA::individual_ptr_type p = ea.make_individual();

                std::random_shuffle(i->population().begin(), i->population().end(), ea.rng());

                int p_size = 0;
                
                typename EA::individual_type::individual_type org;

                for(typename EA::individual_type::population_type::iterator j=i->population().begin(); j!=i->population().end(); ++j) {
                    typename EA::individual_type::individual_type& prop_org=**j;
                    org = prop_org;
                    org.repr().resize(org.hw().original_size());
                    org.hw().initialize();
                    
                    // mutate it:
                    configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
                    mutate(org,m,*p);
                    
                    if (get<PROP_COMPOSITION>(*i) == 0) {
                        for (int k=0; k<get<PROP_SIZE>(*i); ++k) {
                            // and fill up the offspring population with copies of the germ:
                            typename EA::individual_type::individual_ptr_type o=p->make_individual(org.repr());
                            if(exists<EPIGENETIC_INFO>(org)) {
                                put<EPIGENETIC_INFO>(get<EPIGENETIC_INFO>(org),*o);
                            }
                            
                            p->append(o);
                        }
                    } else {
                        // and fill up the offspring population with copies of the germ:
                        typename EA::individual_type::individual_ptr_type o=p->make_individual(org.repr());
                        if(exists<EPIGENETIC_INFO>(org)) {
                            put<EPIGENETIC_INFO>(get<EPIGENETIC_INFO>(org),*o);
                        }
                        p->append(o);
                        
                    }
                    
                    ++p_size;
                    if (p_size >= num_parents) break;
                }
                
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

template <typename EA>
struct propagule_size_tracking : end_of_update_event<EA> {
    propagule_size_tracking(EA& ea) : end_of_update_event<EA>(ea), _df("prop_size.dat") {
        _df.add_field("update")
        .add_field("mean_prop_size");
        
    }
    
    //! Destructor.
    virtual ~propagule_size_tracking() {
    }
    
    //! Track how many task-switches are being performed!
    virtual void operator()(EA& ea) {
        if ((ea.current_update() % 100) == 0) {
            double mean_ps = 0.0;
            double sub_pop_size = 0.0;
            
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                ++sub_pop_size;
                mean_ps += get<ACTUAL_PROP_SIZE>(*i,1.0);
            }
            mean_ps /= sub_pop_size;
            _df.write(ea.current_update())
            .write(mean_ps)
            .endl();
            
        }
        
    }
    datafile _df;
    
};



#endif
