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

#include "subpopulation_lod_analysis.h"
#include "lod_knockouts.h"
#include "multi_birth_selfrep_not_ancestor.h"

#include <ea/digital_evolution/population_founder.h>
#include <ea/line_of_descent.h>


//! Configuration object for an EA.
struct gls_configuration : public default_configuration {
    
//    typedef typename EA::task_library_type::task_ptr_type task_ptr_type;
//    typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
//    
    
    template <typename EA>
    void after_construction(EA& ea) {
        using namespace ealib::instructions;
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
        append_isa<fixed_input>(ea);
        append_isa<output>(ea);//25
        append_isa<become_soma>(ea);
        append_isa<if_germ>(ea);
        append_isa<if_soma>(ea);
        append_isa<donate_res_to_group>(ea);
        append_isa<get_xy>(ea);
        append_isa<apoptosis>(ea);
        
        add_event<task_mutagenesis>(ea);
        add_event<gs_apoptosis_event>(ea);
//        add_event<gs_inherit_event>(ea);
        add_event<task_resource_consumption>(ea);
        
    }

    template <typename EA>
    void initialize(EA& ea) {
        typedef typename EA::task_library_type::task_ptr_type task_ptr_type;
        typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
        
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
        
        put<TASK_MUTATION_MULT>(get<NOT_MUTATION_MULT>(ea), *task_not);
        put<TASK_MUTATION_MULT>(get<NAND_MUTATION_MULT>(ea), *task_nand);
        put<TASK_MUTATION_MULT>(get<AND_MUTATION_MULT>(ea), *task_and);
        put<TASK_MUTATION_MULT>(get<ORNOT_MUTATION_MULT>(ea), *task_ornot);
        put<TASK_MUTATION_MULT>(get<OR_MUTATION_MULT>(ea), *task_or);
        put<TASK_MUTATION_MULT>(get<ANDNOT_MUTATION_MULT>(ea), *task_andnot);
        put<TASK_MUTATION_MULT>(get<NOR_MUTATION_MULT>(ea), *task_nor);
        put<TASK_MUTATION_MULT>(get<XOR_MUTATION_MULT>(ea), *task_xor);
        put<TASK_MUTATION_MULT>(get<EQUALS_MUTATION_MULT>(ea), *task_equals);
        
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
    }
    
    //! Called to generate the initial EA population.
    template <typename EA>
    void initial_population(EA& ea) {
        generate_ancestors(selfrep_not_ancestor(), 1, ea);
    }
};


/*! Artificial life simulation definition.
 */
/*
typedef digital_evolution<
gls_configuration, spatial, empty_neighbor, round_robin
> ea_type;
*/
typedef digital_evolution
< gls_configuration
, organism< >
, multibirth_selfrep_not_ancestor
, recombination::asexual
, round_robin
, empty_neighbor
> ea_type;

/*
 
template <typename EA>
struct mp_configuration : public abstract_configuration<EA> {
    void initial_population(EA& ea) {
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            (*i).founder() = (**(*i).population().begin());
        }
    }
};
*/

/*
//! Meta-population definition.
typedef meta_population<
population_lod<population_founder<ea_type> >, 
 mp_configuration> mea_type;
 */

typedef metapopulation
< subpopulation<ea_type>
> mea_type;



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

        add_option<ANALYSIS_INPUT>(this);
        
        // gls specific options
        add_option<TASK_MUTATION_PER_SITE_P>(this);
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        add_option<GROUP_REP_THRESHOLD>(this);
        
        // gls mutation multipliers
        add_option<NOT_MUTATION_MULT>(this);
        add_option<NAND_MUTATION_MULT>(this);
        add_option<AND_MUTATION_MULT>(this);
        add_option<ORNOT_MUTATION_MULT>(this);
        add_option<OR_MUTATION_MULT>(this);
        add_option<ANDNOT_MUTATION_MULT>(this);
        add_option<NOR_MUTATION_MULT>(this);
        add_option<XOR_MUTATION_MULT>(this);
        add_option<EQUALS_MUTATION_MULT>(this);
        
    }
    
    virtual void gather_tools() {

        //add_tool<ealib::analysis::lod_knockouts>(this);
        //add_tool<ealib::analysis::lod_gls_circle_square_plot>(this);
        //add_tool<ealib::analysis::lod_gls_germ_soma_mean_var>(this);
        //add_tool<ealib::analysis::lod_gls_aging_res_over_time>(this);
        //add_tool<ealib::analysis::lod_gls_aging_res_over_time_compact>(this);
        //add_tool<ealib::analysis::lod_gls_task_count>(this);

        
    }
    
    virtual void gather_events(EA& ea) {
        add_event<gls_replication>(ea);
        add_event<task_performed_tracking>(ea);
        add_event<apoptosis_tracking>(ea);
//        add_event<datafiles::mrca_lineage>(ea);
//        add_event<population_founder_event>(ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
