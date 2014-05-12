#include <ea/digital_evolution/population_founder.h>
#include <ea/line_of_descent.h>
#include <ea/digital_evolution.h>
#include <ea/cmdline_interface.h>
using namespace ealib;


#include "ts.h"
#include "stripes.h"
#include "resource_phenotype.h"
#include "propagule_size.h"
#include "lod_knockouts.h"
#include "multi_birth_selfrep_not_ancestor.h"
#include "multi_birth_selfrep_not_nand_ancestor.h"


//! Configuration object for an EA.
struct configuration : public default_configuration {
    
    //! Called as the final step of EA construction (must not depend on configuration parameters)
    template <typename EA>
    void after_construction(EA& ea) {
        using namespace instructions;
        append_isa<nop_a>(0,ea);
        append_isa<nop_b>(0,ea);
        append_isa<nop_c>(0,ea);
        append_isa<nop_x>(ea);
        append_isa<mov_head>(ea);
        append_isa<if_label>(ea);
        append_isa<h_search>(ea);
        append_isa<nand>(ea);
        append_isa<push>(ea);
        append_isa<pop>(ea);
        append_isa<swap>(ea);
        append_isa<inc>(ea);
        append_isa<dec>(ea);
        append_isa<tx_msg_check_task>(ea);
        append_isa<tx_msg>(ea);
        append_isa<rx_msg>(ea);
        append_isa<bc_msg>(ea);
        append_isa<rotate>(ea);
        append_isa<rotate_cw>(ea);
        append_isa<rotate_ccw>(ea);
        append_isa<if_less>(ea);
        append_isa<h_alloc>(ea);
        append_isa<h_copy>(ea);
        append_isa<h_divide_soft_parent_reset>(ea);
        append_isa<fixed_input>(ea);
        append_isa<output>(ea);
        append_isa<donate_res_to_group>(ea);
        //        append_isa<get_xy>(ea);
        append_isa<if_equal>(ea);
        append_isa<if_not_equal>(ea);
        append_isa<jump_head>(ea);
        append_isa<is_neighbor>(ea);
        append_isa<is_origin>(ea);
        append_isa<get_epigenetic_info>(ea);
        append_isa<set_epigenetic_info>(ea);
        
        // SOMA
        //        append_isa<inc_propagule_size>(ea);
        //        append_isa<dec_propagule_size>(ea);
        //        append_isa<get_propagule_size>(ea);
        //
        //        append_isa<become_soma>(ea);
        //        append_isa<if_soma>(ea);
        //        append_isa<if_germ>(ea);
        
        add_event<task_resource_consumption>(ea);
        add_event<task_switching_cost>(ea);
        add_event<ts_birth_event>(ea);
    }
    
    //! Initialize! Things are live and are mostly setup. All the objects are there, but they
    // may not have the parameters that they need.
    template <typename EA>
    void initialize(EA& ea) {
        typedef typename EA::task_library_type::task_ptr_type task_ptr_type;
        typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
        
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", ea);
        
        resource_ptr_type resA = make_resource("resA", ea);
        resource_ptr_type resB = make_resource("resB", ea);
        
        task_not->consumes(resA);
        task_nand->consumes(resB);
        
        
    }
    
    //! Called to generate the initial EA population.
    template <typename EA>
    void initial_population(EA& ea) {
        
        int ancest = get<ANCESTOR>(ea, 0);
        switch (ancest) {
            case 0:
                generate_ancestors(multibirth_selfrep_not_ancestor(), 1, ea);
                break;
            case 1:
                generate_ancestors(multibirth_selfrep_not_nand_ancestor(), 1, ea);
                break;
        }
        
        
    }
    
};


/*! Artificial life simulation definition.
 */
typedef digital_evolution
< configuration
, organism< >
, multibirth_selfrep_not_ancestor
, recombination::asexual
, round_robin
, empty_facing_neighbor
> ea_type;


////template <typename EA>
//struct mp_configuration : public default_configuration {
//    template <typename EA>
//    void initial_population(EA& ea) {
//        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
//            (*i).founder() = (**(*i).population().begin());
//        }
//    }
//};


//
////! Meta-population definition.
//typedef meta_population<
//population_lod<population_founder<ea_type> >,
//mp_configuration> mea_type;
/*
< typename Subpopulation
, typename AncestorGenerator=ancestors::default_representation
, typename MutationOperator=mutation::operators::no_mutation
, typename RecombinationOperator=recombination::no_recombination
, typename GenerationalModel=generational_models::isolated_subpopulations
, typename EarlyStopCondition=dont_stop
, typename UserDefinedConfiguration=default_configuration
, typename PopulationGenerator=fill_metapopulation
> class metapopulation {
    */

// hjg - ask about lod...
//! Metapopulation definition:
typedef metapopulation
< subpopulation<ea_type>,
ancestors::default_representation,
mutation::operators::subpopulation_mutator<ea_type::mutation_operator_type>
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
        add_option<NUM_PROPAGULE_GERM>(this);
        add_option<NUM_PROPAGULE_CELL>(this);
        
        // ts specific options
        add_option<TASK_SWITCHING_COST>(this);
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        
        // stripes
        add_option<ANCESTOR>(this);
        add_option<METAPOP_COMPETITION_PERIOD>(this);
        add_option<TOURNAMENT_SELECTION_N>(this);
        add_option<TOURNAMENT_SELECTION_K>(this);
        
        
    }
    
    virtual void gather_tools() {
        
    }
    
    virtual void gather_events(EA& ea) {
        //        add_event<ts_replication_propagule>(this,ea);
        //        add_event<ps_size_propagule2>(this,ea);
        //        add_event<datafiles::fitness_dat>(ea);
        add_event<permute_stripes>(ea);
        add_event<task_performed_tracking>(ea);
        add_event<task_switch_tracking>(ea);
        //add_event<propagule_size_tracking>(ea);
        //        add_event<population_founder_event>(this,ea);
        //add_event<reward_tracking>(ea);
    }
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);