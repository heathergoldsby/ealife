#include "age_poly.h"
#include "subpopulation_founder.h"

#include <ea/line_of_descent.h>


//! Configuration object for an EA.
template <typename EA>
struct ts_configuration : public abstract_configuration<EA> {
    
    typedef typename EA::tasklib_type::task_ptr_type task_ptr_type;
    typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
    
    
    //! Called as the final step of EA construction.
    void construct(EA& ea) {
        using namespace ea::instructions;
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
        append_isa<get_xy>(ea);
        
        // Add tasks
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", ea);

        
        put<TASK_LETHALITY_PROB>(0, *task_not);
        put<TASK_LETHALITY_PROB>(0.25, *task_nand);
     
        
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
ts_configuration, spatial, empty_neighbor, round_robin
> ea_type;

template <typename EA>
struct mp_configuration : public abstract_configuration<EA> {
};


//! Meta-population definition.
typedef meta_population<
subpopulation_founder <ea_type> 
, mp_configuration> mea_type;


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
        
        // ts specific options
        add_option<TASK_LETHALITY_PROB>(this);
        add_option<LAST_TASK>(this);
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        add_option<EACH_TASK_THRESH>(this);
        
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<ape_two_task_replication>(this,ea);
        add_event<task_performed_tracking>(this,ea);
        add_event<task_lethality>(this,ea);
        add_event<founder_event>(this,ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
