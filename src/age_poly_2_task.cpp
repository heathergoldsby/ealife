#include "age_poly.h"
#include <ea/digital_evolution/population_founder.h>
#include "multi_birth_selfrep_not_nand_ancestor.h"
#include "multi_birth_selfrep_nand_not_ancestor.h"

#include "selfrep_not_ancestor.h"


#include <ea/line_of_descent.h>


//! Configuration object for an EA.
template <typename EA>
struct ts_configuration : public abstract_configuration<EA> {
    
    typedef typename EA::task_library_type::task_ptr_type task_ptr_type;
    typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
    
    
    //! All code. No data. Never use 
    // parameters, depend on input, read a file, etc. It is done *BEFORE* initialize. 
    // If you depend on metadata, you MUST go into intialize. 
    void configure(EA& ea) {
        using namespace ealib::instructions;
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
        append_isa<h_divide>(ea);
        append_isa<h_divide_soft_parent_reset>(ea);
        append_isa<fixed_input>(ea);
        append_isa<output>(ea);
        append_isa<donate_res_to_group>(ea);
        append_isa<get_xy>(ea);
        append_isa<if_equal>(ea);
        append_isa<if_not_equal>(ea);
        append_isa<jump_head>(ea);
        append_isa<get_age>(ea);
        
        
        add_event<task_resource_consumption>(this,ea);
        add_event<task_lethality>(this,ea);
        add_event<task_first_age>(this,ea);

        
    }

    //! Initialize! Things are live and are mostly setup. All the objects are there, but they
    // may not have the parameters that they need. 
    void initialize(EA& ea) {
        // Add tasks
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<0> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<0> >("nand", ea);
        
        resource_ptr_type resA = make_resource("resA", ea);
        resource_ptr_type resB = make_resource("resB", ea);
        
        task_not->consumes(resA);
        task_nand->consumes(resB);

        put<TASK_LETHALITY_PROB>(get<NOT_LETHALITY_PROB>(ea), *task_not);
        put<TASK_LETHALITY_PROB>(get<NAND_LETHALITY_PROB>(ea), *task_nand);
    }

    
    //! Called to generate the initial EA population.
    void initial_population(EA& ea) {
        int ancest = get<ANCESTOR>(ea, 0);
        if (ancest == 1) {
            generate_ancestors(multibirth_selfrep_nand_not_ancestor(), 1, ea);
        } else {
            generate_ancestors(multibirth_selfrep_not_nand_ancestor(), 1, ea);
        }

    }
};


/*! Artificial life simulation definition.
 */
typedef digital_evolution<
ts_configuration, spatial, empty_neighbor, round_robin
> ea_type;

template <typename EA>
struct mp_configuration : public abstract_configuration<EA> {
    void initial_population(EA& ea) {
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            (*i).founder() = (**(*i).population().begin());
        }
    }

};


//! Meta-population definition.
typedef meta_population<
population_lod<population_founder<ea_type> >,
mp_configuration> mea_type;

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
        
        // age poly specific options
        add_option<TASK_LETHALITY_PROB>(this);
        add_option<LAST_TASK>(this);
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        add_option<EACH_TASK_THRESH>(this);
        add_option<ANCESTOR>(this);
        
        add_option<NOT_LETHALITY_PROB>(this);
        add_option<NAND_LETHALITY_PROB>(this);
        
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<ape_two_task_replication>(this,ea);
        add_event<task_performed_tracking>(this,ea);
        add_event<population_founder_event>(this,ea);
        add_event<task_first_age_tracking>(this,ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
