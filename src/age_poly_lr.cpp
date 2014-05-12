#include "age_poly.h"
#include <ea/digital_evolution/population_founder.h>
#include "multi_birth_selfrep_not_ancestor.h"



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
        task_ptr_type task_and = make_task<tasks::task_and,catalysts::additive<0> >("and", ea);
        task_ptr_type task_ornot = make_task<tasks::task_ornot,catalysts::additive<0> >("ornot", ea);
        task_ptr_type task_or = make_task<tasks::task_or,catalysts::additive<0> >("or", ea);
        task_ptr_type task_andnot = make_task<tasks::task_andnot,catalysts::additive<0> >("andnot", ea);
        task_ptr_type task_nor = make_task<tasks::task_nor,catalysts::additive<0> >("nor", ea);
        task_ptr_type task_xor = make_task<tasks::task_xor,catalysts::additive<0> >("xor", ea);
        task_ptr_type task_equals = make_task<tasks::task_equals,catalysts::additive<0> >("equals", ea);

        double init_amt = get<RES_INITIAL_AMOUNT>(ea, 0);
        double inflow = get<RES_INFLOW_AMOUNT>(ea,0);
        double outflow = get<RES_OUTFLOW_FRACTION>(ea,0);
        double frac = get<RES_FRACTION_CONSUMED>(ea,0);
        
        // initial amount (unit), inflow (unit), outflow (percentage), percent consumed, ea
        resource_ptr_type resA = make_resource("resA", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resB = make_resource("resB", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resC = make_resource("resC", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resD = make_resource("resD", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resE = make_resource("resE", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resF = make_resource("resF", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resG = make_resource("resG", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resH = make_resource("resH", init_amt, inflow, outflow, frac, ea);
        resource_ptr_type resI = make_resource("resI", init_amt, inflow, outflow, frac, ea);

        
        
        task_not->consumes(resA);
        task_nand->consumes(resB);
        task_and->consumes(resC);
        task_ornot->consumes(resD);
        task_or->consumes(resE);
        task_andnot->consumes(resF);
        task_nor->consumes(resG);
        task_xor->consumes(resH);
        task_equals->consumes(resI);
        
        put<TASK_LETHALITY_PROB>(get<NOT_LETHALITY_PROB>(ea), *task_not);
        put<TASK_LETHALITY_PROB>(get<NAND_LETHALITY_PROB>(ea), *task_nand);
        put<TASK_LETHALITY_PROB>(get<AND_LETHALITY_PROB>(ea), *task_and);
        put<TASK_LETHALITY_PROB>(get<ORNOT_LETHALITY_PROB>(ea), *task_ornot);
        put<TASK_LETHALITY_PROB>(get<OR_LETHALITY_PROB>(ea), *task_or);
        put<TASK_LETHALITY_PROB>(get<ANDNOT_LETHALITY_PROB>(ea), *task_andnot);
        put<TASK_LETHALITY_PROB>(get<NOR_LETHALITY_PROB>(ea), *task_nor);
        put<TASK_LETHALITY_PROB>(get<XOR_LETHALITY_PROB>(ea), *task_xor);
        put<TASK_LETHALITY_PROB>(get<EQUALS_LETHALITY_PROB>(ea), *task_equals);
        
    }
    
    
    //! Called to generate the initial EA population.
    void initial_population(EA& ea) {
        generate_ancestors(multibirth_selfrep_not_ancestor(), 1, ea);
        
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
        
        // ape lr specific options
        add_option<TASK_LETHALITY_PROB>(this);
        add_option<LAST_TASK>(this);
        add_option<GERM_MUTATION_PER_SITE_P>(this);
        add_option<GROUP_REP_THRESHOLD>(this);
        
        add_option<NOT_LETHALITY_PROB>(this);
        add_option<NAND_LETHALITY_PROB>(this);
        add_option<AND_LETHALITY_PROB>(this);
        add_option<ORNOT_LETHALITY_PROB>(this);
        add_option<OR_LETHALITY_PROB>(this);
        add_option<ANDNOT_LETHALITY_PROB>(this);
        add_option<NOR_LETHALITY_PROB>(this);
        add_option<XOR_LETHALITY_PROB>(this);
        add_option<EQUALS_LETHALITY_PROB>(this);
        
        // initial amount (unit), inflow (unit), outflow (percentage), percent consumed
        add_option<RES_INITIAL_AMOUNT>(this);
        add_option<RES_INFLOW_AMOUNT>(this);
        add_option<RES_OUTFLOW_FRACTION>(this);
        add_option<RES_FRACTION_CONSUMED>(this);
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<ape_lr_replication>(this,ea);
        add_event<task_performed_tracking>(this,ea);
        add_event<population_founder_event>(this,ea);
        add_event<task_first_age_tracking>(this,ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);

