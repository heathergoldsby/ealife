#include <ea/artificial_life/artificial_life.h>
#include <ea/artificial_life/hardware.h>
#include <ea/artificial_life/isa.h>
#include <ea/artificial_life/topology.h>
#include <ea/artificial_life/task_library.h>
#include <ea/mutation.h>
#include <ea/cmdline_interface.h>
using namespace ea;


/*! Artificial life simulation definition.
 */
typedef artificial_life<
hardware, isa
> al_type;


/*! 
 */
template <typename EA>
class ealife : public cmdline_interface<EA> {
public:
    virtual void configure(EA& ea) {
        add_task<tasks::task_nand,resources::unlimited,catalysts::power>("nand", ea);
    }
    
    virtual void gather_options() {
        add_option<POPULATION_SIZE>(this);
        add_option<REPLACEMENT_RATE_P>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_UNIFORM_INT_MAX>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<MUTATION_DUPLICATION_P>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);        
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
    };
};
LIBEA_CMDLINE_INSTANCE(al_type, ealife);
