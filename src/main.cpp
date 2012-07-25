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
#include <ea/artificial_life/artificial_life.h>
#include <ea/artificial_life/hardware.h>
#include <ea/artificial_life/isa.h>
#include <ea/artificial_life/spatial.h>
#include <ea/artificial_life/datafiles/reactions.h>
#include <ea/artificial_life/datafiles/generation_priority.h>
#include <ea/cmdline_interface.h>

using namespace ea;

/*! Artificial life simulation definition.
 */
typedef artificial_life<
hardware, isa, spatial
> al_type;


/*! 
 */
template <typename EA>
class ealife : public cmdline_interface<EA> {
public:
    virtual void configure(EA& ea) {
        add_task<tasks::task_not,resources::unlimited,catalysts::additive<2> >("not", ea); // 1
        add_task<tasks::task_nand,resources::unlimited,catalysts::additive<2> >("nand", ea); //1
        add_task<tasks::task_and,resources::unlimited,catalysts::additive<4> >("and", ea); // 2
        add_task<tasks::task_ornot,resources::unlimited,catalysts::additive<4> >("ornot", ea); // 2
        add_task<tasks::task_or,resources::unlimited,catalysts::additive<8> >("or", ea); // 3
        add_task<tasks::task_andnot,resources::unlimited,catalysts::additive<8> >("andnot", ea); // 3
        add_task<tasks::task_nor,resources::unlimited,catalysts::additive<16> >("nor", ea); // 4
        add_task<tasks::task_xor,resources::unlimited,catalysts::additive<16> >("xor", ea); // 4
        add_task<tasks::task_equals,resources::unlimited,catalysts::additive<32> >("equals", ea); // 5
    }
    
    virtual void gather_options() {
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
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::record_reactions_event>(this,ea);
        add_event<datafiles::generation_priority>(this,ea);
    };
};
LIBEA_CMDLINE_INSTANCE(al_type, ealife);
