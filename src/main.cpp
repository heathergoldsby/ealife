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
#include <ea/artificial_life/hardware.h>
#include <ea/artificial_life/isa.h>
#include <ea/artificial_life/spatial.h>
#include <ea/artificial_life/datafiles/reactions.h>
#include <ea/artificial_life/datafiles/generation_priority.h>
#include <ea/cmdline_interface.h>
using namespace ea;


/*! Group-tracking spatial environment.
 */
template <typename EA>
struct grouping : spatial<EA> {
    typedef spatial<EA> base_type;
    typedef typename EA::individual_ptr_type individual_ptr_type;
    
    //! Properties for vertices in the group graph.
    struct grouping_vertex_properties {
        individual_ptr_type p;
    };
    
    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, grouping_vertex_properties> graph_type;
    typedef std::map<long,typename graph_type::vertex_descriptor> nv_map_type;
    
    //! Constructor.
    grouping() {
    }

    //! Destructor.
    virtual ~grouping() {
    }
    
    //! Initialize this environment.
    void initialize(EA& ea) {
        base_type::initialize(ea);
        _inheritance = ea.events().inheritance.connect(boost::bind(&grouping::inheritance, this, _1, _2, _3));
        _death =  ea.events().death.connect(boost::bind(&grouping::death, this, _1, _2));
    }
    
    //! Called when an offspring inherits from parents.
    void inheritance(typename EA::population_type& parents,
                     typename EA::individual_type& offspring,
                     EA& ea) {
    }

    //! Called when an individual dies.
    void death(typename EA::individual_type& individual, EA& ea) {
    }

    graph_type _g; //!< Group graph.
    nv_map_type _map; //!< Map of individual's name to vertex in the group graph.
    boost::signals::scoped_connection _inheritance; //<! Connection to inheritance event.
    boost::signals::scoped_connection _death; //!< Connection to death event.
};
    

/*! Artificial life simulation definition.
 */
typedef artificial_life<
hardware, isa, spatial
> al_type;


/*! 
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    typedef typename EA::tasklib_type::task_ptr_type task_ptr_type;
    typedef typename EA::environment_type::resource_ptr_type resource_ptr_type;
    
    virtual void configure(EA& ea) {
        task_ptr_type task_not = make_task<tasks::task_not,catalysts::additive<2> >("not", ea);
        task_ptr_type task_nand = make_task<tasks::task_nand,catalysts::additive<2> >("nand", ea);

        resource_ptr_type resA = make_resource("resA", ea);
        resource_ptr_type resB = make_resource("resB", 100.0, 1.0, 0.01, 0.05, ea);
        
        task_nand->consumes(resB);
    }
    
    virtual void gather_options() {
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
        add_option<SPATIAL_X>(this);
        add_option<SPATIAL_Y>(this);
    }
    
    virtual void gather_tools() {
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::record_reactions_event>(this,ea);
        add_event<datafiles::generation_priority>(this,ea);
    };
};
LIBEA_CMDLINE_INSTANCE(al_type, cli);
