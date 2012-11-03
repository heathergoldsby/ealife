//
//  group_lod_analysis.h
//  ealife
//
//  Created by Heather Goldsby on 10/18/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_SUBPOPULATION_LOD_ANALYSIS_H_
#define _EALIFE_SUBPOPULATION_LOD_ANALYSIS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/interface.h>
#include <ea/line_of_descent.h>
#include <ea/analysis/tool.h>
#include <ea/digital_evolution/isa.h>
#include <ea/digital_evolution/spatial.h>



namespace ea {
    namespace analysis {
        
        /* Rerun a LoD, where a LoD is the sequential list of eas that contributed to a final dominant
         ea. Each ea is tagged with a founder, which is used for the rerun.*/
        
        /*! lod_gls_aging is designed to work for the gls project. It reruns each subpopulation along a
         line of descent and records how the subpopulation 'ages' as measured by how rapidly it
         accumulates enough resources to produce a new group.
         */
        template <typename EA>
        struct lod_gls_aging : public ea::analysis::unary_function<EA> {
            static const char* name() { return "lod_gls_aging"; }
            
            virtual void operator()(EA& ea) {
                using namespace ea;
                using namespace ea::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                datafile df("gls_aging.dat");
                df.add_field("lod depth [depth]");
                
                int lod_depth = 0;
                // skip def ancestor (that's what the +1 does)
                for( ; i!=lod.end(); ++i) {
                    
                    df.write(lod_depth);
                    
                    // **i is the EA, AS OF THE TIME THAT IT DIED!
                    
                    // To replay, need to create a new ea
                    // setup the population (really, an ea):
                    typename EA::individual_ptr_type p = ea.make_individual();
                    p->rng().reset(get<RNG_SEED>(**i));
                    
                    // setup the founder
                    typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
                    o->hw().initialize();
                    p->append(o);
                    //                    p->env().reset_resources();
                    
                    // replay!
                    int res_reset_inc = 500;
                    int res_resource_thresh = 500;
                    int cur_update = 0;
                    // and run till the group amasses the right amount of resources
                    while ((get<GROUP_RESOURCE_UNITS>(*p,0) < 5000) &&
                           (cur_update < 10000)){
                        double blah = get<GROUP_RESOURCE_UNITS>(*p);
                        p->update();
                        ++cur_update;
                        if (get<GROUP_RESOURCE_UNITS>(*p,0) > res_resource_thresh) {
                            df.write(cur_update);
                            p->env().reset_resources();
                            
                            res_resource_thresh += res_reset_inc;
                        }
                    }
                    
                    
                    // catch missing data
                    int cur_res = get<GROUP_RESOURCE_UNITS>(*p);
                    while (cur_res < 5000) {
                        df.write("99999");
                        cur_res += res_reset_inc;
                    }
                    
                    df.endl();
                    
                    ++lod_depth;
                }
            }
            
        };
        
        
        /*! lod_knockouts reruns each subpopulation along a line of descent and records how the subpopulation
         fares with key coordination instructions removed.
         */
        template <typename EA>
        struct lod_knockouts : public ea::analysis::unary_function<EA> {
            static const char* name() { return "lod_knockouts"; }
            
            virtual void operator()(EA& ea) {
                using namespace ea;
                using namespace ea::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                datafile df("lod_knockouts.dat");
                df.add_field("lod depth [depth]")
                .add_field("no knockouts")
                .add_field("rx knockedout")
                .add_field("location knockedout");
                
                
                int lod_depth = 0;
                // skip def ancestor (that's what the +1 does)
                for( ; i!=lod.end(); ++i) {
                    
                    df.write(lod_depth);
                    
                    // **i is the EA, AS OF THE TIME THAT IT DIED!
                    
                    // To replay, need to create new eas for each knockout exper.
                    // setup the population (really, an ea):
                    typename EA::individual_ptr_type control_ea = ea.make_individual();
                    control_ea->rng().reset(get<RNG_SEED>(**i));
                    
                    typename EA::individual_ptr_type knockout_rx_ea = ea.make_individual();
                    knockout_rx_ea->rng().reset(get<RNG_SEED>(**i));
                    knockout_isa<instructions::nop_x>("rx_msg", *knockout_rx_ea);
                    
                    typename EA::individual_ptr_type knockout_location_ea = ea.make_individual();
                    knockout_location_ea->rng().reset(get<RNG_SEED>(**i));
                    knockout_isa<instructions::nop_x>("get_xy", *knockout_location_ea);
                    
                    // setup the founder
                    typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
                    o->hw().initialize();
                    control_ea->append(o);
                    
                    // setup send knockout founder
                    typename EA::individual_type::individual_ptr_type ko_s = (*i)->make_individual((*i)->founder().repr());
                    ko_s->hw().initialize();
                    knockout_rx_ea->append(ko_s);
                    
                    // setup location knockout founder
                    typename EA::individual_type::individual_ptr_type ko_l = (*i)->make_individual((*i)->founder().repr());
                    ko_l->hw().initialize();
                    knockout_location_ea->append(ko_l);
                    
                    
                    // replay! till the group amasses the right amount of resources
                    // or exceeds its window...
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                           (cur_update < update_max)){
                        control_ea->update();
                        ++cur_update;
                    }
                    df.write(cur_update);
                    
                    cur_update = 0;
                    while ((get<GROUP_RESOURCE_UNITS>(*knockout_rx_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_rx_ea)) &&
                           (cur_update < update_max)){
                        knockout_rx_ea->update();
                        ++cur_update;
                    }
                    df.write(cur_update);
                    
                    cur_update = 0;
                    while ((get<GROUP_RESOURCE_UNITS>(*knockout_location_ea,0) < get<GROUP_REP_THRESHOLD>(*knockout_location_ea)) &&
                           (cur_update < update_max)){
                        knockout_location_ea->update();
                        ++cur_update;
                    }
                    df.write(cur_update);
                        
                    
                    df.endl();
                    
                    ++lod_depth;
                }
            }
            
        };
        
        
        /*! lod_knockouts reruns each subpopulation along a line of descent - setup for circle / square plot
         */
        template <typename EA>
        struct lod_gls_rerun : public ea::analysis::unary_function<EA> {
            static const char* name() { return "lod_gls_rerun"; }
            
            virtual void operator()(EA& ea) {
                using namespace ea;
                using namespace ea::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                datafile df("lod_gls_rerun.dat");
                df.add_field("lod depth [depth]");
                
                
                int lod_depth = 0;
                // skip def ancestor (that's what the +1 does)
                for( ; i!=lod.end(); ++i) {
                    
                    df.write(lod_depth);
                    
                    // **i is the EA, AS OF THE TIME THAT IT DIED!
                    
                    // To replay, need to create new eas for each knockout exper.
                    // setup the population (really, an ea):
                    typename EA::individual_ptr_type control_ea = ea.make_individual();
                    control_ea->rng().reset(get<RNG_SEED>(**i));
                    
                    // setup the founder
                    typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
                    o->hw().initialize();
                    control_ea->append(o);
                    
                    // replay! till the group amasses the right amount of resources
                    // or exceeds its window...
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<GROUP_RESOURCE_UNITS>(*control_ea,0) < get<GROUP_REP_THRESHOLD>(*control_ea)) &&
                           (cur_update < update_max)){
                        control_ea->update();
                        ++cur_update;
                    }
                    df.write(cur_update);
                    
                    // grab info based on location...
                    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                            typename EA::individual_type::environment_type::location_type l = control_ea->env().location(x,y);
                            if (l.occupied()) {
                                bool gs_status = get<GERM_STATUS>(*l.inhabitant());
                                df.write(get<GERM_STATUS>(*l.inhabitant(), 0))
                                .write(get<WORKLOAD>(*l.inhabitant(),0));
                            } else {
                                df.write("2")
                                .write("0");
                            }
                            
                        }
                    }
                    
                    df.endl();
                    
                    ++lod_depth;
                }
            }
            
        };

    }
    
}

#endif
