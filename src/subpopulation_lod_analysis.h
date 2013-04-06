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
#include <ea/line_of_descent.h>
#include <ea/analysis/tool.h>
#include <ea/digital_evolution/isa.h>
#include <ea/digital_evolution/spatial.h>



namespace ealib {
    namespace analysis {
        
        /* Rerun a LoD, where a LoD is the sequential list of eas that contributed to a final dominant
         ea. Each ea is tagged with a founder, which is used for the rerun.*/
        


        
        /*! lod_gls_aging_res_over_time */
        template <typename EA>
        struct lod_gls_aging_res_over_time_compact : public ealib::analysis::unary_function<EA> {
            static const char* name() { return "lod_gls_aging_res_over_time_compact"; }
            
            virtual void operator()(EA& ea) {
                using namespace ealib;
                using namespace ealib::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                
                datafile df("lod_gls_aging_res_over_time_compact.dat");
                df.add_field("lod_depth")
                .add_field("res1000")
                .add_field("res2000")
                .add_field("res3000")
                .add_field("res4000")
                .add_field("res5000")
                .add_field("res6000")
                .add_field("res7000")
                .add_field("res8000")
                .add_field("res9000")
                .add_field("res10000");
                
                
                int lod_depth = 0;
                // skip def ancestor (that's what the +1 does)
                for( ; i!=lod.end(); ++i) {
                    
                    
                    if ((lod_depth % 10) != 0) {
                        lod_depth++;
                        continue;
                    }
                    
                    
                    
                    // **i is the EA, AS OF THE TIME THAT IT DIED!
                    
                    // To replay, need to create a new ea
                    // setup the population (really, an ea):
                    typename EA::individual_ptr_type p = ea.make_individual();
                    p->rng().reset(get<RNG_SEED>(**i));
                    
                    // setup the founder
                    typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
                    o->hw().initialize();
                    p->append(o);
                    
                    // replay!
                    int res_reset_inc = 500;
                    int res_resource_thresh = 500;
                    int cur_update = 0;
                    int max_update = 10000;
                    int num_rep = 0;
                    int prev_res = 0;
                    // and run till the group amasses the right amount of resources
                    df.write(lod_depth);
                    
                    while (cur_update < max_update){
                        p->update();
                        ++cur_update;
                        
                        if (get<GROUP_RESOURCE_UNITS>(*p,0) >= res_resource_thresh) {
                            p->env().reset_resources();
                            
                            res_resource_thresh += res_reset_inc;
                            num_rep++;
                            
                        }
                        if ((cur_update % 1000) == 0) {
                            double cur_res = get<GROUP_RESOURCE_UNITS>(*p,0);
                            double res = cur_res - prev_res;
                            prev_res = cur_res;
                            df.write(res);
                        }
                    }
                    
                    df.endl();
                    ++lod_depth;
                }
            }
            
        };


        /*! lod_gls_aging_res_over_time */
        template <typename EA>
        struct lod_gls_aging_res_over_time : public ealib::analysis::unary_function<EA> {
            static const char* name() { return "lod_gls_aging_res_over_time"; }
            
            virtual void operator()(EA& ea) {
                using namespace ealib;
                using namespace ealib::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                
                datafile df("lod_gls_aging_res_over_time.dat");
                df.add_field("lod_depth")
                .add_field("update")
                .add_field("res");

                
                int lod_depth = 0;
                // skip def ancestor (that's what the +1 does)
                for( ; i!=lod.end(); ++i) {
                    
                    
                    if ((lod_depth % 10) != 0) {
                        lod_depth++;
                        continue;
                    }
                    
                    
                    
                    // **i is the EA, AS OF THE TIME THAT IT DIED!
                    
                    // To replay, need to create a new ea
                    // setup the population (really, an ea):
                    typename EA::individual_ptr_type p = ea.make_individual();
                    p->rng().reset(get<RNG_SEED>(**i));
                    
                    // setup the founder
                    typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
                    o->hw().initialize();
                    p->append(o);
                    
                    // replay!
                    int res_reset_inc = 500;
                    int res_resource_thresh = 500;
                    int cur_update = 0;
                    int max_update = 10000;
                    int num_rep = 0;
                    int prev_res = 0;
                    // and run till the group amasses the right amount of resources
                    while (cur_update < max_update){
                        p->update();
                        ++cur_update;
                        
                        if (get<GROUP_RESOURCE_UNITS>(*p,0) >= res_resource_thresh) {
                            p->env().reset_resources();
                            
                            res_resource_thresh += res_reset_inc;
                            num_rep++;
                            
                        }
                        if ((cur_update % 1000) == 0) {
                            double cur_res = get<GROUP_RESOURCE_UNITS>(*p,0);
                            double res = cur_res - prev_res;
                            prev_res = cur_res;
                            df.write(lod_depth);
                            df.write(cur_update);
                            df.write(res);
                            df.endl();
                        }
                    }
                    
//                    df.write(num_rep);
                    
                    ++lod_depth;
                }
            }
            
        };
        
        
        /*! lod_knockouts reruns each subpopulation along a line of descent and records how the subpopulation
         fares with key coordination instructions removed.
         */
        template <typename EA>
        struct lod_knockouts : public ealib::analysis::unary_function<EA> {
            static const char* name() { return "lod_knockouts"; }
            
            virtual void operator()(EA& ea) {
                using namespace ealib;
                using namespace ealib::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                datafile df("lod_knockouts.dat");
                df.add_field("lod_depth")
                .add_field("no_knockouts")
                .add_field("rx_knockedout")
                .add_field("location_knockedout");
                
                
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
                    int update_max = 1000;
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
        struct lod_gls_circle_square_plot : public ealib::analysis::unary_function<EA> {
            static const char* name() { return "lod_gls_circle_square_plot"; }
            
            virtual void operator()(EA& ea) {
                using namespace ealib;
                using namespace ealib::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                datafile df("lod_gls_circle_square_plot.dat");
                df.add_field("lod_depth");
                
                
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

        
        /*! lod_gls_germ_soma_mean_var reruns each subpopulation along a line of descent - setup for circle / square plot
         */
        template <typename EA>
        struct lod_gls_germ_soma_mean_var : public ealib::analysis::unary_function<EA> {
            static const char* name() { return "lod_gls_germ_soma_mean_var"; }
            
            virtual void operator()(EA& ea) {
                using namespace ealib;
                using namespace ealib::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                datafile df("lod_gls_germ_soma_mean_var.dat");
                df.add_field("lod_depth")
                .add_field("time_to_first_rep")
                .add_field("num_types_of_tasks")
                .add_field("num_germ")
                .add_field("num_pop")
                .add_field("germ_percent")
                .add_field("mean_germ_workload")
                .add_field("mean_germ_workload_var")
                .add_field("mean_soma_workload")
                .add_field("mean_soma_workload_var");
                
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
                    
                    double germ_count = 0;
                    double pop_count = 0;
                    accumulator_set<double, stats<tag::mean, tag::variance> > germ_workload_acc;
                    accumulator_set<double, stats<tag::mean, tag::variance> > soma_workload_acc;
                    
                    for(typename EA::individual_type::population_type::iterator j=control_ea->population().begin(); j!=control_ea->population().end(); ++j) {
                        typename EA::individual_type::individual_type& org=**j;
                        if (get<GERM_STATUS>(org, true)) {
                            ++germ_count;
                            germ_workload_acc(get<WORKLOAD>(org, 0.0));
                        } else {
                            soma_workload_acc(get<WORKLOAD>(org, 0.0));
                        }
                        ++pop_count;
                    }
                    
                    
                    // How many different types of tasks does the group do?
                    int task_type_count = 0;
                    if (get<TASK_NOT>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_NAND>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_AND>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_ORNOT>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_OR>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_ANDNOT>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_NOR>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_XOR>(*control_ea,0.0)) ++task_type_count;
                    if (get<TASK_EQUALS>(*control_ea,0.0)) ++task_type_count;

                    
                    
                    
                    double germ_percent = (germ_count/pop_count);
                    df.write(cur_update)
                    .write(task_type_count)
                    .write(germ_count)
                    .write(pop_count)
                    .write(germ_percent)
                    .write(mean(germ_workload_acc))
                    .write(variance(germ_workload_acc));
                    
                    if (germ_count != pop_count){
                        df.write(mean(soma_workload_acc))
                        .write(variance(soma_workload_acc));
                    } else {
                        df.write(0)
                        .write(0); 
                    }
                    
                    
                    df.endl();
                    
                    ++lod_depth;
                }
            }
            
            
                        
        };
        
        
        /*! lod_gls_task_count reruns each subpopulation along a line of descent - prints how many of each task were done.
         */
        template <typename EA>
        struct lod_gls_task_count : public ealib::analysis::unary_function<EA> {
            static const char* name() { return "lod_gls_task_count"; }
            
            virtual void operator()(EA& ea) {
                using namespace ealib;
                using namespace ealib::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                datafile df("lod_tasks.dat");
                df.add_field("lod_depth")
                .add_field("not")
                .add_field("nand")
                .add_field("and")
                .add_field("ornot")
                .add_field("or")
                .add_field("andnot")
                .add_field("nor")
                .add_field("xor")
                .add_field("equals");
                
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
                    
                   
                    // How many different types of tasks does the group do?
                    df.write(get<TASK_NOT>(*control_ea,0.0))
                    .write(get<TASK_NAND>(*control_ea,0.0))
                    .write(get<TASK_AND>(*control_ea,0.0))
                    .write(get<TASK_ORNOT>(*control_ea,0.0))
                    .write(get<TASK_OR>(*control_ea,0.0))
                    .write(get<TASK_ANDNOT>(*control_ea,0.0))
                    .write(get<TASK_NOR>(*control_ea,0.0))
                    .write(get<TASK_XOR>(*control_ea,0.0))
                    .write(get<TASK_EQUALS>(*control_ea,0.0));
                    
                    
                    df.endl();
                    
                    ++lod_depth;
                }
            }
            
            
            
        };
    }

    
    
}

#endif
