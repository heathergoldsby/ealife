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



namespace ea {
    namespace analysis {
        
        /*! Rerun a LoD, where a LoD is the sequential list of eas that contributed to a final dominant
         ea. Each ea is tagged with a founder, which is used for the rerun.
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
                    
                    // To replay, need to reset original configuration.
                    (*i)->md() = ea.md();
                    (*i)->rng().reset(get<RNG_SEED>(**i));
                    (*i)->initialize();
                    typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
                    (*i)->append(o);
                    
                    
                    int res_reset_inc = 500;
                    int res_resource_thresh = 500; 
                    int cur_update = 0;
                    // and run till the group amasses the right amount of resources
                    while ((get<GROUP_RESOURCE_UNITS>(**i) < 5000) && 
                           (cur_update < 10000)){
                        (*i)->update();
                        ++cur_update;
                        if (get<GROUP_RESOURCE_UNITS>(**i) > res_resource_thresh) {
                            df.write(cur_update); 
                            (*i)->env().reset_resources();
                            res_resource_thresh += res_reset_inc;
                        }
                    }
                    
                    
                    // catch missing data
                    int cur_res = get<GROUP_RESOURCE_UNITS>(**i); 
                    while (cur_res < 5000) {
                        df.write("99999"); 
                        cur_res += res_reset_inc;
                    }

                    df.endl();
                    
                   /* for (typename EA::individual_type::environment_type::location_ptr_type j=(*i)->env().locations().begin(); 
                     j!=(*i)->env().locations().end(); ++j) {
                     
                    }*/
                    ++lod_depth;
                }            
            }
            
            
            
            
       
            
            // analyze the constituents
            //                    for(typename EA::individual_type::population_type::iterator j=(*i)->population().begin(); j!=(*i)->population().end(); ++j) {
            //                        //                //    typename EA::individual_type::individual_type& org=**j;
            //                    
            //                    }
//
//
//    /*! Analyse final subpopulation variables
//     */
//    template <typename EA>
//    struct subpop_lod_gls : public ea::analysis::unary_function<EA> {
//        static const char* name() { return "lod_data"; }
//        
//        virtual void operator()(EA& ea) {
//            using namespace ea;
//            using namespace ea::analysis;
//            using namespace boost::accumulators;
//            
//            //line_of_descent<typename EA::individual_type> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea[0]);
//            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
//            //lod.runiq();
//            
//            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
//            
//            // skip def ancestor (that's what the +1 does)
//            for( ; i!=lod.end(); ++i) { 
//                // **i is the EA, AS OF THE TIME THAT IT DIED!
//                
//                // To replay, need to reset original configuration.
//                (*i)->md() = ea.md();
//                (*i)->rng().reset(get<RNG_SEED>(**i));
//                (*i)->initialize();
//                typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
//                (*i)->append(o);
//                
//                
//                // and run
//                (*i)->update();
//                
//                
//                // need to iterate through locations, not individuals.
//                //for(typename EA::individual_type::population_type::iterator j=(*i)->population().begin(); j!=(*i)->population().end(); ++j) {
//                //    typename EA::individual_type::individual_type& org=**j;
//                
//                // WORKLOAD // GS status
//                
//            }
//            
//        }
//        
//        //for(; i!=lod.end(); ++i, ++j, ++k) {
//        // rerun here...
//        //    (*i)->update();
//        //}
//        
//    };

        
        };
    };
}
        
#endif
