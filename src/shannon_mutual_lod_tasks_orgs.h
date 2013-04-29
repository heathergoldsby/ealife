//
//  shannon_mutual_lod_tasks_orgs.h
//  ealife
//
//  Created by Heather Goldsby on 4/16/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_SHANNON_MUTUAL_LOD_TASKS_ORGS_H
#define _EALIFE_SHANNON_MUTUAL_LOD_TASKS_ORGS_H

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
        
        /* Rerun a LoD, track shannon mutual information among tasks and individuals along LoD. */
        
        
        
        
        /*! lod_shannon_tasks_orgs 
         
         This particular measurement is the same one we used in the PNAS paper. It does not take into account how much time
         an org spends on a particular activity.
         
         Also, this configuration assumes orgs cannot replicate over one another.
         
         pi = 1 / # active orgs
         pj = # of times task done / # tasks done total 
         pij = (# times org did task j / # tasks done by org) / # active orgs
         
         */
        template <typename EA>
        struct lod_shannon_tasks_orgs : public ealib::analysis::unary_function<EA> {
            static const char* name() { return "lod_shannon_tasks_orgs"; }
            
            virtual void operator()(EA& ea) {
                using namespace ealib;
                using namespace ealib::analysis;
                
                line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
                
                typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
                
                
                datafile df("lod_shannon_tasks_orgs.dat");
                df.add_field("lod_depth")
                .add_field("shannon")
                .add_field("shannon norm");
                

                
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
                    

                    // replay! till the group amasses the right amount of resources
                    // or exceeds its window...
                    int cur_update = 0;
                    int update_max = 10000;
                    // and run till the group amasses the right amount of resources
                    while ((get<GROUP_RESOURCE_UNITS>(**i,0) < get<GROUP_REP_THRESHOLD>(**i)) &&
                           (cur_update < update_max)){
                        (*i)->update();
                        ++cur_update;
                    }
                    
                    df.write(lod_depth);
                    
                    std::vector< std::vector<double> > pij;
                    std::vector<double> pj (9);
                    double pi = 0.0;
                    double pop_count = 0;
                    double active_pop = 0;

                    
                    // cycle through orgs and create matrix for shannon mutual information.
                    for(typename EA::individual_type::population_type::iterator j=(*i)->population().begin(); j!=(*i)->population().end(); ++j) {
                        typename EA::individual_type::individual_type& org=**j;
                        ++pop_count;
                        std::vector<double> porg (9);
                        porg[0] = get<TASK_NOT>(org,0.0);
                        porg[1] = get<TASK_NAND>(org,0.0);
                        porg[2] = get<TASK_AND>(org,0.0);
                        porg[3] = get<TASK_ORNOT>(org,0.0);
                        porg[4] = get<TASK_OR>(org,0.0);
                        porg[5] = get<TASK_ANDNOT>(org,0.0);
                        porg[6] = get<TASK_NOR>(org,0.0);
                        porg[7] = get<TASK_XOR>(org,0.0);
                        porg[8] = get<TASK_EQUALS>(org,0.0);
                        
                        double total_num_tasks = std::accumulate(porg.begin(), porg.end(), 0);
                        
                        // Normalize the tasks and add to matrix
                        if(total_num_tasks > 0) {
                            for (int k=0; k<porg.size(); ++k) {
                                porg[k] /= total_num_tasks;
                            }
                            ++active_pop;
                            pij.push_back(porg);
                        }
                    }
                    
                    double shannon_sum = 0.0;
                    double shannon_norm = 0.0;
                    if (active_pop > 1) {
                    
                        // figure out pj
                        for (int k=0; k<pj.size(); ++k) {
                            for (int m=0; m<active_pop; ++m) {
                                pj[k] += pij[m][k];
                            }
                            pj[k] /= active_pop;
                        }

                        // compute shannon mutual information based on matrix...

                        double shannon_change = 0.0;
                        double t_pij = 0.0;
                        double t_pi = 1.0/active_pop;
                        double t_pj = 0;
                        double pij_sum = 0.0;
                        // calculate shannon mutual information
                        for (int i=0; i<active_pop; i++) {
                            for (int j=0; j<pj.size(); j++) {
                                t_pij = pij[i][j]/active_pop;
                                t_pj = pj[j];
                                pij_sum += t_pij;
                                if (t_pi && t_pj && t_pij) {
                                    shannon_change= (t_pij * log(t_pij / (t_pi * t_pj)));
                                    shannon_sum += shannon_change;
                                }
                            }
                        }
                        shannon_norm = shannon_sum / log((double)active_pop);

                    }
                    df.write(shannon_sum)
                    .write(shannon_norm);
                    df.endl();
                     
                     
                    ++lod_depth;
                }
 
            }
            
        };
    }
    
    
    
}


#endif
