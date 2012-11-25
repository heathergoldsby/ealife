//
//  multi_birth_selfrep_not_nand_ancestor.h
//  ealife
//
//  Created by Heather Goldsby on 10/16/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef _EALIFE_MULTI_BIRTH_SELFREP_NOT_NAND_ANCESTOR_H_
#define _EALIFE_MULTI_BIRTH_SELFREP_NOT_NAND_ANCESTOR_H_


#include <ea/digital_evolution.h>
#include <ea/meta_data.h>


namespace ea {
    /*! Generates a self-replicating ancestor that performs not.
     */
    struct multibirth_selfrep_not_nand_ancestor {
        template <typename EA>
        typename EA::representation_type operator()(EA& ea) {
            typename EA::representation_type repr;
            repr.resize(get<REPRESENTATION_SIZE>(ea));
            std::fill(repr.begin(), repr.end(), ea.isa()["nop_x"]);
            
            // Must use representation size of 100.
            assert(repr.size() == 100);
            
            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
            repr[1] =  ea.isa()["nop_c"]; // nopc
            repr[2] =  ea.isa()["nop_a"]; // nopa
            repr[3] =  ea.isa()["h_search"]; // hsearch
            repr[4] =  ea.isa()["nop_c"]; // nopc
            repr[5] =  ea.isa()["mov_head"]; // movhead
            
            // not
            repr[24] = ea.isa()["fixed_input"]; // input
            repr[25] = ea.isa()["fixed_input"]; // input
            repr[26] = ea.isa()["push"]; // push
            repr[27] = ea.isa()["nop_c"]; // nopc
            repr[28] = ea.isa()["pop"]; // pop
            repr[29] = ea.isa()["nand"]; // nand
            repr[30] = ea.isa()["output"]; //output
            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            // nand
            repr[74] = ea.isa()["fixed_input"]; // input
            repr[75] = ea.isa()["nop_c"]; // nopc
            repr[76] = ea.isa()["fixed_input"]; // input
            repr[77] = ea.isa()["nand"]; // nand
            repr[78] = ea.isa()["output"]; //output
            repr[79] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
            
            repr[91] =  ea.isa()["h_search"]; // hsearch
            repr[92] =  ea.isa()["h_copy"]; // hcopy
            repr[93] =  ea.isa()["nop_c"]; // nopc
            repr[94] =  ea.isa()["nop_a"]; // nopa
            repr[95] =  ea.isa()["if_label"]; // iflabel
            repr[96] =  ea.isa()["h_divide_soft_parent_reset"]; // hdivide
            repr[97] =  ea.isa()["mov_head"]; // movhead
            repr[98] =  ea.isa()["nop_a"]; // nopa
            repr[99] =  ea.isa()["nop_b"]; // nopb
            
//            repr[0] =  ea.isa()["h_alloc"]; // h_alloc
//            repr[1] =  ea.isa()["nop_c"]; // nopc
//            repr[2] =  ea.isa()["nop_a"]; // nopa
//            repr[3] =  ea.isa()["h_search2"]; // hsearch - finds a, b
//            repr[4] =  ea.isa()["nop_c"]; // nopc
//            repr[5] =  ea.isa()["mov_head"]; // movhead
//            
//            // not
//            repr[24] = ea.isa()["fixed_input"]; // input
//            repr[25] = ea.isa()["fixed_input"]; // input
//            repr[26] = ea.isa()["push"]; // push
//            repr[27] = ea.isa()["nop_c"]; // nopc
//            repr[28] = ea.isa()["pop"]; // pop
//            repr[29] = ea.isa()["nand"]; // nand
//            repr[30] = ea.isa()["output"]; //output
//            repr[31] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
//            
//            // nand
//            repr[74] = ea.isa()["fixed_input"]; // input
//            repr[75] = ea.isa()["nop_c"]; // nopc
//            repr[76] = ea.isa()["fixed_input"]; // input
//            repr[77] = ea.isa()["nand"]; // nand
//            repr[78] = ea.isa()["output"]; //output
//            repr[79] = ea.isa()["donate_res_to_group"]; // donate_res_to_group
//            
//            repr[85] =  ea.isa()["h_search2"]; // hsearch // no complement
//            repr[86] =  ea.isa()["h_copy"]; // hcopy
//            repr[87] =  ea.isa()["nop_c"]; // nopc
//            repr[88] =  ea.isa()["nop_a"]; // nopa
//            repr[89] =  ea.isa()["if_label"]; // iflabel
//            repr[90] =  ea.isa()["nop_a"];
//            repr[91] =  ea.isa()["h_search2"];
//            repr[92] =  ea.isa()["nop_a"];
//            repr[93] =  ea.isa()["nop_c"];
//            repr[94] =  ea.isa()["mov_head"]; // movhead
//            repr[95] =  ea.isa()["nop_b"];
//            repr[96] =  ea.isa()["nop_x"];
//            repr[97] =  ea.isa()["h_divide_soft_parent_reset"]; // hdivide
//            repr[98] =  ea.isa()["nop_a"]; // nopa
//            repr[99] =  ea.isa()["nop_b"]; // nopa

            // if this doesnt work flip the hdivide and last nopa
            return repr;
        }
    };
    
}

#endif
