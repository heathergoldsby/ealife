//
//  subpopulation_founder.h
//  ealife
//
//  Created by Heather Goldsby on 10/4/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//
#ifndef _EALIFE_SUBPOPULATION_FOUNDER_H_
#define _EALIFE_SUBPOPULATION_FOUNDER_H_


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>


using namespace ea;
using namespace boost::accumulators;


template <typename Individual>
class subpopulation_founder : public Individual {
public:
    typedef Individual base_type;
    typedef typename Individual::individual_type founder_type;
    typedef typename Individual::individual_ptr_type founder_ptr_type;

    
    //! Constructor.
    subpopulation_founder() : base_type() {
    }
    
    //! Copy constructor.
    subpopulation_founder(const subpopulation_founder& that) : base_type(that) {
        _founder = that._founder;
    }
    
    //! Assignment operator.
    subpopulation_founder& operator=(const subpopulation_founder& that) {
        if(this != &that) {
            base_type::operator=(that);
            _founder = that._founder;
        }
        return *this;
    }
    
    //! Destructor.
    virtual ~subpopulation_founder() {
    }
    
    founder_type& founder() { return _founder; }
    
protected:
    founder_type _founder;
    
private:
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::make_nvp("founder", _founder);
        ar & boost::serialization::make_nvp("individual", boost::serialization::base_object<base_type>(*this));
    }
};


/*! Chains together offspring and their parents, called for every inheritance event.
 */
template <typename EA>
struct founder_event : inheritance_event<EA> {
    
    //! Constructor.
    founder_event(EA& ea) : inheritance_event<EA>(ea) {
    }
    
    //! Destructor.
    virtual ~founder_event() {
    }
    
    //! Called for every inheritance event.
    virtual void operator()(typename EA::population_type& parents,
                            typename EA::individual_type& offspring,
                            EA& ea) {
        offspring.founder() = *offspring.population().front();
    }
};



#endif
