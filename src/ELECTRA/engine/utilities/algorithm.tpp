/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#ifndef ELECTRA_UTILITIES_ALGORITHM_TPP_
#define ELECTRA_UTILITIES_ALGORITHM_TPP_


#include "ELECTRA/engine/utilities/algorithm.hpp"

namespace ELECTRA {

namespace ALGORITHM {


template<typename TYPE>
TYPE RushLarsen(TYPE val_steady, TYPE val_old, TYPE dt, TYPE time_steady)
{
    return val_steady - (val_steady - val_old) * std::exp(-dt/time_steady);
}


template<typename TYPE>
TYPE ForwardEuler(TYPE func_old, TYPE dt, TYPE func_derivative)
{
    return func_old + dt*func_derivative;
}


template<typename TYPE>
TYPE KahanSum(std::initializer_list<TYPE> values) {

    // Initialize sum and error.
    TYPE sum = static_cast<TYPE>(0.);
    TYPE c = static_cast<TYPE>(0.);

    for (auto val : values) {
        TYPE y = val - c;
        TYPE t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}


} //end of namespace ALGORITHM


} //end of namespace ELECTRA

#endif //ELECTRA_UTILITIES_ALGORITHM_TPP_