/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file algorithm.hpp
   \brief General algorithm collection header file.
   \author Konstantinos A. Mountris
   \date 20/04/2019
*/

#ifndef ELECTRA_UTILITIES_ALGORITHM_HPP_
#define ELECTRA_UTILITIES_ALGORITHM_HPP_


#include <cmath>
#include <initializer_list>
#include <string>
#include <iostream>
#include <cctype>
#include <algorithm>
		

namespace ELECTRA {

/** \addtogroup Utilities \{ */


/**
 * \namespace ELECTRA::ALGORITHM
 * \brief Collects general algorithmic functions for mathematical manipulations and more.
 */
namespace ALGORITHM {


/**
 * \brief Temporal integration of differential equation using the Rush-Larsen method \cite rush1978. 
 * \tparam TYPE The data type. A standard floating point type should be preferred. 
 * \param [in] val_steady The value of the differential equation at steady state. 
 * \param [in] val_old The value of the differential equation at the previous time step.
 * \param [in] dt The time increment between time steps.
 * \param [in] time_steady The final time at steady state.
 * \return [TYPE] The new value of the differential equation at the next time step.
 */
template<typename TYPE>
inline TYPE RushLarsen(TYPE val_steady, TYPE val_old, TYPE dt, TYPE time_steady);


/**
 * \brief Temporal integration of differential equation using the Forward Euler method.
 * \tparam TYPE The data type. A standard floating point type should be preferred.
 * \param [in] func_old The value of the differential equation at the previous time step.
 * \param [in] dt The time increment between time steps.
 * \param [in] func_derivative The value of the used function to update the value of the previous time step.
 * \return [TYPE] The new value of the differential equation at the next time step.
 */
template<typename TYPE>
inline TYPE ForwardEuler(TYPE func_old, TYPE dt, TYPE func_derivative);


/**
 * \brief Kahan summation algorithm for precise summation \cite kahan1965. 
 * \tparam TYPE The data type. A standard floating point type should be preferred.
 * \param [in] values The list of values to be summed.
 * \return [TYPE] The sum of the given values. 
 */
template<typename TYPE>
inline TYPE KahanSum(std::initializer_list<TYPE> values);


/**
 * \brief 
 * \param str1 
 * \param str2 
 * \return true 
 * \return false 
 */
bool StringCompCaseInsensitive(const std::string &str1, const std::string &str2);


/**
 * \brief 
 * \param phrase 
 * \param word 
 * \return true 
 * \return false 
 */
bool ExistExactWord(const std::string &phrase, const std::string &word);


bool ExistWord(const std::string &phrase, const std::string &word);


} //end of namespace ALGORITHM



/** \} End of Doxygen Groups*/

} //end of namespace ELECTRA

#endif //ELECTRA_UTILITIES_ALGORITHM_HPP_

#include "ELECTRA/engine/utilities/algorithm.tpp"