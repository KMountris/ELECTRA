/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
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
 *
 */

/**
   \file measure_units.hpp
   \brief MeasureUnits class header file.
   \author Konstantinos A. Mountris
   \date 11/03/2019
*/

#ifndef ELECTRA_UTILITIES_MEASURE_UNITS_HPP_
#define ELECTRA_UTILITIES_MEASURE_UNITS_HPP_

#include <string>
#include <unordered_map>
#include <iostream>
#include <stdexcept>
#include <exception>


namespace ELECTRA {

/** \addtogroup Utilities \{ */


/**
 * \class MeasureUnits
 * \brief Class implemmenting measure units conversion on pre-defined reference units in the SI system.
 * \author Konstantinos A. Mountris
 */
class MeasureUnits
{

private:

    double ref_time_si_value_;                          /**< The reference value of the time unit with respect to the SI system */

    double ref_length_si_value_;                        /**< The reference value of the length unit with respect to the SI system */

    double ref_conductance_si_value_;                   /**< The reference value of the conductance unit with respect to the SI system */

    double ref_capacitance_si_value_;                   /**< The reference value of the capacitance unit with respect to the SI system */

    double ref_current_si_value_;                       /**< The reference value of the current unit with respect to the SI system */

    std::unordered_map<std::string, double> units_;     /**< The measure units values with respect to the available reference unit values */


protected:

    /**
     * \brief Set the corresponding values of all the measure units considering the given reference values.
     * 
     * Reference unit values by default are expressed in SI system (i.e., length=m, conductance=S, capacitance=F, etc.)
     */
    void SetTimeUnits() noexcept;

    void SetLengthUnits() noexcept;

    void SetConductanceUnits() noexcept;

    void SetCapacitanceUnits() noexcept;

    void SetCurrentUnits() noexcept;


public:

    /**
     * \brief The MeasureUnits default constructor. 
     */
    MeasureUnits();


    /**
     * \brief The MeasureUnits destructor.
     * 
     */
    virtual ~MeasureUnits();


    /**
     * \brief Set the reference time unit value with respect to the SI system (e.g., s=1., ds=0.1, ms=0.01, etc.).
     * \param [in] val The value of the reference time unit with respect to the SI system.
     */
    void SetRefTimeSIValue(double val) noexcept;


    /**
     * \brief Set the reference length unit value with respect to the SI system (e.g., m=1., dm=0.1, cm=0.01, etc.).
     * \param [in] val The value of the reference length unit with respect to the SI system.
     */
    void SetRefLengthSIValue(double val) noexcept;


    /**
     * \brief Set the reference conductance unit value with respect to the SI system (e.g., S=1., dS=0.1, cS=0.01, etc.).
     * \param [in] val The value of the reference conductance unit with respect to the SI system.
     */
    void SetRefConductanceSIValue(double val) noexcept;


    /**
     * \brief Set the reference capacitance unit value with respect to the SI system (e.g., F=1., dF=0.1, cF=0.01, etc.).
     * \param [in] val The value of the reference capacitance unit with respect to the SI system.
     */
    void SetRefCapacitanceSIValue(double val) noexcept;


    /**
     * \brief Set the reference current unit value with respect to the SI system (e.g., A=1., dA=0.1, cA=0.01, etc.).
     * \param [in] val The value of the reference current unit with respect to the SI system.
     */
    void SetRefCurrentSIValue(double val) noexcept;


    /**
     * \brief Operator to access values from the units container by key. 
     * 
     * \note The returned value is expressed with respect to the reference unit values.
     *       For example: MeasureUnits["cm"] = 0.01 if the reference length unit value is in meters
     *       and MeasureUnits["cm"] = 1 if the reference length unit value is in centimeters.
     * 
     * \param [in] key The name of the unit to retrieve its value. 
     * \return [const double&] The value of the unit corresponding to the given key. 
     */
    const double & operator [] (const std::string &key) const;

};



/** \} End of Doxygen Groups */

} // End of namespace ELECTRA


#endif // ELECTRA_UTILITIES_MEASURE_UNITS_HPP_