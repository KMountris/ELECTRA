/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_units.hpp
   \brief ConfigUnits class header file.
   \author Konstantinos A. Mountris
   \date 13/05/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_UNITS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_UNITS_HPP_

#include "parser.hpp"

#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"

#include <termcolor/termcolor.hpp>

#include <string>
#include <iostream>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigUnits
 * \brief Class to configure the measure units of a simulation.
 */
class ConfigUnits {

public:

    /**
     * \brief ConfigUnits object constructor.
     */
    ConfigUnits();


    /**
     * \brief ConfigUnits object destructor.
     */
    virtual ~ConfigUnits();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void SetReferenceScale(const Parser &parser, ELECTRA::MeasureUnits &units, std::ostream &stream) const;

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_UNITS_HPP_