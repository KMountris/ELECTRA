/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_stimuli.hpp
   \brief ConfigStimuli class header file.
   \author Konstantinos A. Mountris
   \date 24/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_STIMULI_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_STIMULI_HPP_

#include "parser.hpp"

#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/conditions/body_loads/stimulus.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"


#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>
#include <set>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigStimuli
 * \brief Class to configure the stimuli of a simulation.
 */
class ConfigStimuli {

protected:

public:

    /**
     * \brief ConfigStimuli object constructor.
     */
    ConfigStimuli();


    /**
     * \brief ConfigStimuli object destructor.
     */
    virtual ~ConfigStimuli();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void SetStimuli(const Parser &parser, const std::string &body_type, const std::unordered_map<std::string, IMP::NodeSet> &node_sets,
                    int nodes_num, const ELECTRA::MeasureUnits &units, std::vector<ELECTRA::Stimulus> &stimuli, std::ostream &stream) const;

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_STIMULI_HPP_