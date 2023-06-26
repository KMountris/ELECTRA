/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_sim.hpp
   \brief ConfigSim class header file.
   \author Konstantinos A. Mountris
   \date 12/05/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_SIM_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_SIM_HPP_

#include "config_approximation.hpp"
#include "config_conduct_sys.hpp"
#include "config_electric.hpp"
#include "config_electrophys.hpp"
#include "config_geo.hpp"
#include "config_output.hpp"
#include "config_physics.hpp"
#include "config_post_process.hpp"
#include "config_stimuli.hpp"
#include "config_units.hpp"
#include "parser.hpp"

#include <CLOUDEA/CLOUDEA>
#include "ELECTRA/ELECTRA"

#include <termcolor/termcolor.hpp>

#include <string>
#include <vector>
#include <algorithm>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigSim
 * \brief Class to configure and execute a simulation wth the ElectraSim app.
 * \tparam DIM The dimensions of the geometry model of the simulation.
 * \tparam CELL_NODES The number of nodes of the geometry model's cells.
 */
template<short DIM, short CELL_NODES>
class ConfigSim {

public:

    /**
     * \brief ConfigSim object constructor.
     */
    inline ConfigSim();


    /**
     * \brief ConfigSim object destructor.
     */
    inline virtual ~ConfigSim();


    /**
     * \brief Check if a valid simulation file has been provided.
     * \param [in] parser The parser of the simulation file.
     * \return [void]
     */
    inline void CheckValid(const Parser &parser);


    /**
     * \brief Launch a electrophysiology simulation for a tissue.
     * Both 2D and 3D tissues are supported.
     * \param [in] parser The parser of the simulation file.
     * \param [out] stream The output logging stream.
     * \return [void].
     */
    inline void Tissue(const Parser &parser, std::ostream &stream);

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_SIM_HPP_

#include "config_sim.tpp"