/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_pre.hpp
   \brief ConfigPre class header file.
   \author Konstantinos A. Mountris
   \date 02/07/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_PRE_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_PRE_HPP_

#include "parser.hpp"

#include "ELECTRA/ELECTRA"

#include <termcolor/termcolor.hpp>

#include "config_approximation.hpp"
#include "config_geo.hpp"
#include "config_fibers.hpp"
#include "config_output.hpp"

#include <string>
#include <vector>
#include <algorithm>

using namespace ELECTRA;

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigPre
 * \brief Class to configure and execute preprocessing tasks wth the ElectraPre app.
 * \tparam DIM The dimensions of the geometry model to be preprocessed.
 * \tparam CELL_NODES The number of nodes of the geometry model's cells.
 */
template<short DIM, short CELL_NODES>
class ConfigPre {

public:

    /**
     * \brief ConfigPre object constructor.
     */
    inline ConfigPre();


    /**
     * \brief ConfigPre object destructor.
     */
    inline virtual ~ConfigPre();


    /**
     * \brief Check if a valid simulation file has been provided.
     * \param [in] parser The parser of the simulation file.
     * \return [void]
     */
    inline void CheckValid(const Parser &parser);


    /**
     * \brief
     * \param [in] parser The parser of the preprocessing configuration file.
     * \param [out] stream The output logging stream.
     * \return [void].
     */
    inline void Preprocess(const Parser &parser, std::ostream &stream);


    /**
     * \brief Compute the fibers distribution in an atrial or ventricular 3D geometry.
     * \param [in] parser The parser of the preprocessing configuration file.
     * \param [out] stream The output logging stream.
     * \return [void].
     */
    inline void Fibers(const Parser &parser, std::ostream &stream);

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_PRE_HPP_

#include "config_pre.tpp"