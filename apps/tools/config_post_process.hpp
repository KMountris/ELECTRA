/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_post_process.hpp
   \brief ConfigPostProcess class header file.
   \author Konstantinos A. Mountris
   \date 27/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_POST_PROCESS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_POST_PROCESS_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/Physics"
#include "ELECTRA/PostProcess"


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
 * \class ConfigPostProcess
 * \brief Class to configure the post processing of simulation results.
 */
template<short DIM, short CELL_NODES>
class ConfigPostProcess {

protected:

public:

    /**
     * \brief ConfigPostProcess object constructor.
     */
    ConfigPostProcess();


    /**
     * \brief ConfigPostProcess object destructor.
     */
    virtual ~ConfigPostProcess();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void PerformPostProcess(const Parser &parser, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff,
                            const std::unordered_map<std::string, IMP::NodeSet> &node_sets, ELECTRA::PostProcess &post_process, std::ostream &stream) const;

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#include "config_post_process.tpp"

#endif //ELECTRA_APPS_TOOLS_CONFIG_POST_PROCESS_HPP_