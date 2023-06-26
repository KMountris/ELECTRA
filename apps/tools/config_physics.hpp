/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_physics.hpp
   \brief ConfigPhysics class header file.
   \author Konstantinos A. Mountris
   \date 25/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/Physics"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <string>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigGeo
 * \brief Class to configure the physics of a simulation.
 */
template <short DIM, short CELL_NODES>
class ConfigPhysics {

private:

    std::unordered_map<std::string, ELECTRA::ReactionDiffusionType> react_diff_map_;        /**< Map of the reaction diffusion model types */


protected:

    void SetBidomain(const Parser &parser, const ELECTRA::MeasureUnits &units,
            std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const;


    void SetMonodomain(const Parser &parser, const ELECTRA::MeasureUnits &units,
            std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const;


public:

    /**
     * \brief ConfigPhysics object constructor.
     */
    ConfigPhysics();


    /**
     * @brief Destroy the Config Physics object
     * 
     */
    virtual ~ConfigPhysics();


    /**
     * @brief 
     * 
     * @param parser 
     * @param react_diff 
     * @param stream 
     */
    void InitializeReactionDiffusion(const Parser &parser, std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff) const;


    /**
     * @brief Set the Reaction Diffusion object
     * 
     * @param parser 
     * @param units 
     * @param react_diff 
     * @param stream 
     */
    void SetReactionDiffusion(const Parser &parser, const ELECTRA::MeasureUnits &units,
            std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_HPP_

#include "config_physics.tpp"