/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_fibers.hpp
   \brief ConfigFibers class header file.
   \author Konstantinos A. Mountris
   \date 02/07/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_FIBERS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_FIBERS_HPP_

#include "parser.hpp"

#include <IMP/IMP>
#include <CLOUDEA/CLOUDEA>

#include "ELECTRA/Utilities"
#include "ELECTRA/Fibers"

#include <termcolor/termcolor.hpp>

#include <string>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <algorithm>

using namespace ELECTRA;

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigFibers
 * \brief Class to configure the preprocessing of cardiac fibers.
 */
template <short DIM, short CELL_NODES>
class ConfigFibers {

public:

    /**
     * \brief ConfigFibers object constructor.
     */
    ConfigFibers();


    virtual ~ConfigFibers();


    void SetVentriFiberRules(const Parser &parser, VentriFiberRules &ventri_fiber_rules, std::ostream &stream) const;


    void SetVentriTags(const Parser &parser, VentriTags &ventri_tags, std::ostream &stream) const;


    void SetAtriFiberRules(const Parser &parser, AtriFiberRules &atri_fiber_rules, std::ostream &stream) const;


    void SetAtriTags(const Parser &parser, AtriTags &atri_tags, std::ostream &stream) const;


    void ComputeVentriFibers(const Parser &parser, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
                             const CLOUDEA::Fpm<DIM> &fpm, const VentriTags &tags, const VentriFiberRules &rules,
                             Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const;


    void ComputeAtriFibers(const Parser &parser, const IMP::Mesh<DIM,CELL_NODES> &mesh, const IMP::Voronoi<DIM> &voro,
                           const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &tags, const AtriFiberRules &rules,
                           Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_FIBERS_HPP_

#include "config_fibers.tpp"