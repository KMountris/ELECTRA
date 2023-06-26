/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_conduct_sys.hpp
   \brief ConfigConductSys class header file.
   \author Konstantinos A. Mountris
   \date 20/10/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/ConductionSystem"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <iostream>
#include <unordered_map>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigConductSys
 * \brief Class to configure the conduction system of the heart.
 */
template <short DIM, short CELL_NODES>
class ConfigConductSys {

protected:

    void SetCsGeometry(const Parser &parser, const ELECTRA::MeasureUnits &units, const std::vector<IMP::Vec<DIM, double>> &tissue_nodes,
            ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const;


    void SetCsDiffusivity(const Parser &parser, const ELECTRA::MeasureUnits &units, ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const;


    void ObtainValuesFromNodesets(const Parser &parser, const std::string &attribute,
            const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int values_num, std::vector<double> &values) const;

public:

    /**
     * \brief ConfigConductSys object constructor.
     */
    ConfigConductSys();


    virtual ~ConfigConductSys();


    void SetConductSystem(const Parser &parser, const ELECTRA::MeasureUnits &units,
            const std::vector<IMP::Vec<int(DIM), double>> &tissue_nodes, ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_HPP_

#include "config_conduct_sys.tpp"