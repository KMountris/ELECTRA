/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_TPP_

#include "config_conduct_sys.hpp"

namespace APP_ELECTRA
{


//!  PUBLIC  //
template<short DIM, short CELL_NODES>
ConfigConductSys<DIM,CELL_NODES>::ConfigConductSys()
{}


template<short DIM, short CELL_NODES>
ConfigConductSys<DIM,CELL_NODES>::~ConfigConductSys()
{}


template<short DIM, short CELL_NODES>
void ConfigConductSys<DIM,CELL_NODES>::SetConductSystem(const Parser &parser, const ELECTRA::MeasureUnits &units,
    const std::vector<IMP::Vec<int(DIM), double>> &tissue_nodes, ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const
{
    this->SetCsGeometry(parser, units, tissue_nodes, conduct_sys, stream);
    this->SetCsDiffusivity(parser, units, conduct_sys, stream);
}


//!  PROTECTED  //
template<short DIM, short CELL_NODES>
void ConfigConductSys<DIM,CELL_NODES>::SetCsGeometry(const Parser &parser, const ELECTRA::MeasureUnits &units,
        const std::vector<IMP::Vec<DIM, double>> &tissue_nodes, ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const
{
    // Get units for conduction system geometry.
    std::string cs_geo_unit = parser.GetValue<std::string>("conduction system.geometry.unit");
    if (cs_geo_unit != parser.GetValue<std::string>("tissue.geometry.unit")) {
        std::string tissue_geo_unit = parser.GetValue<std::string>("tissue.geometry.unit");
        std::string error_msg = "Could not configure conduction system. The conduction system length unit ["+cs_geo_unit+"] is different from the tissue length unit["+tissue_geo_unit+"].";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
    }

    // Load the conduction system tree.
    stream << ELECTRA::Logger::Message("Loading conduction system tree geometry... ");
    std::string cs_filename = parser.ResolvePath(parser.GetValue<std::string>("conduction system.geometry.file"));
    conduct_sys.LoadTree(cs_filename);
    stream << "OK\n";

    // Set the conduction system nodesets.
    conduct_sys.SetAvNodeId(parser.GetValue<std::string>("conduction system.nodesets.av node"));
    conduct_sys.SetTerminalNodeIds(parser.GetValue<std::string>("conduction system.nodesets.terminal nodes"));

    // Compute the Purkinje-myocardial junctions.
    stream << ELECTRA::Logger::Message("Computing Purkinje-myocardial junctions... ");
    conduct_sys.ComputePmjs(tissue_nodes, parser.GetValue<double>("conduction system.geometry.pmj radius")*units[cs_geo_unit]);
    stream << "OK\n";

}


template<short DIM, short CELL_NODES>
void ConfigConductSys<DIM,CELL_NODES>::SetCsDiffusivity(const Parser &parser, const ELECTRA::MeasureUnits &units,
        ELECTRA::ConductionSystem<DIM> &conduct_sys, std::ostream &stream) const
{

    // Obtain diffusivity values of the conduction system's tree.
    stream << ELECTRA::Logger::Message("Obtaining diffusivity data for conduction system... ");
    std::vector<double> diff_values;
    if (parser.IsSingleValue("conduction system.diffusivity.tree diffusivity")) {
        diff_values.clear();
        diff_values.assign(conduct_sys.NodesNum(), parser.GetValue<double>("conduction system.diffusivity.tree diffusivity"));

    } else if (parser.IsArray("conduction system.diffusivity.tree diffusivity")) {
        auto object = parser.GetObject("conduction system.diffusivity.tree diffusivity");
        diff_values = object.get<std::vector<double>>();

        // Check if the number of the parsed diffusivity values is the expected.
        if (static_cast<int>(diff_values.size()) != conduct_sys.NodesNum()) {
            std::string error_msg = "Could not assign diffusivity to conduction system. The conduction system has [" + std::to_string(conduct_sys.NodesNum()) +
                                    "] and the parsed diffusivity values are [" + std::to_string(diff_values.size()) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
        }
    } else if (parser.IsMultiArray("conduction system.diffusivity.tree diffusivity")) {
        this->ObtainValuesFromNodesets(parser, "conduction system.diffusivity.tree diffusivity", conduct_sys.NodeSets(), conduct_sys.NodesNum(), diff_values);
    } else {
        std::string error_msg = "Could not assign diffusivity values to conduction system. Diffusivity attribute has unexpected format.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
    }

    // Get units of diffusivity.
    std::string nominator = parser.GetValue<std::string>("conduction system.diffusivity.unit");
    std::string denominator = "";
    std::size_t pos = nominator.find('/');
    if (pos != std::string::npos) {
        denominator = nominator.substr(pos+1);
        nominator = nominator.substr(0, pos);
    } else {
        throw std::invalid_argument(ELECTRA::Logger::Error("Could not extract conduction system diffusivity unit. Expected unit format: area / time"));
    }

    // Apply diffusivity unit.
    std::transform(std::begin(diff_values), std::end(diff_values), std::begin(diff_values),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, units[nominator]/units[denominator]));

    // Set the conduction system tree diffusivity.
    conduct_sys.SetDiffusionCoeffs(diff_values);
    conduct_sys.SetPmjDiffusionCoeffs(parser.GetValue<double>("conduction system.diffusivity.pmj diffusivity") * (units[nominator]/units[denominator]));
    stream << "OK\n";

    // Compute the diffusivity decay at the terminal nodes of the conduction system.
    conduct_sys.ComputeDiffuseTransition(parser.GetValue<double>("conduction system.diffusivity.sigmoid kappa"));

}


template<short DIM, short CELL_NODES>
void ConfigConductSys<DIM, CELL_NODES>::ObtainValuesFromNodesets(const Parser &parser, const std::string &attribute,
            const std::unordered_map<std::string, IMP::NodeSet> &nodesets, int values_num, std::vector<double> &values) const
{
    values.clear();
    values.resize(values_num, 0.);
    int counter = 0;
    double value = 0.;
    auto object = parser.GetObject(attribute);
    if (object[0][0].is_string() && object[0][1].is_string()) {
        // Case where values are given as strings together with nodeset names.
        auto entries = object.get<std::vector<std::pair<std::string,std::string>>>();

        // Assign diffusivity values to the nodes of each nodeset.
        for (const auto &entry : entries) {
            value = std::stod(entry.first);
            for (const auto &id : nodesets.at(entry.second).NodeIds()) {
                values[id] = value;
                counter++;
            }
        }

        // Check that all nodes obtained a diffusivity value.
        if (counter != values_num) {
            std::string error_str = "Could not obtain values from node sets correctly. The total number of nodes is [" + std::to_string(values_num) +
                                    "] while the number of the defined nodes in the node sets is [" + std::to_string(counter) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }

    } else if (object[0][0].is_number() && object[0][1].is_string()) {
        // Case where values are given as numbers together with nodeset names.
        auto entries = object.get<std::vector<std::pair<double,std::string>>>();

        // Assign diffusivity values to the nodes of each nodeset.
        for (const auto &entry : entries) {
            value = entry.first;
            for (const auto &id : nodesets.at(entry.second).NodeIds()) {
                values[id] = value;
                counter++;
            }
        }

        // Check that all nodes obtained a diffusivity value.
        if (counter != values_num) {
            std::string error_str = "Could not obtain values from node sets correctly. The total number of nodes is [" + std::to_string(values_num) +
                                    "] while the number of the defined nodes in the node sets is [" + std::to_string(counter) + "]";
            throw std::runtime_error(ELECTRA::Logger::Error(error_str));
        }
    }

}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_CONDUCT_SYS_TPP_