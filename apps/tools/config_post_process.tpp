/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_CONFIG_POST_PROCESS_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_POST_PROCESS_TPP_

#include "config_post_process.hpp"

namespace APP_ELECTRA
{

template<short DIM, short CELL_NODES>
ConfigPostProcess<DIM, CELL_NODES>::ConfigPostProcess()
{}


template<short DIM, short CELL_NODES>
ConfigPostProcess<DIM, CELL_NODES>::~ConfigPostProcess()
{}


template<short DIM, short CELL_NODES>
void ConfigPostProcess<DIM, CELL_NODES>::PerformPostProcess(const Parser &parser, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff,
        const std::unordered_map<std::string, IMP::NodeSet> &node_sets, ELECTRA::PostProcess &post_process, std::ostream &stream) const
{
    // Requested post processing actions.
    stream << "Requested post processing: \n";

    // Apply action potential generation post processing.
    if (parser.HasAttribute("post process.action potential")) {
        std::string ep_nset_name = parser.GetValue<std::string>("post process.action potential");

        if (ELECTRA::ALGORITHM::StringCompCaseInsensitive(ep_nset_name, "all") ||
            ELECTRA::ALGORITHM::StringCompCaseInsensitive(ep_nset_name, "")) {

            std::vector<int> ep_node_ids(react_diff->Vout(0).size());
            std::iota (std::begin(ep_node_ids), std::end(ep_node_ids), 0);
            post_process.ComputeActionPotential(react_diff->Vout(), ep_node_ids, react_diff->Dt()*react_diff->OutputSteps());
        } else {
            post_process.ComputeActionPotential(react_diff->Vout(), node_sets.at(ep_nset_name).NodeIds(), react_diff->Dt()*react_diff->OutputSteps());
        }

        stream << "        - action potential computation \n";
    } // End of Apply action potential generation post processing.

    // Apply local activation time generation.
    if (parser.HasAttribute("post process.local activation time")) {

        // Get the LAT time range.
        double time_start = parser.GetValue<double>("post process.local activation time.time start");
        double time_end = parser.GetValue<double>("post process.local activation time.time end");

        post_process.ComputeLocalActivationTimes(react_diff->Vout(), react_diff->Dt()*react_diff->OutputSteps(), time_start, time_end);

        stream << "        - local activation time computation \n";
    } // End of Apply local activation time generation.


    // Compute action potential duration.
    if (parser.HasAttribute("post process.action potential duration")) {

        // Get the repolarization percentage.
        int apd_lvl = parser.GetValue<int>("post process.action potential duration.repolarization percentage");

        // Get the APD calculation time range.
        double time_start = parser.GetValue<double>("post process.action potential duration.time start");
        double time_end = parser.GetValue<double>("post process.action potential duration.time end");

        post_process.ComputeAPD(react_diff->Vout(), react_diff->Dt()*react_diff->OutputSteps(), time_start, time_end, apd_lvl);

        stream << "        - action potential duration computation \n";
    } // End of Apply local activation time generation.
}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_POST_PROCESS_TPP_