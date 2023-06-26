/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_TPP_

#include "config_approximation.hpp"

namespace APP_ELECTRA
{

template <short DIM, short CELL_NODES>
ConfigApproximation<DIM, CELL_NODES>::ConfigApproximation()
{}


template <short DIM, short CELL_NODES>
ConfigApproximation<DIM, CELL_NODES>::~ConfigApproximation()
{}


template <short DIM, short CELL_NODES>
void ConfigApproximation<DIM, CELL_NODES>::SetFemApproximation(std::ostream &stream) const
{
    stream << ELECTRA::Logger::Message("Numerical approximation: finite element method\n");
    stream << ELECTRA::Logger::Message("Approximation type: linear\n");
}


template <short DIM, short CELL_NODES>
void ConfigApproximation<DIM, CELL_NODES>::SetFpmApproximation(const Parser &parser, const IMP::Voronoi<DIM> &voro, CLOUDEA::Fpm<DIM> &fpm_approx, std::ostream &stream) const
{
    if(!parser.HasAttribute("numerical approximation.fpm")) {
        std::string error_msg = "Could not configure FPM numerical approximation. Provide [numerical approximation.fpm] section in the configuration .json file.";
        throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
    }

    stream << ELECTRA::Logger::Message("Numerical approximation: fragile points method\n");

    // Set penalty of the fpm approximation.
    fpm_approx.SetPenalty(parser.GetValue<double>("numerical approximation.fpm.penalty"));

    // Compute fpm approximants.
    stream << ELECTRA::Logger::Message("Computing FPM numerical approximation...\n");
    fpm_approx.Support().IdentifyInfluenceNodesVoronoi(voro);
    fpm_approx.Compute(voro);
}


template <short DIM, short CELL_NODES>
void ConfigApproximation<DIM, CELL_NODES>::SetMcmApproximation(const Parser &parser, const IMP::Grid<DIM, CELL_NODES> &grid, std::unique_ptr<CLOUDEA::Mfree<DIM>> &mcm_approx,
                                                               IMP::NodeSet &neumann_nset, std::ostream &stream) const
{
    if(!parser.HasAttribute("numerical approximation.mcm")) {
        std::string error_msg = "Could not configure MCM numerical approximation. Provide [numerical approximation.mcm] section in the configuration .json file.";
        throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
    }

    stream << ELECTRA::Logger::Message("Numerical approximation: mixed collocation method\n");

    // Set map of entries of mcm approximant types.
    std::unordered_map<std::string, CLOUDEA::MfreeType> mfree_map;
    mfree_map["mls"] = CLOUDEA::MfreeType::mls;
    mfree_map["rpi"] = CLOUDEA::MfreeType::rpi;
    mfree_map["mki"] = CLOUDEA::MfreeType::mki;

    // MCM approximant type.
    std::string mfree_type = parser.GetValue<std::string>("numerical approximation.mcm.type");
    std::transform(std::begin(mfree_type), std::end(mfree_type), std::begin(mfree_type), ::tolower);

    // Set up the MCM approximant.
    mcm_approx = CLOUDEA::MfreeFactory<DIM>::Create(mfree_map.at(mfree_type));
    mcm_approx->EditSupport().SetFieldNodesNum(grid.NodesNum());
    if (mcm_approx->Type() == CLOUDEA::MfreeType::mls) { stream << ELECTRA::Logger::Message("Approximation type: MLS - Moving Least Squares\n"); }
    else if (mcm_approx->Type() == CLOUDEA::MfreeType::rpi) { stream << ELECTRA::Logger::Message("Approximation type: RPI - Radial Point Interpolation\n"); }
    else if (mcm_approx->Type() == CLOUDEA::MfreeType::mki) { stream << ELECTRA::Logger::Message("Approximation type: MKI - Moving Kriging Interpolation\n"); }

    // Set up nodeset of neumann boundary condition nodes.
    neumann_nset.Clear();
    std::string neumann_nset_name = parser.GetValue<std::string>("numerical approximation.mcm.neumann nodeset");
    for (const auto &nset : grid.NodeSets()) {
        if (nset.second.Name() == neumann_nset_name) { neumann_nset.Set(nset.second.Name(), nset.second.NodeIds()); }
    }
    if (neumann_nset.NodeIds().size() == 0) {
        throw std::runtime_error(ELECTRA::Logger::Error("Could not extract neumann boundary node set. Check given name."));
    }
    stream << ELECTRA::Logger::Message("Neumann node set name: " + neumann_nset.Name() + "\n");


    // Optional set support domain dilatation coefficient.
    if (parser.HasAttribute("numerical approximation.mcm.support dilatation")) {
        double dilate_coeff = parser.GetValue<double>("numerical approximation.mcm.support dilatation");
        std::vector<double> coeffs(grid.NodesNum(), dilate_coeff);

        if (parser.HasAttribute("numerical approximation.mcm.support dilatation surface")) {
            double dilate_coeff_surf = parser.GetValue<double>("tissue.meshfree approximant.support dilatation surface");
            for (const auto &id : neumann_nset.NodeIds())  coeffs[id] = dilate_coeff_surf;
        }

        mcm_approx->SetSupportDilatation(coeffs);
    }

    // Get support domain spacing type.
    std::string support_spacing_type = parser.GetValue<std::string>("numerical approximation.mcm.support spacing");
    std::transform(std::begin(support_spacing_type), std::end(support_spacing_type), std::begin(support_spacing_type), ::tolower);

    // Compute support domains.
    if (support_spacing_type == "regular") {
        mcm_approx->EditSupport().ComputeRadiusFromRegularGrid(grid.Nodes());
    }
    else if (support_spacing_type == "irregular") {
        mcm_approx->EditSupport().ComputeRadiusFromIrregularGrid(grid);
    }
    else if (support_spacing_type == "immersed") {
        mcm_approx->EditSupport().ComputeRadiusFromImmersedGrid(grid);
    } else {
        std::string error_msg = "Could not configure MCM approximation. "
                                "Supported [support spacing] attribute is \"regular\", \"irregular\", or \"immersed\".";
        throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
    }

    // Get influence nodes number if available.
    int influence_nodes_num = 0;
    if (parser.HasAttribute("numerical approximation.mcm.influence nodes number"))
        influence_nodes_num = parser.GetValue<int>("numerical approximation.mcm.influence nodes number");

    // Shape parameter coefficient.
    double theta = 1.0; ;
    if (parser.HasAttribute("numerical approximation.mcm.shape parameter"))
        theta = parser.GetValue<double>("numerical approximation.mcm.shape parameter");

    // Exponent coefficient.
    double beta = 1.0;
    if (parser.HasAttribute("numerical approximation.mcm.exponent"))
        beta = parser.GetValue<double>("numerical approximation.mcm.exponent");

    // Identify influence nodes in support domains.
    stream << ELECTRA::Logger::Message("Identifying influence domain neighbor nodes...\n");
    if (influence_nodes_num > 0) {
        mcm_approx->EditSupport().IdentifyNearestInfluenceNodes(grid.Nodes(), neumann_nset, influence_nodes_num);
    } else {
        mcm_approx->EditSupport().IdentifyInfluenceNodesInRange(grid.Nodes(), grid.Nodes());
    }
    stream << ELECTRA::Logger::Message("Influence domain nodes number: ") << mcm_approx->Support().MinInfluenceNodesNum() << " - " << mcm_approx->Support().MaxInfluenceNodesNum() << "\n";

    // Compute average support domain size.
    double avg_rad = std::accumulate(std::begin(mcm_approx->Support().Radius()),
                                     std::end(mcm_approx->Support().Radius()), 0.0) / mcm_approx->Support().Radius().size();
    stream << ELECTRA::Logger::Message("Influence domain average radius: ") << avg_rad << " " << parser.GetValue<std::string>("tissue.geometry.unit") << "\n";

    // Compute average number of influence nodes for boundary and internal nodes.
    std::vector<int> bound_flag(grid.NodesNum());
    for (const auto &id : neumann_nset.NodeIds())  bound_flag[id] = 1;

    int mean_neighs_surf = 0; int cnt_surf = 0;
    int mean_neighs_in = 0; int cnt_in = 0;
    for (const auto &neighs : mcm_approx->Support().InfluenceNodeIds()) {
        auto i = &neighs - &mcm_approx->Support().InfluenceNodeIds()[0];

        if (bound_flag[i] == 1) { mean_neighs_surf += neighs.size(); cnt_surf++;}
        else { mean_neighs_in += neighs.size(); cnt_in++;}
    }
    mean_neighs_surf /= cnt_surf;
    mean_neighs_in /= cnt_in;

    stream << ELECTRA::Logger::Message("Surface influence nodes average number: ") << mean_neighs_surf << "\n";
    stream << ELECTRA::Logger::Message("Internal influence nodes average number: ") << mean_neighs_in << "\n";

    // Weight function type.
    CLOUDEA::WeightType weight_type = CLOUDEA::WeightType::unknown;
    if (parser.HasAttribute("numerical approximation.mcm.weight function")) {
        // Create map of rbf types.
        std::unordered_map<std::string, CLOUDEA::WeightType> weight_map;
        weight_map["cubic"] = CLOUDEA::WeightType::cubic;
        weight_map["quartic"] = CLOUDEA::WeightType::quartic;
        weight_map["gaussian"] = CLOUDEA::WeightType::gaussian;
        weight_map["multiquadric"] = CLOUDEA::WeightType::multiquadric;
        weight_map["polyharmonic"] = CLOUDEA::WeightType::polyharmonic;

        // Set rbf type to the mfree parameters.
        std::string temp = parser.GetValue<std::string>("numerical approximation.mcm.weight function");
        std::transform(std::begin(temp), std::end(temp), std::begin(temp), ::tolower);
        weight_type = weight_map.at(temp);
    }

    // Set tissue.meshfree approximant weight function.
    mcm_approx->SetWeightFunction(weight_type, theta, beta);

    // Monomial basis type.
    CLOUDEA::MonomialType mbasis_type = CLOUDEA::MonomialType::unknown;
    if (parser.HasAttribute("numerical approximation.mcm.monomial basis")) {
        // Create map of monomial basistypes.
        std::unordered_map<std::string, CLOUDEA::MonomialType> mbasis_map;
        mbasis_map["linear"] = CLOUDEA::MonomialType::linear;
        mbasis_map["quadratic"] = CLOUDEA::MonomialType::quadratic;
        mbasis_map["cubic"] = CLOUDEA::MonomialType::cubic;

        // Set monomial basis type to the mfree parameters.
        std::string temp = parser.GetValue<std::string>("numerical approximation.mcm.monomial basis");
        std::transform(std::begin(temp), std::end(temp), std::begin(temp), ::tolower);
        mbasis_type = mbasis_map.at(temp);
    }

    // Set tissue.meshfree approximant monomial basis.
    mcm_approx->SetMonomialBasis(mbasis_type);

    stream << ELECTRA::Logger::Message("Computing basis functions and gradients...");
    mcm_approx->Compute(grid.Nodes(), grid.Nodes());
    stream << " OK\n";
}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_