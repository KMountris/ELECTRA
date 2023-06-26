/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_CONDUCTION_SYSTEM_CONDUCTION_SYSTEM_TPP_
#define ELECTRA_CONDUCTION_SYSTEM_CONDUCTION_SYSTEM_TPP_

#include "ELECTRA/engine/conduction_system/conduction_system.hpp"


namespace ELECTRA {


template<short DIM>
ConductionSystem<DIM>::ConductionSystem() : tree_(), pmj_(), diffusion_coeffs_(), pmj_diffusion_coeffs_(),
        purkinje_node_ids_(), terminal_node_ids_(), his_node_ids_(), boltz_coeff_(0.), av_node_id_(-1)

{}


template<short DIM>
ConductionSystem<DIM>::~ConductionSystem()
{}


template<short DIM>
void ConductionSystem<DIM>::LoadTree(const std::string& filename)
{
    // Load the cardiac conduction system's tree.
    this->tree_.LoadFrom(filename);
}


template<short DIM>
void ConductionSystem<DIM>::SetDiffusionCoeffs(double diffusion_coeff)
{
    if (this->tree_.NodesNum() == 0) {
        auto error_str = "Could not set diffusion coefficients. First, load the cardiac conduction system tree.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    this->diffusion_coeffs_.clear();
    this->diffusion_coeffs_.resize(this->tree_.NodesNum(), diffusion_coeff);
}


template<short DIM>
void ConductionSystem<DIM>::SetDiffusionCoeffs(const std::vector<double> &diffusion_coeffs)
{
    if (this->tree_.NodesNum() == 0) {
        auto error_str = "Could not set diffusion coefficients. First, load the cardiac conduction system tree.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    if (this->tree_.NodesNum() != static_cast<int>(diffusion_coeffs.size()) ) {
        auto error_str = "Could not set diffusion coefficients. Diffusion coefficients vector is not equal to conduction sytem's nodes number.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    this->diffusion_coeffs_ = diffusion_coeffs;
}


template<short DIM>
void ConductionSystem<DIM>::SetPmjDiffusionCoeffs(double diffusion_coeff)
{
    if (this->pmj_.size() == 0) {
        auto error_str = "Could not set diffusion coefficients for the PMJ. First, compute the PMJ of the cardiac conduction system.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    this->pmj_diffusion_coeffs_.clear();
    this->pmj_diffusion_coeffs_.resize(this->pmj_.size(), diffusion_coeff);
}


template<short DIM>
void ConductionSystem<DIM>::SetPurkinjeNodeIds(const std::string &nset_name)
{
    // Find the node set in the conduction system's tree.
    bool found = false;
    for (const auto &nset : this->tree_.NodeSets()) {
        if (nset.Name() == nset_name) { this->purkinje_node_ids_ = nset.NodeIds(); found = true; break; }
    }

    if (!found) {
        auto error_str = "Could not set Purkinje node ids. The node set was not found in the conduction system's tree.";
        throw std::invalid_argument(Logger::Error(error_str));
    }
}


template<short DIM>
void ConductionSystem<DIM>::SetTerminalNodeIds(const std::string &nset_name)
{
    // Find the node set in the conduction system's tree.
    bool found = false;
    for (const auto &nset : this->tree_.NodeSets()) {
        if (nset.second.Name() == nset_name) { this->terminal_node_ids_ = nset.second.NodeIds(); found = true; break; }
    }

    if (!found) {
        auto error_str = "Could not set terminal node ids. The node set was not found in the conduction system's tree.";
        throw std::invalid_argument(Logger::Error(error_str));
    }
}


template<short DIM>
void ConductionSystem<DIM>::SetHisNodeIds(const std::string &nset_name)
{
    // Find the node set in the conduction system's tree.
    bool found = false;
    for (const auto &nset : this->tree_.NodeSets()) {
        if (nset.Name() == nset_name) { this->his_node_ids_ = nset.NodeIds(); found = true; break; }
    }

    if (!found) {
        auto error_str = "Could not set His-bundle node ids. The node set was not found in the conduction system's tree.";
        throw std::invalid_argument(Logger::Error(error_str));
    }
}


template<short DIM>
void ConductionSystem<DIM>::SetAvNodeId(const std::string &nset_name)
{
    // Find the node set in the conduction system's tree.
    bool found = false;
    for (const auto &nset : this->tree_.NodeSets()) {
        if (nset.second.Name() == nset_name) {
            if (nset.second.NodeIds().size() > 1) {
                auto error_str = "Could not set atrioventricular node id. The node set contains more than one node indices.";
                throw std::invalid_argument(Logger::Error(error_str));
            }
            this->av_node_id_ = nset.second.NodeIds()[0];
            found = true;
            break;
        }
    }

    if (!found) {
        auto error_str = "Could not set atrioventricular node id. The node set was not found in the conduction system's tree.";
        throw std::invalid_argument(Logger::Error(error_str));
    }
}


template<short DIM>
void ConductionSystem<DIM>::ComputePmjs(const std::vector<IMP::Vec<DIM, double>> &tissue_nodes, double pmj_radius)
{
    // Boost definitions.
    const short dimss = DIM;
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;
    typedef bg::model::point<double, dimss, bg::cs::cartesian> point;
    typedef std::pair<point, std::size_t> value;

    // Check that terminal node indices are available.
    if (this->terminal_node_ids_.size() <= 0) {
        auto err_msg = "Could not compute PMJs. Set terminal node ids first.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Reset the pmj container.
    this->pmj_.clear();
    this->pmj_.resize(this->terminal_node_ids_.size(), std::vector<int>());

    // Populate rtree with tissue nodes.
    auto rtree = bgi::rtree<value, bgi::quadratic<16>>{};
    auto p = point{};
    for (std::size_t i = 0; i != tissue_nodes.size(); ++i) {
        if constexpr (DIM == 1) { 
            bg::set<0>(p, tissue_nodes[i][0]);
        } else if constexpr (DIM == 2) { 
            bg::set<0>(p, tissue_nodes[i][0]);
            bg::set<1>(p, tissue_nodes[i][1]);
        } else if constexpr (DIM == 3) { 
            bg::set<0>(p, tissue_nodes[i][0]);
            bg::set<1>(p, tissue_nodes[i][1]);
            bg::set<2>(p, tissue_nodes[i][2]);
        }
         
        rtree.insert(std::make_pair(p,i));
    }

    // Search for the pmj in the tissue nodes for each terminal node of the conduction system
    auto connected_nodes = std::vector<value>{};
    auto term_nid = int{0};
    auto terminal_point = point{};
    for (std::size_t i = 0; i != this->pmj_.size(); ++i) {

        // Get the index of the terminal node.
        term_nid = this->terminal_node_ids_[i];

        // Coordinates of the current terminal point.
        if constexpr (DIM == 1) { 
            bg::set<0>(terminal_point, this->tree_.Nodes(term_nid)[0]);
        } else if constexpr (DIM == 2) { 
            bg::set<0>(terminal_point, this->tree_.Nodes(term_nid)[0]);
            bg::set<1>(terminal_point, this->tree_.Nodes(term_nid)[1]);
        } else if constexpr (DIM == 3) { 
            bg::set<0>(terminal_point, this->tree_.Nodes(term_nid)[0]);
            bg::set<1>(terminal_point, this->tree_.Nodes(term_nid)[1]);
            bg::set<2>(terminal_point, this->tree_.Nodes(term_nid)[2]);
        }
        
        // Perform the search.
        rtree.query(bgi::satisfies(
            [&](value const& v) {return bg::distance(v.first, terminal_point) <= pmj_radius;}), std::back_inserter(connected_nodes)
        );

        // Store the connected tissue nodes to the terminal point in the pmj container.
        for (const auto &cnode : connected_nodes) {
            this->pmj_[i].emplace_back(cnode.second);
        }

        // Clear connected nodes container for the next iteration.
        connected_nodes.clear();
    }

}


template<short DIM>
void ConductionSystem<DIM>::ComputeDiffuseTransition(double kappa)
{
    if (this->diffusion_coeffs_.size() == 0 || this->pmj_diffusion_coeffs_.size() == 0) {
        auto err_msg = "Could not compute diffusivity transition. First, set diffusivity for conduction system and Purkinje-myocardial junctions.";
        throw std::runtime_error(Logger::Error(err_msg));
    }

    // Create a mapping from nodes to connected segments.
    auto node_to_seg_mapping = std::vector<std::vector<int>>(this->tree_.NodesNum());
    for (const auto &seg : this->tree_.Cells()) {
        auto seg_id = &seg - &this->tree_.Cells()[0];

        // Map segment id to its nodes.
        node_to_seg_mapping[seg.N(0)].emplace_back(seg_id);
        node_to_seg_mapping[seg.N(1)].emplace_back(seg_id);
    }

    // Get the diffusivity value at the Purkinje-myocardial junction.
    auto dpmj = *std::min_element(std::begin(this->pmj_diffusion_coeffs_), std::end(this->pmj_diffusion_coeffs_));
    auto dp = double{0.};

    // Adjust the diffusivity of the transition terminals nodes
    auto trans_node_ids = std::vector<int>{};
    trans_node_ids.reserve(10);
    auto center = IMP::Vec<DIM, double>{};
    auto conn_nid = int{0};
    for (const auto &term_nid : this->terminal_node_ids_) {

        // Get first connected node to terminal.
        conn_nid = this->tree_.Cells(node_to_seg_mapping[term_nid][0]).N(0);
        trans_node_ids.emplace_back(conn_nid);

        // Get the rest transition nodes. Max number is 10.
        for (int iter = 0; iter != 9; ++iter) {
            // Stop the transition nodes if a bifurcation is found.
            if (node_to_seg_mapping[conn_nid].size() > 2) {
                break;
            }

            // Otherwise iterate segments of the current transition node.
            for (const auto &seg_id : node_to_seg_mapping[conn_nid]) {
                // Search the segment ending at the current transition node.
                if (this->tree_.Cells(seg_id).N(1) == conn_nid) {
                    // Add the starting node of the segment in the transition nodes.
                    conn_nid = this->tree_.Cells(seg_id).N(0);
                    trans_node_ids.emplace_back(conn_nid);
                    break;
                }
            }
        }

        // Compute coordinates of the transition area's center.
        if (trans_node_ids.size() == 2) {
            center = this->tree_.Nodes(trans_node_ids[0]) + this->tree_.Nodes(trans_node_ids[1]);
            center *= 0.5;
        } else {
            auto tcntr_id = std::ceil(0.5*trans_node_ids.size())-1.;
            center = this->tree_.Nodes(static_cast<int>(tcntr_id));
        }

        // Apply transition diffusivity coefficient
        for (const auto &id : trans_node_ids) {
            double x = std::sqrt(this->tree_.Nodes(id).Distance2(this->tree_.Nodes(term_nid)));
            double xc = std::sqrt(center.Distance2(this->tree_.Nodes(term_nid)));

            // Update the diffusivity of the node.
            dp = this->diffusion_coeffs_[id];
            this->diffusion_coeffs_[id] = dpmj + ( (dp-dpmj) / (1. + std::exp(kappa*(x-xc))) );
        }

        // Clear transition nodes container for next iteration.
        trans_node_ids.clear();

    }
}


template<short DIM>
void ConductionSystem<DIM>::SavePmjs(const std::string &filename)
{
    std::ofstream out(filename, std::ios::out);

    for (const auto &pmj : this->pmj_) {
        for (const auto &nid : pmj) {
            out << nid << " ";
        }
        out << "\n\n";
    }

    out.close();

}



} //end of namespace ELECTRA


#endif //ELECTRA_CONDUCTION_SYSTEM_CONDUCTION_SYSTEM_TPP_