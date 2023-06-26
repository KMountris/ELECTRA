/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file conduction_system.hpp
   \brief ConductionSystem class header file.
   \author Konstantinos A. Mountris
   \date 18/02/2020
*/

#ifndef ELECTRA_CONDUCTION_SYSTEM_CONDUCTION_SYSTEM_HPP_
#define ELECTRA_CONDUCTION_SYSTEM_CONDUCTION_SYSTEM_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"

#include <IMP/IMP>
#include <Eigen/Dense>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <vector>
#include <string>
#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <fstream>

namespace ELECTRA {

/** \addtogroup Conduction-System \{ */

/**
 * \class ConductionSystem
 * \brief Class implemmenting the cardiac conduction system.
 * \tparam DIM The spatial dimensions of the cardiac conduction system.
 * \author Konstantinos A. Mountris
 */
template<short DIM>
class ConductionSystem
{
private:

    IMP::Mesh<DIM, 2> tree_;                       /**< The mesh of the conduction system tree */

    std::vector<std::vector<int>> pmj_;            /**< The Purkinje-myocardium junctions for each terminal node of the conduction system */

    std::vector<double> diffusion_coeffs_;         /**< The diffusion coefficient value for each node of the conduction system */

    std::vector<double> pmj_diffusion_coeffs_;     /**< The diffusion coefficient value for each Purkinje-myocardium junction */

    std::vector<int> purkinje_node_ids_;           /**< The indices of the conduction system's Purkinje nodes */

    std::vector<int> terminal_node_ids_;           /**< The indices of the conduction system's terminal nodes */

    std::vector<int> his_node_ids_;                /**< The indices of the conduction system's His-bundle nodes */

    double boltz_coeff_;                           /**< The coefficient controlling the Boltzmann curve's slope */

    int av_node_id_;                               /**< The index of the conduction system's atrio-ventricular node */


public:

    /**
     * \brief The default constructor of the ConductionSystem class.
     */
    ConductionSystem();


    /**
     * \brief The default destructor of the ConductionSystem class.
     */
    virtual ~ConductionSystem();


    /**
     * \brief Load the cardiac conduction system's tree mesh from a file.
     * File expected in Abaqus format (.inp).
     * \param [in] filename The filename of the mesh to be loaded.
     * \return [void]
     */
    void LoadTree(const std::string& filename);


    /**
     * \brief Set the same diffusion coefficients for all the conduction system nodes.
     * \param [in] diffusion_coeff The diffusion coefficient.
     * \return [void]
     */
    void SetDiffusionCoeffs(double diffusion_coeff);


    /**
     * \brief Set the diffusion coefficients to the conduction system nodes from a vector container.
     * \param [in] diffusion_coeffs The vector container with the diffusion coefficients for each conduction system node.
     * \return [void]
     */
    void SetDiffusionCoeffs(const std::vector<double> &diffusion_coeffs);


    /**
     * \brief Set the same diffusion coefficients for all the Purkinje-myocardial junctions.
     * \param [in] diffusion_coeff The diffusion coefficient.
     * \return [void]
     */
    void SetPmjDiffusionCoeffs(double diffusion_coeff);


    /**
     * \brief Set the indices that correspond to the Purkinje nodes of the cardiac conduction system.
     * \param [in] nset_name The name of the nodeset of the Purkinje nodes in the tree of the cardiac conduction system.
     * \return [void]
     */
    void SetPurkinjeNodeIds(const std::string &nset_name);


    /**
     * \brief Set the indices that correspond to the terminal nodes of the cardiac conduction system.
     * \param [in] nset_name The name of the nodeset of the terminal nodes in the tree of the cardiac conduction system.
     * \return [void]
     */
    void SetTerminalNodeIds(const std::string &nset_name);


    /**
     * \brief Set the indices that correspond to the His-bundle nodes of the cardiac conduction system.
     * \param [in] nset_name The name of the nodeset of the His-bundle nodes in the tree of the cardiac conduction system.
     * \return [void]
     */
    void SetHisNodeIds(const std::string &nset_name);


    /**
     * \brief Set the indices that correspond to the atrio-ventricular node of the cardiac conduction system.
     * \param [in] nset_name The name of the nodeset of the atrio-ventricular node in the tree of the cardiac conduction system.
     * \return [void]
     */
    void SetAvNodeId(const std::string &nset_name);


    /**
     * \brief Set the Boltzmann Coeff object
     * \param [in] coeff The value of the coefficient of the Boltzmann curve's slope.
     * \return [void]  
     */
    inline void SetBoltzmannCoeff(double coeff) { this->boltz_coeff_ = coeff; }


    /**
     * \brief Compute the Purkinje-myocardial junctions for each terminal node of the conduction system.
     * \param [in] tissue_nodes The nodes of the tissue geometry to be tested if belong in a Purkinje-myocardial junction.
     * \param [in] pmj_radius The radius of a ND sphere for Purkinje-myocardial junctions search.
     * \return [void]
     */
    void ComputePmjs(const std::vector<IMP::Vec<DIM, double>> &tissue_nodes, double pmj_radius);


    /**
     * \brief Computes the transition of diffusivity from the conduction system towards the Purkinje-myocardial junctions.
     * 
     * The transition is performed at the last 10 segments of the conduction system. It is modelled using the sigmoid curve:
     * \f[ d_i = \frac{d_p-d_0}{1+e^{k(x_c-x_i)}} + d_0 \f]
     * where \f$d_i\f$ is the diffusivity of the transition node \f$i\f$, \f$d_p\f$ is the diffusivity of the conduction system tree, 
     * \f$d_0\f$ is the diffusivity of the Purkinje-myocardial junction, \f$k\f$ is the coefficient that controls the slope of the curve,
     * \f$x_c\f$ is the distance of the transition area center from the terminal node, and \f$x_i\f$ is the distance of the transition 
     * node \f$i\f$ from the terminal node.
     * 
     * \param [in] kappa The slope coefficient of the sigmoid curve. 
     */
    void ComputeDiffuseTransition(double kappa);


    /**
     * \brief Save the myocardial node indices for each Purkinje-myocardial junction in a file.
     * \param [in] filename The name of the file to save the data.
     * \return [void]
     */
    void SavePmjs(const std::string &filename);


    /**
     * \brief Get the diffusion coefficients of the conduction system's nodes. 
     * \return [const std::vector<double>&] The diffusion coefficients of the conduction system's nodes.
     */
    inline const std::vector<double> & DiffusionCoeffs() const { return this->diffusion_coeffs_; }


    /**
     * \brief Get the diffusion coefficient of a specific node of the conduction system.
     * Fast access with no range check. 
     * \return [double] The diffusion coefficient of a specific node of the conduction system.
     */
    inline double DiffusionCoeffs(std::size_t id) const { return this->diffusion_coeffs_[id]; }


    /**
     * \brief Get the diffusion coefficient of a specific node of the conduction system.
     * Slower access with range check. 
     * \return [double] The diffusion coefficient of a specific node of the conduction system.
     */
    inline double DiffusionCoeffsAt(std::size_t id) const { return this->diffusion_coeffs_.at(id); }


        /**
     * \brief Get the Purkinje-myocardial junction diffusion coefficients of the conduction system's nodes. 
     * \return [const std::vector<double>&] The Purkinje-myocardial junction diffusion coefficients of the conduction system's nodes.
     */
    inline const std::vector<double> & PmjDiffusionCoeffs() const { return this->pmj_diffusion_coeffs_; }


    /**
     * \brief Get the Purkinje-myocardial junction diffusion coefficient of a specific node of the conduction system.
     * Fast access with no range check. 
     * \return [double] The Purkinje-myocardial junction diffusion coefficient of a specific node of the conduction system.
     */
    inline double PmjDiffusionCoeffs(std::size_t id) const { return this->pmj_diffusion_coeffs_[id]; }


    /**
     * \brief Get the Purkinje-myocardial junction diffusion coefficient of a specific node of the conduction system.
     * Slower access with range check. 
     * \return [double] The Purkinje-myocardial junction diffusion coefficient of a specific node of the conduction system.
     */
    inline double PmjDiffusionCoeffsAt(std::size_t id) const { return this->pmj_diffusion_coeffs_.at(id); }


    /**
     * \brief Get the Purkinje-myocardial junctions of the conduction system. 
     * \return [const std::vector<std::vector<int>>&] The Purkinje-myocardial junctions of the conduction system.
     */
    inline const std::vector<std::vector<int>> & Pmj() const { return this->pmj_; } 


    /**
     * \brief Get the terminal node indices of the conduction system.
     * \return [const std::vector<int>&] The terminal node indices of the conduction system. 
     */
    inline const std::vector<int> & TerminalNodeIds() const { return this->terminal_node_ids_; }


    /**
     * \brief Get the nodes of the conduction system.
     * \return [const std::vector<IMP::Vec<DIM, double>>&] The nodes of the conduction system.
     */
    inline const std::vector<IMP::Vec<DIM, double>> & Nodes() const { return this->tree_.Nodes(); }


    /**
     * \brief Get a specific node of the conduction system.
     * Fast access without range check.
     * \return [const IMP::Vec<DIM, double>&] The requested node of the conduction system.
     */
    inline const IMP::Vec<DIM, double> & Nodes(std::size_t id) const { return this->tree_.Nodes(id); }


    /**
     * \brief Get a specific node of the conduction system.
     * Slower access with range check.
     * \return [const IMP::Vec<DIM, double>&] The requested node of the conduction system.
     */
    inline const IMP::Vec<DIM, double> & NodesAt(std::size_t id) const { return this->tree_.NodesAt(id); }


    /**
     * \brief Get the segments of the conduction system.
     * \return [const std::vector<IMP::Cell<DIM, 2>>&] The segments of the conduction system.
     */
    inline const std::vector<IMP::Cell<DIM, 2>> & Segments() const { return this->tree_.Cells(); }


    /**
     * \brief Get the nodesets of the conduction system.
     * \return [const std::unordered_map<std::string,IMP::NodeSet>&] The nodesets of the conduction system.
     */
    inline const std::unordered_map<std::string,IMP::NodeSet> & NodeSets() const { return this->tree_.NodeSets(); }


    /**
     * \brief Get the node set of the conduction system with the specified name.
     * \param [in] name The name of the node set.
     * \return [const IMP::NodeSet&] The node set of the conduction system with the specified name.
     */
    inline const IMP::NodeSet & NodeSets(const std::string &name) const { return this->tree_.NodeSets(name); }


    /**
     * \brief Get the number of the nodes in the conduction system tree. 
     * \return [int] The number of the nodes in the conduction system tree.
     */
    inline int NodesNum() const { return this->tree_.NodesNum(); }


    /**
     * \brief Get the number of the segments in the conduction system tree. 
     * \return [int] The number of the segments in the conduction system tree.
     */
    inline int SegmentsNum() const { return this->tree_.CellsNum(); }

};


/** \} End of Doxygen Groups */

} //end of namespace ELECTRA

#endif //ELECTRA_CONDUCTION_SYSTEM_CONDUCTION_SYSTEM_HPP_


#include "ELECTRA/engine/conduction_system/conduction_system.tpp"