/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
   \file ldrbm.hpp
   \brief Ldrbm class header file.
   \author Konstantinos A. Mountris
   \date 15/03/2021
*/

#ifndef ELECTRA_FIBERS_LDRBM_HPP_
#define ELECTRA_FIBERS_LDRBM_HPP_

#include "ELECTRA/engine/fibers/atri_fiber_rules.hpp"
#include "ELECTRA/engine/fibers/atri_tags.hpp"
#include "ELECTRA/engine/fibers/ventri_fiber_rules.hpp"
#include "ELECTRA/engine/fibers/ventri_tags.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <CLOUDEA/CLOUDEA>
#include <IMP/IMP>
#include <Eigen/Eigen>

#include <string>
#include <vector>
#include <list>
#include <thread>
#include <stdexcept>
#include <exception>
#include <memory>
#include <limits>

namespace ELECTRA
{

/** \addtogroup Fibers \{ */

/**
 * \class Ldrbm
 * \brief Class implemmenting the Laplace-Dirichlet rule based method for estimation of cardiac fibers' direction according to \cite bayer2012 \cite doste2019 \cite piersanti2021.
 * \tparam DIM the spatial dimension of the Ldrbm model.
 * \tparam CELL_NODES the number of nodes in each cell of the Ldrbm model.
 */
template<short DIM, short CELL_NODES>
class Ldrbm
{
private:

    Eigen::SparseMatrix<double, Eigen::RowMajor> laplacian_mat_;        // The Laplacian operator matrix

    Eigen::MatrixXd transmural_direction_;                              // The transmural direction of each fiber

    Eigen::MatrixXd appendage_veins_direction_;                         // The appendage to veins direction of each atrium fiber

    Eigen::MatrixXd inter_veins_direction_;                             // The direction between the veins of each atrium fiber

    Eigen::MatrixXd valve_veins_direction_;                             // The direction between the valve and the veins of each atrium fiber

    Eigen::MatrixXd atri_tricuspid_direction_;                          // The direction between the septal and free part of the tricuspid valve of each atrium fiber

    Eigen::MatrixXd apicobasal_direction_;                              // The apicobasal direction of each ventricle fiber

    Eigen::MatrixXd fiber_direction_;                                   // The direction of the fibers

    Eigen::MatrixXd sheet_direction_;                                   // The direction of the fiber sheets

    Eigen::VectorXd transmural_distance_;                               // The transmural distance of each fiber

    Eigen::VectorXd appendage_veins_distance_;                          // The appendage to veins distance of each atrium fiber

    Eigen::VectorXd inter_veins_distance_;                              // The distance between the veins of each atrium fiber

    Eigen::VectorXd valve_veins_distance_;                              // The distance between the valve and the veins of each atrium fiber

    Eigen::VectorXd atri_tricuspid_distance_;                           // The direction between the septal and free part of the tricuspid valve of each atrium fiber

    Eigen::VectorXd septal_distance_;                                   // The septal distance of each ventricle fiber

    Eigen::VectorXd apicobasal_distance_;                               // The apicobasal distance of each ventricle fiber

    Eigen::VectorXd intraventricular_function_;                         // The intraventricular function of each ventricle fiber

    std::vector<int> la_node_ids_;                                      // The indices of nodes in the left atrium

    std::vector<int> ra_node_ids_;                                      // The indices of nodes in the right atrium

    std::vector<int> lv_node_ids_;                                      // The indices of nodes in the left ventricle

    std::vector<int> rv_node_ids_;                                      // The indices of nodes in the right ventricle

    std::vector<int> ventricle_interface_node_ids_;                     // The indices of nodes at the interface of the two ventricles

    std::vector<int> septum_node_ids_;                                  // The indices of nodes in the septum


protected:

    /**
     * \brief Slice a vector by extracting specific entries.
     * \param [in] eigen_vec The original vector from where entries are extracted.
     * \param [in] row_ids The indices of rows that correspond to the desired entries.
     * \return [Eigen::VectorXd] The sliced subvector containing only the entries that are specified by row_ids.
     */
    Eigen::VectorXd SliceVector(const Eigen::VectorXd &eigen_vec, const std::vector<int> &row_ids);


    /**
     * \brief Slice a sparse matrix with row major ordering by extracting specific rows and columns.
     * \param [in] eigen_spmat The original sparse matrix in row major ordering to extract the desired rows and columns.
     * \param [in] row_ids The indices of the desired rows.
     * \param [in] col_ids The indices of the desired columns.
     * \return [Eigen::SparseMatrix<double, Eigen::RowMajor>] The sliced submatrix containing only the specified rows by row_ids and columns by col_ids.
     */
    Eigen::SparseMatrix<double, Eigen::RowMajor> SliceSparseMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor> &eigen_spmat,
                                                                   const std::vector<int> &row_ids, const std::vector<int> &col_ids);


    /**
     * \brief Solve the laplacian equation for given Dirichlet boundary conditions.
     * \param [in] nsets
     * \param [in] dirichlet_tags
     * \param [in] dirichlet_values
     * \param [out] solution
     * \return [void]
     */
    void SolveLaplacian(const std::unordered_map<std::string, IMP::NodeSet> &nsets, const std::vector<std::string> &dirichlet_tags,
                        const std::vector<double> &dirichlet_values, Eigen::VectorXd &solution);


    /**
     * @brief Compute the gradient of a scalar field in a given mesh.
     *
     * The gradient is computed using the node-based Green-Gauss gradient scheme at the centroid of each cell
     * and mapped to the nodes by computing the weighted sum of the gradient in the attached cells to each node.
     *
     * @param [in] mesh The mesh on which the gradient is computed.
     * @param [in] scalars The scalar field for which we compute the gradient.
     * @param [in] gradient The matrix of the scalar field's gradient.
     */
    void ComputeGradient(const IMP::Mesh<DIM,CELL_NODES> &mesh, const Eigen::VectorXd &scalars, Eigen::MatrixXd &gradient);


    /**
     * @brief Compute the gradient of a scalar field in a given mesh.
     *
     * The gradient is computed using the node-based Green-Gauss gradient scheme at the centroid of each cell
     * and mapped to the nodes by computing the weighted sum of the gradient in the attached cells to each node.
     *
     * @param [in] mesh The mesh on which the gradient is computed.
     * @param [in] scalars The scalar field for which we compute the gradient.
     * @param [in] gradient The matrix of the scalar field's gradient.
     */
    void ComputeGradient(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const Eigen::VectorXd &scalars, Eigen::MatrixXd &gradient);


public:

    /**
     * \brief Default constructor of Ldrbm
     */
    Ldrbm();


    /**
     * @brief Default destructor of Ldrbm
     *
     */
    virtual ~Ldrbm();


    /**
     * \brief Assemble the Laplacian operator using FEM.
     * \param [in] mesh The mesh discretization of the cardiac geometry.
     */
    void AssembleLaplacian(const IMP::Mesh<DIM,CELL_NODES> &mesh);


    /**
     * \brief Assemble the Laplacian operator using FPM.
     * \param [in] voro The Voronoi tesselation of the cardiac geometry.
     */
    void AssembleLaplacian(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, double penalty);


    /**
     * \brief
     * \param nsets
     */
    void ComputeAtriTransmural(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags);


    /**
     * \brief
     * \param nsets
     */
    void ComputeAtriTransmural(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriAppendageVeins(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriAppendageVeins(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriInterVeins(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriInterVeins(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriValveVeins(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriValveVeins(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriTricuspid(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeAtriTricuspid(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags);


    /**
     * \brief
     * \param nsets
     */
    void ComputeVentriTransmural(const IMP::Mesh<DIM,CELL_NODES> &mesh, const VentriTags &ventri_tags);


    /**
     * @brief Construct a new Compute Transmural object
     * @param voro
     * @param fpm
     * @param ventri_tags
     */
    void ComputeVentriTransmural(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const VentriTags &ventri_tags);


    /**
     * \brief
     * \param nsets
     */
    void ComputeVentriSeptal(const IMP::Mesh<DIM, CELL_NODES> &mesh, const VentriTags &ventri_tags, double threshold);


    /**
     * \brief
     * \param nsets
     */
    void ComputeVentriSeptal(const IMP::Voronoi<DIM> &voro, const VentriTags &ventri_tags, double threshold);


    /**
     * @brief
     * @param nsets
     */
    void ComputeVentriApicoBasal(const IMP::Mesh<DIM,CELL_NODES> &mesh, const VentriTags &ventri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeVentriApicoBasal(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const VentriTags &ventri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeVentriIntraventricular(const IMP::Mesh<DIM,CELL_NODES> &mesh, const VentriTags &ventri_tags);


    /**
     * @brief
     * @param nsets
     */
    void ComputeVentriIntraventricular(const IMP::Voronoi<DIM> &voro, const VentriTags &ventri_tags);


    /**
     * @brief
     *
     * @param mesh
     */
    void ExtractAtriumNodalPartitions(const IMP::Mesh<DIM,CELL_NODES> &mesh);


    /**
     * @brief
     *
     * @param mesh
     */
    void ExtractVentricleNodalPartitions(const IMP::Mesh<DIM,CELL_NODES> &mesh);


    /**
     * @brief
     *
     */
    void ComputeAtriFibers(const AtriFiberRules &rules);


    /**
     * @brief
     *
     */
    void ComputeVentriFibers(const VentriFiberRules &rules);


    //////
    void TestSaveAtriFibers(const IMP::Mesh<DIM,CELL_NODES> &mesh, const std::string &outname);

    void TestSaveVentriFibers(const IMP::Mesh<DIM,CELL_NODES> &mesh, const std::string &outname);


    /**
     * @brief 
     * 
     * @return const std::vector<int>& 
     */
    inline auto & LeftAtriNodeIds() const { return this->la_node_ids_; }


    /**
     * @brief 
     * 
     * @return const std::vector<int>& 
     */
    inline auto & RightAtriNodeIds() const { return this->ra_node_ids_; }


    /**
     * @brief 
     * 
     * @return const std::vector<int>& 
     */
    inline auto & LeftVentriNodeIds() const { return this->lv_node_ids_; }


    /**
     * @brief 
     * 
     * @return const std::vector<int>& 
     */
    inline auto & RightVentriNodeIds() const { return this->rv_node_ids_; }


    /**
     * @brief 
     * 
     * @return const std::vector<int>& 
     */
    inline auto & VentriInterfaceNodeIds() const { return this->ventricle_interface_node_ids_; }


    /**
     * @brief 
     * 
     * @return const std::vector<int>& 
     */
    inline auto & SeptumNodeIds() const { return this->septum_node_ids_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::MatrixXd&] 
     */
    inline auto & LongFiberDirection() const { return this->fiber_direction_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::MatrixXd&] 
     */
    inline auto & SheetFiberDirection() const { return this->sheet_direction_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & TransmuralDistance() const { return this->transmural_distance_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & AppendageVeinsDistance() const { return this->appendage_veins_distance_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & InterVeinsDistance() const { return this->inter_veins_distance_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & ValveVeinsDistance() const { return this->valve_veins_distance_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & AtriTricuspidDistance() const { return this->atri_tricuspid_distance_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & ApicobasalDistance() const { return this->apicobasal_distance_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & SeptalDistance() const { return this->septal_distance_; }


    /**
     * @brief 
     * 
     * @return [const Eigen::VectorXd&] 
     */
    inline auto & IntraventricularFunction() const { return this->intraventricular_function_; }


    /**
     * @brief 
     * 
     * @return const Eigen::MatrixXd &
     */
    inline auto & TransmuralDirection() const { return this->transmural_direction_; }


    /**
     * @brief 
     * 
     * @return const Eigen::MatrixXd &
     */
    inline auto & AppendageVeinsDirection() const { return this->appendage_veins_direction_; }


    /**
     * @brief 
     * 
     * @return const Eigen::MatrixXd &
     */
    inline auto & InterVeinsDirection() const { return this->inter_veins_direction_; }


    /**
     * @brief 
     * 
     * @return const Eigen::MatrixXd &
     */
    inline auto & ValveVeinsDirection() const { return this->valve_veins_direction_; }


    /**
     * @brief 
     * 
     * @return const Eigen::MatrixXd &
     */
    inline auto & AtriTricuspidDirection() const { return this->atri_tricuspid_direction_; }


    /**
     * @brief 
     * 
     * @return const Eigen::MatrixXd &
     */
    inline auto & ApicobasalDirection() const { return this->apicobasal_direction_; }

};


} // End of namespace ELECTRA

#endif //ELECTRA_FIBERS_LDRBM_HPP_

#include "ELECTRA/engine/fibers/ldrbm.tpp"