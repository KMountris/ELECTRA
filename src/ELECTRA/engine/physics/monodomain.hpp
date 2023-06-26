/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


/**
   \file monodomain.hpp
   \brief Monodomain class header file.
   \author Konstantinos A. Mountris
   \date 15/03/2019
*/

#ifndef ELECTRA_ENGINE_PHYSICS_MONODOMAIN_HPP_
#define ELECTRA_ENGINE_PHYSICS_MONODOMAIN_HPP_


#include "ELECTRA/engine/physics/reaction_diffusion.hpp"
#include "ELECTRA/engine/materials/electric_basic.hpp"
#include "ELECTRA/engine/electrophysiology/ep_basic.hpp"
#include "ELECTRA/engine/electrophysiology/ep_varying_params.hpp"
#include "ELECTRA/engine/utilities/algorithm.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"

#include <CLOUDEA/CLOUDEA>
#include <IMP/Tesselations>
#include <IMP/Topology>

#include <Eigen/Eigen>

#include <vector>
#include <string>
#include <memory>
#include <numeric>
#include <utility>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <thread>
#include <mutex>
#include <cmath>



namespace ELECTRA {

/** \addtogroup Physics \{ */


/**
 * \class Monodomain
 * \brief Class implemmenting the monodomain model for reaction-diffusion in cardiac electrophysiology problems.
 * The monodomain model is implemented using the Operator-Splitting method \cite qu1999.
 * \tparam DIM The dimensions of the monodomain model's domain. Supported: [1 | 2 | 3].
 * \tparam CELL_NODES The number of nodes composing the cells of the monodomain model's domain.
 */
template<short DIM, short CELL_NODES=1>
class Monodomain : public ReactionDiffusion<DIM,CELL_NODES>
{

private:

    std::vector<Eigen::VectorXd> vout_;                         /**< The nodal values of the potential at the selected output time instants */

    Eigen::SparseMatrix<double, Eigen::RowMajor> stiff_mat_;    /**< The stiffness matrix of the monodomain model's diffusion term */

    Eigen::VectorXd mass_vec_;                                  /**< The lamped mass vector of the monodomain model's diffusion term */

    CLOUDEA::ThreadLoopManager thread_loop_manager_;            /**< The managing object for the multithreaded loop execution */

    std::size_t threads_number_;                                /**< The number of threads for parallel execution of the monodomain model */

    std::mutex thread_mutex_;                                   /**< The mutex for the multithreaded execution */

protected:

    inline void ApplyFpmCorrection2D(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
            const CLOUDEA::Fpm<DIM> &fpm, double penalty);


    inline void ApplyFpmCorrection3D(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
            const CLOUDEA::Fpm<DIM> &fpm, double penalty);


    /**
     * @brief
     * @param material
     * @param stimulus
     * @param v_old
     * @param step
     */
    inline void ComputeReaction(const std::vector<Stimulus> &stimuli, int step, Eigen::VectorXd &v_nodal,
            std::vector<std::unique_ptr<EpBasic>> &cells);


    /**
     * @brief
     *
     * @param thread_id
     * @param material
     * @param stimulus
     * @param v_old
     * @param step
     */
    inline void ReactionCallbackStandard(std::size_t thread_id, const std::vector<Stimulus> &stimuli,
            int step, Eigen::VectorXd &v_nodal, std::vector<std::unique_ptr<EpBasic>> &cells);


    /**
     * @brief 
     * 
     * \bug The peak is moved forward in time depending on the dt value
     * 
     * @param material 
     * @param stimulus 
     * @param v_old 
     * @param step 
     */
    inline void ReactionCallbackAdaptive(std::size_t thread_id, const std::vector<Stimulus> &stimuli,
            int step, Eigen::VectorXd &v_nodal, std::vector<std::unique_ptr<EpBasic>> &cells);


    /**
     * @brief 
     * 
     * @param conductivity_matrix_trans 
     * @param v_old 
     * @param v_new 
     */
    inline void ComputeDiffusionExplicit(const Eigen::SparseMatrix<double, Eigen::RowMajor> &stiff_mat, const Eigen::VectorXd &mass_vec, 
            double dt, Eigen::VectorXd &v_nodal);


    /**
     * @brief 
     * 
     * @param stiffness_mat 
     * @param mass_mat 
     * @param v_old 
     * @param v_new 
     */
    inline void ComputeDiffusionImplicit(const Eigen::SparseMatrix<double, Eigen::RowMajor> &stiff_mat, 
                                         const Eigen::SparseMatrix<double, Eigen::RowMajor> &mass_mat, 
                                         const Eigen::VectorXd &v_old, Eigen::VectorXd &v_new);


    /**
     * @brief 
     * 
     * @param thread_id 
     * @param conductivity_matrix_trans 
     * @param v_old 
     * @param v_new 
     */
    inline void DiffusionCallback(std::size_t thread_id, const Eigen::SparseMatrix<double, Eigen::RowMajor> &stiff_mat, 
                                  const Eigen::VectorXd &mass_vec, double dt, Eigen::VectorXd &v_nodal);

public:

    /**
     * \brief The default constructor of the Monodomain.
     */
    inline Monodomain();


    /**
     * \brief The destructor of the Monodomain.
     */
    inline ~Monodomain();


    /**
     * \brief Assembles the conductivity matrix and the capacitance vector for the geometry of the monodomain model for Finite Elements Method solution.
     * \param [in] mesh The mesh of the domain of interest.
     * \param [in] material The electrical material assigned to the domain of interest.
     * \return [void]
     */
    inline void AssembleMatrices(const IMP::Mesh<DIM, CELL_NODES> &mesh, const std::shared_ptr<ElectricBasic<DIM>> &material);


    /**
     * \overload
     * \brief Assemble the conductivity matrix and the capacitance vector for the geometry of the monodomain model for Meshfree Method solution.
     * \param [in] grid The point cloud of the domain of interest.
     * \param [in] material The electrical material assigned to the domain of interest.
     * \param [in] mfree_approx The meshfree approximant (e.g., MLS, RPI, etc.).
     * \return [void]
     */
    inline void AssembleMatrices(const IMP::Grid<DIM, CELL_NODES> &grid, const std::shared_ptr<ElectricBasic<DIM>> &material,
                                 const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx, const IMP::NodeSet &neumann_set);


    /**
     * \overload
     * \brief Assemble the conductivity matrix and the capacitance vector for the geometry of the monodomain model for Fragile Points Method solution.
     * \param [in] voro The voronoi tesselation of the domain of interest.
     * \param [in] material The electrical material assigned to the domain of interest.
     * \param [in] fpm The Fragile Points Method approximant.
     * \return [void]
     */
    inline void AssembleMatrices(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material, const CLOUDEA::Fpm<DIM> &fpm);


    /**
     * \brief Compute the stable time step for the diffusion term of the monodomain model.
     * \return [void]
     */
    inline void ComputeCriticalTimeStep();


    /**
     * \brief Compute the monodomain solution of the reaction-diffusion model with the operator-split method. 
     * \param [in] mesh The mesh of the domain of interest. 
     * \param [in] material The electrical material assigned to the domain of interest.
     * \param [in] fem_mats The FEM matrices (stiffness, mass) of the domain of interest.
     * \param [in] cell_model The cellular model to simulate the reaction term.
     * \param [in] stimulus The stimulus current to be applied on the stimulated region of the domain.
     */
    inline void Compute(const std::shared_ptr<ElectricBasic<DIM>> &material, const std::vector<Stimulus> &stimuli);


    /**
     * \brief 
     * \param [in] mfree_approx
     * \return [void] 
     */
    inline void FictitiousValuesToReal(const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx);


    /**
     * \brief Get the conductivity matrix of the monodomain model.
     * \return [const Eigen::SparseMatrix<double, Eigen::RowMajor>&] The conductivity matrix of the monodomain model.
     */
    inline const Eigen::SparseMatrix<double, Eigen::RowMajor> & ConductivityMat() const { return this->stiff_mat_; }


    /**
     * \brief Get the capacitance lumped vector of the monodomain model.
     * \return [const Eigen::VectorXd&] The capacitance lumped vector of the monodomain model.
     */
    inline const Eigen::VectorXd & CapacitanceVec() const { return this->mass_vec_; }


    /**
     * \brief Get the stored potential values for all the recorded time steps.
     * \return [const std::vector<Eigen::VectorXd>&] The stored potential values for all the recorded time steps. 
     */
    inline const std::vector<Eigen::VectorXd> & Vout() const { return this->vout_; }


    /**
     * \brief Get the stored potential values at a given recorded time step.
     * Fast access without range check.
     * \param [in] step_id The index of the recorded time step to get the stored potential values.  
     * \return [const Eigen::VectorXd&] The stored potential values at a given recorded time step.
     */
    inline const Eigen::VectorXd & Vout(std::size_t i) const { return this->vout_[i]; }


    /**
     * \brief Get the stored potential values at a given recorded time step.
     * Slower access with range check.
     * \param [in] step_id The index of the recorded time step to get the stored potential values.  
     * \return [const Eigen::VectorXd&] The stored potential values at a given recorded time step.
     */
    inline const Eigen::VectorXd & VoutAt(std::size_t i) const { return this->vout_.at(i); }


    /**
     * \brief Get the number of the used threads for the monodomain model solution.
     * \return [const std::size_t&] The number of the used threads for the monodomain model solution.
     */
    inline const std::size_t & ThreadsNumber() const { return this->threads_number_; }


};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_ENGINE_PHYSICS_MONODOMAIN_HPP_

#include "ELECTRA/engine/physics/monodomain.tpp"