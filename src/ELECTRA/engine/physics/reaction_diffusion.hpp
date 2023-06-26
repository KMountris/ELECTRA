/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


/**
   \file reaction_diffusion.hpp
   \brief Base class of reaction diffusion physics models.
   \author Konstantinos A. Mountris
   \date 15/03/2019
*/

#ifndef ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_HPP_
#define ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_HPP_

#include "ELECTRA/engine/materials/electric_basic.hpp"
#include "ELECTRA/engine/conduction_system/conduction_system.hpp"
#include "ELECTRA/engine/electrophysiology/ep_basic.hpp"
#include "ELECTRA/engine/electrophysiology/ep_varying_params.hpp"

#include <CLOUDEA/CLOUDEA>
#include <IMP/Tesselations>
#include <IMP/Topology>

#include <Eigen/Eigen>

#include <string>
#include <memory>
#include <stdexcept>
#include <exception>

namespace ELECTRA {

/** \addtogroup Physics \{ */

/**
 * \enum ReactionDiffusionType
 * \author Konstantinos A. Mountris
 * \brief Enumeration to declare the type of the reaction diffusion model.
 */
enum class ReactionDiffusionType {
    unknown,         /**< Unknown type of reaction diffusion model */
    bidomain,       /**< Bidomain model */
    monodomain     /**< Monodomain model */
};


/**
 * \enum IntegrationType
 * \author Konstantinos A. Mountris
 * \brief Enumeration to declare the type of the time integration.
 */
enum class IntegrationType {
    unknown,          /**< Unknown time integration */
    implicitly,       /**< Implicit time integration */
    explicitly        /**< Implicit time integration */
};


/**
 * \enum DiffusionSolverType
 * \author Konstantinos A. Mountris
 * \brief Enumeration to declare the type of the diffusion term solver.
 */
enum class DiffusionSolverType {
    unknown,        /**< Unknown solver */
    fe,             /**< Forward Euler solver */
    be,             /**< Backward Euler solver */
    cn,             /**< Crank-Nicolson solver */
    daeti           /**< Dual adaptive explicit time integration solver */
};


/**
 * \class ReactionDiffusion
 * \brief Base class of the general reaction-diffusion model for cardiac electrophysiology problems.
 * \tparam DIM The dimensions of the reaction diffusion model's domain. Supported: [1 | 2 | 3].
 * \tparam CELL_NODES The number of nodes composing the cells of the reaction diffusion's domain.
 */
template<short DIM, short CELL_NODES=1>
class ReactionDiffusion
{

private:

    std::vector<std::unique_ptr<EpBasic>> cells_;               /**< The container of the cells assigned to the geometry of the monodomain model */

    std::vector<EpVaryingParams> cell_var_param_groups_;        /**< The groups of parameters values for cells with time-dependent parameters */

    std::vector<int> cell_var_param_flags_;                     /**< A flag indicating which cells have time-dependent parameters handling */

    ConductionSystem<DIM> conduct_sys_;                         /**< The cardiac conduction system of the model */

    double simulation_time_;                                    /**< The time of the reaction-diffusion simulation */

    double dt_critical_;                                        /**< The critical time step for the diffusion term of the reaction-diffusion simulation */

    double dt_min_;                                             /**< The minimum time step for the reaction-diffusion simulation */

    double dt_;                                                 /**< The used time step for the reaction-diffusion simulation */

    int output_steps_;                                          /**< The number of steps to output the potential nodal values */

    int simulation_steps_;                                      /**< The total number of the simulation's time integration steps */

    IntegrationType integration_type_;                          /**< The time integration type */

    DiffusionSolverType diff_solver_type_;                      /**< The diffusion term solver type */

    bool adaptive_reaction_;                                    /**< Boolean for the selection of adaptive explicit integration in time for the reaction term */

    bool adaptive_diffusion_;                                   /**< Boolean for the selection of adaptive explicit integration in time for the diffusion term */


public:

    /**
     * \brief The default constructor of the ReactionDiffusion class.
     */
    explicit ReactionDiffusion();


    /**
     * \brief The destructor of the ReactionDiffusion class.
     */
    virtual ~ReactionDiffusion();


    /**
     * \brief Set the simulation time of the reaction diffusion model.
     * \param [in] time The simulation time of the reaction diffusion model.
     * \return [void]
     */
    inline void SetSimulationTime(double time);


    /**
     * \brief Set the minimum time step of time integration.
     * \param [in] dt_min The minimum time step of time integration.
     * \return [void]
     */
    inline void SetDtMin(double dt_min);


    /**
     * \brief Set the critical time step of time integration.
     * \param [in] dt_crit The critical time step of time integration.
     * \return [void]
     */
    inline void SetDtCritical(double dt_crit);


    /**
     * \brief Set the simulation time step of time integration.
     * \param [in] dt_crit The simulation time step of time integration.
     * \return [void]
     */
    inline void SetDt(double dt);


    /**
     * \brief Set the time integration type of the simulation.
     * \param [in] integration_type The time integration type of the simulation.
     * \return [void]
     */
    inline void SetIntegrationType(IntegrationType integration_type);


    /**
     * \brief Set the diffusion term solver type.
     * \param [in] diff_solver_type The type of the diffusion term solver.
     * \return [void]
     */
    inline void SetDiffusionSolverType(DiffusionSolverType diff_solver_type);


    /**
     * \brief Set the status for adaptive integration of the reaction term.
     * \param [in] status The status for adaptive integration of the reaction term.
     * \return [void]
     */
    inline void SetAdaptiveReaction(bool status);


    /**
     * \brief Set the status for adaptive integration of the diffusion term.
     * \param [in] status The status for adaptive integration of the diffusion term.
     * \return [void]
     */
    inline void SetAdaptiveDiffusion(bool status);


    /**
     * \brief Set the status for adaptive integration of both the reaction and diffusion terms.
     * \param [in] status The status for adaptive integration of both the reaction and diffusion terms.
     * \return [void]
     */
    inline void SetAdaptiveIntegration(bool status);


    /**
     * \brief Set the number of steps for output of the stored results.
     * \param [in] output_steps The number for output of the stored results.
     * \return [void]
     */
    inline void SetOutputSteps(int output_steps);


    /**
     * \brief Set up the time step to be used during the time integration.
     * \return [void]
     */
    inline void SetupTimeStepping();


    /**
     * \brief Set the container of the nodal cells by moving another container.
     * \param [in] nodal_cells The nodal cells to be moved. They can be of any cell model given that is inheriting the EpBasic.
     * \return [void]
     */
    inline void SetNodalCells(std::vector<std::unique_ptr<EpBasic>> nodal_cells);


    /**
     * \brief Set the time varying parameters of different nodal cells.
     * \param [in] cardiac_cell_var_param_groups The time varying parameters of different nodal cells.
     * \return [void]
     */
    inline void SetCellVarParamGroups(const std::vector<EpVaryingParams> &cardiac_cell_var_param_groups);


    /**
     * \brief Assembles the conductivity matrix and the capacitance vector for the geometry of the monodomain model for Finite Elements Method solution.
     * \param [in] mesh The mesh of the domain of interest.
     * \param [in] material The electrical material assigned to the domain of interest.
     * \return [void]
     */
    virtual void AssembleMatrices(const IMP::Mesh<DIM, CELL_NODES> &mesh, const std::shared_ptr<ElectricBasic<DIM>> &material) = 0;


    /**
     * \overload
     * \brief Assemble the conductivity matrix and the capacitance vector for the geometry of the monodomain model for Meshfree Method solution.
     * \param [in] grid The point cloud of the domain of interest.
     * \param [in] material The electrical material assigned to the domain of interest.
     * \param [in] mfree_approx The meshfree approximant (e.g., MLS, RPI, etc.).
     * \return [void]
     */
    virtual void AssembleMatrices(const IMP::Grid<DIM, CELL_NODES> &grid, const std::shared_ptr<ElectricBasic<DIM>> &material,
            const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx, const IMP::NodeSet &neumann_set) = 0;


    /**
     * \overload
     * \brief Assemble the conductivity matrix and the capacitance vector for the geometry of the monodomain model for Fragile Points Method solution.
     * \param [in] voro The voronoi tesselation of the domain of interest.
     * \param [in] material The electrical material assigned to the domain of interest.
     * \param [in] fpm The Fragile Points Method approximant.
     * \return [void]
     */
    virtual void AssembleMatrices(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
            const CLOUDEA::Fpm<DIM> &fpm) = 0;


    /**
     * \brief Compute the stable time step for the diffusion term of the reaction diffusion model.
     * \return [void]
     */
    virtual void ComputeCriticalTimeStep() = 0;


    /**
     * \brief Compute the solution of the reaction-diffusion model with the operator-split method.
     * \param [in] material The electrical material assigned to the domain of interest.
     * \param [in] stimuli The applied stimuli on different regions of the domain.
     */
    virtual void Compute(const std::shared_ptr<ElectricBasic<DIM>> &material, const std::vector<Stimulus> &stimuli) = 0;


    /**
     * \brief Convert fictitious voltage values that were computed with the Mixed Collocation method to real.
     * \param [in] mfree_approx The meshfree approximation to be used for the conversion.
     * \return [void]
     */
    virtual void FictitiousValuesToReal(const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx) = 0;


    /**
     * \brief Get the stable time step for the diffusion term of the reaction diffusion model.
     * \return [double] The stable time step for the diffusion term of the reaction diffusion model.
     */
    inline double DtCritical() const { return this->dt_critical_; }


    /**
     * \brief Get the simulation time of the reaction diffusion model.
     * \return [double] The simulation time of the reaction diffusion model.
     */
    inline double SimulationTime() const { return this->simulation_time_; }


    /**
     * \brief The minimum time step of the reaction diffusion model.
     * \return [double] The minimum time step of the reaction diffusion model.
     */
    inline double DtMin() const { return this->dt_min_; }


    /**
     * \brief The used time step of the reaction diffusion model.
     * \return [double] The used time step of the reaction diffusion model.
     */
    inline double Dt() const { return this->dt_; }


    /**
     * \brief Get the status of the adaptive integration of the reaction term of the reaction diffusion model.
     * \return [true] Adaptive integration of the reaction term of the reaction diffusion model is enabled.
     * \return [false] Adaptive integration of the reaction term of the reaction diffusion model is disabled.
     */
    inline bool AdaptiveReaction() const { return this->adaptive_reaction_; }


    /**
     * \brief Get the status of the adaptive integration of the diffusion term of the reaction diffusion model.
     * \return [true] Adaptive integration of the reaction term of the diffusion diffusion model is enabled.
     * \return [false] Adaptive integration of the reaction term of the diffusion diffusion model is disabled.
     */
    inline bool AdaptiveDiffusion() const { return this->adaptive_diffusion_; }


    /**
     * \brief Get the status of the adaptive integration of both terms of the reaction diffusion model.
     * \return [true] Adaptive integration of the reaction term of both terms of the reaction diffusion model is enabled.
     * \return [false] Adaptive integration of the reaction term of both terms of the reaction diffusion model is disabled.
     */
    inline bool AdaptiveIntegration() const { return this->adaptive_reaction_ && this->adaptive_diffusion_; }


    /**
     * \brief Get the number of output steps.
     * \return [int] The number of output steps.
     */
    inline int OutputSteps() const { return this->output_steps_; }


    /**
     * \brief Get the number of simulation steps.
     * \return [int] The number of simulation steps.
     */
    inline int SimulationSteps() const { return this->simulation_steps_; }


    /**
     * @brief 
     * 
     * @return IntegrationType 
     */
    inline IntegrationType Integration() const { return this->integration_type_; }


    /**
     * @brief 
     * 
     * @return DiffusionSolverType 
     */
    inline DiffusionSolverType DiffusionSolver() const { return this->diff_solver_type_; }


    /**
     * \brief Get the conduction system of the reaction diffusion model.
     * \return [const ConductionSystem<DIM>&] The conduction system of the reaction diffusion model.
     */
    inline const ConductionSystem<DIM> & ConductSystem() const { return this->conduct_sys_; }


    /**
     * \brief Get the conduction system of the reaction diffusion model.
     * \return [ConductionSystem<DIM>&] The conduction system of the reaction diffusion model.
     */
    inline ConductionSystem<DIM> & EditConductSystem() { return this->conduct_sys_; }


    /**
     * \brief Get the nodal cells container.
     * \return [const std::vector<std::unique_ptr<BaseCell> >&] The nodal cells container.
     */
    inline const std::vector<std::unique_ptr<EpBasic>> & Cells() const { return this->cells_; }


    /**
     * \brief Get the nodal cells container.
     * \return [std::vector<std::unique_ptr<BaseCell> >&] The nodal cells container.
     */
    inline std::vector<std::unique_ptr<EpBasic>> & EditCells() { return this->cells_; }


    /**
     * @brief 
     * @return const std::vector<EpVaryingParams>& 
     */
    inline const std::vector<EpVaryingParams> & CellVarParamGroups() const { return this->cell_var_param_groups_; }


    /**
     * @brief 
     * @return const std::vector<int>& 
     */
    inline const std::vector<int> & CellVarParamFlags() const { return this->cell_var_param_flags_; }


    /**
     * \brief Get the conductivity matrix of the monodomain model.
     * \return [const Eigen::SparseMatrix<double, Eigen::RowMajor>&] The conductivity matrix of the monodomain model.
     */
    virtual const Eigen::SparseMatrix<double, Eigen::RowMajor> & ConductivityMat() const = 0;


    /**
     * \brief Get the capacitance lumped vector of the monodomain model.
     * \return [const Eigen::VectorXd&] The capacitance lumped vector of the monodomain model.
     */
    virtual const Eigen::VectorXd & CapacitanceVec() const = 0;


    /**
     * \brief Get the stored potential values for all the recorded time steps.
     * \return [const std::vector<Eigen::VectorXd>&] The stored potential values for all the recorded time steps.
     */
    virtual const std::vector<Eigen::VectorXd> & Vout() const = 0;


    /**
     * \brief Get the stored potential values at a given recorded time step.
     * Fast access without range check.
     * \param [in] step_id The index of the recorded time step to get the stored potential values.
     * \return [const Eigen::VectorXd&] The stored potential values at a given recorded time step.
     */
    virtual const Eigen::VectorXd & Vout(std::size_t i) const = 0;


    /**
     * \brief Get the stored potential values at a given recorded time step.
     * Slower access with range check.
     * \param [in] step_id The index of the recorded time step to get the stored potential values.
     * \return [const Eigen::VectorXd&] The stored potential values at a given recorded time step.
     */
    virtual const Eigen::VectorXd & VoutAt(std::size_t i) const = 0;


    /**
     * \brief Get the number of the used threads for the monodomain model solution.
     * \return [const std::size_t&] The number of the used threads for the monodomain model solution.
     */
    virtual const std::size_t & ThreadsNumber() const = 0;

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_HPP_

#include "ELECTRA/engine/physics/reaction_diffusion.tpp"