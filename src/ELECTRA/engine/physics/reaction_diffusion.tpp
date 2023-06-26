/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#ifndef ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_TPP_
#define ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_TPP_


#include "ELECTRA/engine/physics/reaction_diffusion.hpp"


namespace ELECTRA {


template<short DIM, short CELL_NODES>
ReactionDiffusion<DIM, CELL_NODES>::ReactionDiffusion() : cells_(), cell_var_param_groups_(), cell_var_param_flags_(),
        conduct_sys_(), simulation_time_(-1.), dt_critical_(std::numeric_limits<double>::max()),
        dt_min_(-1.), dt_(-1.), output_steps_(-1), simulation_steps_(-1),
        integration_type_(IntegrationType::unknown), diff_solver_type_(DiffusionSolverType::unknown),
        adaptive_reaction_(false), adaptive_diffusion_(false)
{}


template<short DIM, short CELL_NODES>
ReactionDiffusion<DIM, CELL_NODES>::~ReactionDiffusion()
{}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetSimulationTime(double time)
{
    if (time <= 0.)
        throw std::invalid_argument(Logger::Error("Simulation time greater than zero is expected by the reaction diffusion model."));
    this->simulation_time_ = time;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetAdaptiveReaction(bool status)
{
    this->adaptive_reaction_ = status;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetAdaptiveDiffusion(bool status)
{
    this->adaptive_diffusion_ = status;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetAdaptiveIntegration(bool status)
{
    this->adaptive_reaction_ = status;
    this->adaptive_diffusion_ = status;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetOutputSteps(int output_steps)
{
    if (output_steps <= 0)
        throw std::invalid_argument(Logger::Error("A number of output steps greater than zero is expected by the reaction diffusion model."));
    this->output_steps_ = output_steps;

}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetupTimeStepping()
{
    if (this->integration_type_ == IntegrationType::unknown) {
        std::string error_msg = "Could not set up time stepping for reaction-diffusion. The integration type is unknown. Supported: [implicitly | explicitly]";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Choose simulation time step for explicit methods. For implicit, we use the user-defined step.
    if (this->integration_type_ == IntegrationType::explicitly) {
        if (!this->adaptive_reaction_) { this->dt_ = this->dt_min_; }

        // Check selected time step stability.
        if (this->dt_critical_ < this->dt_ && !this->adaptive_diffusion_) {
            this->dt_ = this->dt_critical_;
            std::cout << Logger::Warning("The provided time step is larger than the critical diffusion time step. "
                    "The critical time step will be used instead. Consider enabling the adaptive diffusion mode.\n");
        }
    }

    // Set number of total time steps.
    this->simulation_steps_ = static_cast<int>(std::ceil(this->simulation_time_ / this->dt_));

}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetDtMin(double dt_min)
{
    if (dt_min < 0.)
        throw std::invalid_argument(Logger::Error("Could not set minimum timestep in reaction diffusion model. A negative time step was provided."));
    this->dt_min_ = dt_min;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetDtCritical(double dt_crit)
{
    if (dt_crit < 0.)
        throw std::invalid_argument(Logger::Error("Could not set critical timestep in reaction diffusion model. A negative time step was provided."));
    this->dt_critical_ = dt_crit;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetDt(double dt)
{
    this->dt_ = dt;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetIntegrationType(IntegrationType integration_type)
{
    this->integration_type_ = integration_type;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetDiffusionSolverType(DiffusionSolverType diff_solver_type)
{
    this->diff_solver_type_ = diff_solver_type;
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetNodalCells(std::vector<std::unique_ptr<EpBasic>> nodal_cells)
{
    this->cells_.clear();
    this->cells_.reserve(nodal_cells.size());
    this->cells_ = std::move(nodal_cells);

    // Initialize flags to check for cells with time varying parameters.
    this->cell_var_param_flags_.clear();
    this->cell_var_param_flags_.resize(this->cells_.size(), -1);
}


template<short DIM, short CELL_NODES>
void ReactionDiffusion<DIM, CELL_NODES>::SetCellVarParamGroups(const std::vector<EpVaryingParams> &cardiac_cell_var_param_groups)
{
    // Set the groups of the varying parameters for the cardiac cells.
    this->cell_var_param_groups_ = cardiac_cell_var_param_groups;

    // Update the flags to check for cells with time varying parameters.
    for (std::size_t flag = 0; flag != this->cell_var_param_groups_.size(); ++flag) {
        for (const auto &cid : this->cell_var_param_groups_[flag].CellIds()) {
            this->cell_var_param_flags_[cid] = flag;
        }
    }

}

} // End of namespace ELECTRA


#endif //ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_TPP_