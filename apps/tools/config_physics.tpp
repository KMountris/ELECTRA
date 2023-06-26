/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_TPP_

#include "config_physics.hpp"

namespace APP_ELECTRA
{

template<short DIM, short CELL_NODES>
ConfigPhysics<DIM,CELL_NODES>::ConfigPhysics() : react_diff_map_()
{
    this->react_diff_map_["bidomain"] = ELECTRA::ReactionDiffusionType::bidomain;
    this->react_diff_map_["monodomain"] = ELECTRA::ReactionDiffusionType::monodomain;
}


template<short DIM, short CELL_NODES>
ConfigPhysics<DIM,CELL_NODES>::~ConfigPhysics()
{}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::InitializeReactionDiffusion(const Parser &parser,
        std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff) const
{
    // Initialize reaction diffusion model.
    std::string physics_model_type = parser.GetValue<std::string>("physics.reaction-diffusion.type");
    std::transform(std::begin(physics_model_type), std::end(physics_model_type), std::begin(physics_model_type), ::tolower);
    react_diff = ELECTRA::ReactionDiffusionFactory<DIM,CELL_NODES>::Create(this->react_diff_map_.at(physics_model_type));
}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SetReactionDiffusion(const Parser &parser, const ELECTRA::MeasureUnits &units,
        std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const
{
    std::string physics_model_type = parser.GetValue<std::string>("physics.reaction-diffusion.type");
    stream << ELECTRA::Logger::Message("Physics model: reaction-diffusion\n");
    stream << ELECTRA::Logger::Message("Reaction-diffusion type: "+physics_model_type+"\n");
    std::transform(std::begin(physics_model_type), std::end(physics_model_type), std::begin(physics_model_type), ::tolower);

    std::string solver_type = parser.GetValue<std::string>("physics.reaction-diffusion.solver");
    std::transform(std::begin(solver_type), std::end(solver_type), std::begin(solver_type), ::tolower);

    if (solver_type == "fe") {
        react_diff->SetIntegrationType(IntegrationType::explicitly);
        react_diff->SetDiffusionSolverType(DiffusionSolverType::fe);
        stream << ELECTRA::Logger::Message("Time integration: explicitly\n");
        stream << ELECTRA::Logger::Message("Solver type: fe - Forward Euler\n");
    } else if (solver_type == "be") {
        react_diff->SetIntegrationType(IntegrationType::implicitly);
        react_diff->SetDiffusionSolverType(DiffusionSolverType::be);
        stream << ELECTRA::Logger::Message("Time integration: implicitly\n");
        stream << ELECTRA::Logger::Message("Solver type: be - Backward Euler\n");
    } else if (solver_type == "cn") {
        react_diff->SetIntegrationType(IntegrationType::implicitly);
        react_diff->SetDiffusionSolverType(DiffusionSolverType::cn);
        stream << ELECTRA::Logger::Message("Time integration: implicitly\n");
        stream << ELECTRA::Logger::Message("Solver type: cn - Crank-Nicolson\n");
    } else if (solver_type == "daeti") {
        react_diff->SetIntegrationType(IntegrationType::explicitly);
        react_diff->SetDiffusionSolverType(DiffusionSolverType::daeti);
        stream << ELECTRA::Logger::Message("Time integration: explicitly\n");
        stream << ELECTRA::Logger::Message("Solver type: daeti - Dual adaptive explicit time integration\n");
    } else {
        std::string error_msg = "Could not set reaction-diffusion due to unknown solver type. Supported: [fe | be | cn | daeti]";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Set the reaction diffusion physics models according its type.
    if (physics_model_type == "bidomain") {
        this->SetBidomain(parser, units, react_diff, stream);
    } else if (physics_model_type == "monodomain") {
        this->SetMonodomain(parser, units, react_diff, stream);
    } else {
        std::string error_msg = "Could not set reaction diffusion physics. Supported model type: [monodomain | bidomain]";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
}



////!  PROTECTED ////


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SetBidomain(const Parser &parser, const ELECTRA::MeasureUnits &units,
        std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const
{

    // Set time step for time integration.
    react_diff->SetDt(parser.GetValue<double>("physics.reaction-diffusion.dt") * units[parser.GetValue<std::string>("physics.reaction-diffusion.dt unit")]);
    react_diff->SetDtMin(parser.GetValue<double>("physics.reaction-diffusion.dt min") * units[parser.GetValue<std::string>("physics.reaction-diffusion.dt unit")]);

    // Set output interval.
    react_diff->SetOutputSteps(parser.GetValue<int>("physics.reaction-diffusion.output interval"));

    // Set simulation time.
    react_diff->SetSimulationTime(parser.GetValue<double>("physics.reaction-diffusion.simulation time") * units[parser.GetValue<std::string>("physics.reaction-diffusion.simulation time unit")]);

}


template<short DIM, short CELL_NODES>
void ConfigPhysics<DIM,CELL_NODES>::SetMonodomain(const Parser &parser, const ELECTRA::MeasureUnits &units,
        std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, std::ostream &stream) const
{

    // Set max and min time step for time integration.
    react_diff->SetDt(parser.GetValue<double>("physics.reaction-diffusion.dt") * units[parser.GetValue<std::string>("physics.reaction-diffusion.dt unit")]);
    react_diff->SetDtMin(parser.GetValue<double>("physics.reaction-diffusion.dt min") * units[parser.GetValue<std::string>("physics.reaction-diffusion.dt unit")]);

    // Set output interval.
    react_diff->SetOutputSteps(parser.GetValue<int>("physics.reaction-diffusion.output interval"));

    // Set simulation time.
    react_diff->SetSimulationTime(parser.GetValue<double>("physics.reaction-diffusion.simulation time") * units[parser.GetValue<std::string>("physics.reaction-diffusion.simulation time unit")]);

    // Set adaptive reaction on/off.
    std::string adaptive_reaction_status = parser.GetValue<std::string>("physics.reaction-diffusion.adaptive reaction");
    if (adaptive_reaction_status == "yes" || adaptive_reaction_status == "true" || adaptive_reaction_status == "1")
        react_diff->SetAdaptiveReaction(true);

    // Set adaptive diffusion on/off.
    std::string adaptive_diffusion_status = parser.GetValue<std::string>("physics.reaction-diffusion.adaptive diffusion");
    if (adaptive_diffusion_status == "yes" || adaptive_diffusion_status == "true" || adaptive_diffusion_status == "1")
        react_diff->SetAdaptiveDiffusion(true);

    stream << ELECTRA::Logger::Message("Adaptive integration of reaction term: ");
    react_diff->AdaptiveReaction() ? stream << "YES\n" : stream << "NO\n";

    stream << ELECTRA::Logger::Message("Adaptive integration of diffusion term: ");
    react_diff->AdaptiveDiffusion() ? stream << "YES\n" : stream << "NO\n";
}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_PHYSICS_TPP_