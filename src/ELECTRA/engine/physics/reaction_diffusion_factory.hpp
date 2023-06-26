/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file cell_factory.hpp
   \brief CellFactory class header file.
   \author Konstantinos A. Mountris
   \date 23/10/2019
*/

#ifndef ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_FACTORY_HPP_
#define ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_FACTORY_HPP_

#include "ELECTRA/engine/physics/reaction_diffusion.hpp"
#include "ELECTRA/engine/physics/monodomain.hpp"
#include "ELECTRA/engine/physics/bidomain.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <memory>
#include <string>


namespace ELECTRA {

/**
 *  \addtogroup Physics \{ */

/**
 * \class ReactionDiffusionFactory
 * \brief Class implementing a factory for reaction diffusion model generation.
 */
template <short DIM, short CELL_NODES>
class ReactionDiffusionFactory
{
public:

    /**
     * \brief Creates a reaction diffusion model according to the provided type.
     * \param [in] rd_model_type The type of the desired reaction diffusion model.
     * \return [std::shared_ptr<ReactionDiffusion>] Shared pointer to the created reaction diffusion model.
     */
    inline static std::shared_ptr<ReactionDiffusion<DIM,CELL_NODES>> Create(ReactionDiffusionType rd_model_type) {

        std::shared_ptr<ReactionDiffusion<DIM,CELL_NODES>> rd_ptr;

        switch (rd_model_type) {
        case ReactionDiffusionType::bidomain :
            rd_ptr = std::make_shared<Bidomain<DIM,CELL_NODES>>();
            break;
        case ReactionDiffusionType::monodomain :
            rd_ptr = std::make_shared<Monodomain<DIM,CELL_NODES>>();
            break;
        default:
            std::string error_msg = Logger::Error("Could not create reaction diffusion model. Not supported type.");
            throw std::invalid_argument(error_msg);
            break;
        }
        return rd_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_FACTORY_HPP_