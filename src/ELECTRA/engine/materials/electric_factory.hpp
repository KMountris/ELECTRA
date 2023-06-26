/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file electric_factory.hpp
   \brief ElectricFactory class header file.
   \author Konstantinos A. Mountris
   \date 08/11/2021
*/

#ifndef ELECTRA_ENGINE_MATERIALS_ELECTRIC_FACTORY_HPP_
#define ELECTRA_ENGINE_MATERIALS_ELECTRIC_FACTORY_HPP_

#include "ELECTRA/engine/materials/electric_basic.hpp"
#include "ELECTRA/engine/materials/electric_transversal.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <memory>
#include <string>


namespace ELECTRA {

/**
 *  \addtogroup Materials \{ */

/**
 * \class ElectricFactory
 * \brief Class implementing a factory for electric material generation.
 */
template <short DIM>
class ElectricFactory
{
public:

    /**
     * \brief Creates an electric material according to the provided type.
     * \param [in] mat_type The type of the desired electric material.
     * \return [std::shared_ptr<ElectricBasic<DIM>>] Shared pointer to the created electric material.
     */
    inline static std::shared_ptr<ElectricBasic<DIM>> Create(ElectricType mat_type) {

        std::shared_ptr<ElectricBasic<DIM>> mat_ptr;

        switch (mat_type) {
        case ElectricType::transversal :
            mat_ptr = std::make_shared<ElectricTransversal<DIM>>();
            break;
        // case ElectricType::orthotropic :
        //     mat_ptr = std::make_shared<Monodomain<DIM,CELL_NODES>>();
        //     break;
        default:
            std::string error_msg = Logger::Error("Could not create electric material. Not supported type.");
            throw std::invalid_argument(error_msg);
            break;
        }
        return mat_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ENGINE_MATERIALS_ELECTRIC_FACTORY_HPP_