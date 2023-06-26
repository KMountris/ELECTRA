/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file ep_factory.hpp
   \brief EpFactory class header file.
   \author Konstantinos A. Mountris
   \date 23/10/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_EP_FACTORY_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_EP_FACTORY_HPP_

#include "ELECTRA/engine/electrophysiology/ep_basic.hpp"
#include "ELECTRA/engine/electrophysiology/bueno.hpp"
#include "ELECTRA/engine/electrophysiology/ohara.hpp"
#include "ELECTRA/engine/electrophysiology/courtemanche.hpp"
#include "ELECTRA/engine/electrophysiology/grandi_atri.hpp"
#include "ELECTRA/engine/electrophysiology/gaur2021.hpp"
#include "ELECTRA/engine/electrophysiology/gong2020.hpp"
#include "ELECTRA/engine/electrophysiology/gong2020m.hpp"
#include "ELECTRA/engine/electrophysiology/maleckar2009.hpp"
#include "ELECTRA/engine/electrophysiology/maccannell.hpp"
#include "ELECTRA/engine/electrophysiology/stewart.hpp"
#include "ELECTRA/engine/electrophysiology/paci_ventri.hpp"
#include "ELECTRA/engine/electrophysiology/tentusscher2006.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <memory>
#include <string>


namespace ELECTRA {

/**
 *  \addtogroup Electrophysiology \{ */

/**
 * \class EpFactory
 * \brief Class implementing a factory for action potential models generation.
 */
class EpFactory
{
public:

    /**
     * \brief Creates an action potential model according to the provided type.
     * \param [in] ep_model_type The type of the desired action potential model.
     * \return [std::unique_ptr<EpBasic>] Unique pointer to the created action potential model.
     */
    inline static std::unique_ptr<EpBasic> Create(EpModelType ep_model_type) {

        std::unique_ptr<EpBasic> ep_ptr;

        switch (ep_model_type) {
        case EpModelType::Bueno :
            ep_ptr = std::make_unique<Bueno>();
            break;
        case EpModelType::OHara :
            ep_ptr = std::make_unique<Ohara>();
            break;
        case EpModelType::Gaur2021 :
            ep_ptr = std::make_unique<Gaur2021>();
            break;
        case EpModelType::Gong2020 :
            ep_ptr = std::make_unique<Gong2020>();
            break;
        case EpModelType::Gong2020m :
            ep_ptr = std::make_unique<Gong2020m>();
            break;
        case EpModelType::PaciVentri :
            ep_ptr = std::make_unique<PaciVentri>();
            break;
        case EpModelType::TenTusscher2006 :
            ep_ptr = std::make_unique<TenTusscher2006>();
            break;
        case EpModelType::MacCannell :
            ep_ptr = std::make_unique<Maccannell>();
            break;
        case EpModelType::Courtemanche :
            ep_ptr = std::make_unique<Courtemanche>();
            break;
        case EpModelType::GrandiAtri :
            ep_ptr = std::make_unique<GrandiAtri>();
            break;
        case EpModelType::Maleckar2009 :
            ep_ptr = std::make_unique<Maleckar2009>();
            break;
        case EpModelType::Stewart :
            ep_ptr = std::make_unique<Stewart>();
            break;
        default:
            std::string error_str = Logger::Error("Could not create electrophysiology model. Not supported type.");
            throw std::invalid_argument(error_str);
            break;
        }
        return ep_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_EP_FACTORY_HPP_