/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file ep_varying_params.hpp
   \brief Header file for EpVaryingParams class.
   \author Konstantinos A. Mountris
   \date 18/03/2020
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_EP_VARYING_PARAMS_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_EP_VARYING_PARAMS_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/conditions/load_curve.hpp"

#include <IMP/IMP>

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <exception>
#include <stdexcept>


namespace ELECTRA {

/** \addtogroup Electrophysiology \{ */

/**
 * \class EpVaryingParams
 * \author Konstantinos A. Mountris
 * \brief A data structure for storing data of time-varying electophysiology parameters for cardiac cells.
 */
class EpVaryingParams {

private:

    std::unordered_map<std::size_t, LoadCurve> loads_;      /**< The load curve data */

    std::vector<double> time_;                              /**< The time that correspond to the entries in the load curve data */

    std::vector<int> cell_ids_;                             /**< The cardiac cell indices for which the loads will be applied */


public:

    /**
     * \brief Default constructor of electrophysiology varying parameters.
     */
    EpVaryingParams() : loads_(), time_(), cell_ids_() {}


    /**
     * \brief Default destructor of electrophysiology varying parameters.
     */
    virtual ~EpVaryingParams() {}


    /**
     * \brief Set the load curves of the electrophysiology varying parameters.
     * \param [in] loads The load curves of the electrophysiology varying parameters.
     * \return [void]
     */
    inline void SetLoads(const std::unordered_map<std::size_t, LoadCurve> &loads) { this->loads_ = loads; }


    /**
     * \brief Set the time that correspond to the load curve entries of the electrophysiology varying parameters.
     * \param [in] time The time that correspond to the load curve entries of the electrophysiology varying parameters.
     * \return [void]
     */
    inline void SetTime(std::vector<double> time) { this->time_ = time; }


    /**
     * \brief Set the indices of the cardiac cells where the load curve data are applied.
     * \param [in] cell_ids The indices of the cardiac cells where the load curve data are applied.
     * \return [void]
     */
    inline void SetCellIds(const std::vector<int> cell_ids) { this->cell_ids_ = cell_ids; }


    /**
     * \brief Clear the data structures of the electrophysiology varying parameters.
     * \return [void]
     */
    inline void Clear() {
        this->loads_.clear();
        this->time_.clear(); this->time_.shrink_to_fit();
        this->cell_ids_.clear(); this->cell_ids_.shrink_to_fit();
    }


    /**
     * \brief Get the load curves of the electrophysiology varying parameters.
     * \return [const std::unordered_map<std::size_t, LoadCurve>&] The load curves of the electrophysiology varying parameters.
     */
    inline const std::unordered_map<std::size_t, LoadCurve> & Loads() const { return this->loads_; }


    /**
     * \brief Get the time corresponding to the entries of the load curves of the electrophysiology varying parameters
     * \return [const std::vector<double>&] The time corresponding to the entries of the load curves of the electrophysiology varying parameters.
     */
    inline const std::vector<double> & Time() const { return this->time_; }


    /**
     * \brief Get the cardiac cell indices where the load curves of the electrophysiology varying parameters are applied.
     * \return [const std::vector<int>&] The cardiac cell indices where the load curves of the electrophysiology varying parameters are applied.
     */
    inline const std::vector<int> & CellIds() const { return this->cell_ids_; }


    /**
     * \brief Get the number of the load curves of the electrophysiology varying parameters.
     * \return [int] The number of the load curves of the electrophysiology varying parameters.
     */
    inline int LoadsNum() const { return static_cast<int>(this->loads_.size()); }

};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_ELECTROPHYSIOLOGY_VARYING_AP_PARAMS_HPP_