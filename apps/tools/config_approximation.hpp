/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_approximation.hpp
   \brief ConfigApproximation class header file.
   \author Konstantinos A. Mountris
   \date 24/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"

#include <IMP/IMP>
#include <CLOUDEA/CLOUDEA>

#include <termcolor/termcolor.hpp>
#include <Eigen/Dense>

#include <string>
#include <filesystem>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <functional>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigApproximation
 * \brief Class to configure the numerical approximation method of a simulation.
 */
template <short DIM, short CELL_NODES>
class ConfigApproximation {

public:

    /**
     * \brief Construct a new Config Approximation object
     */
    ConfigApproximation();


    /**
     * \brief Destroy the Config Approximation object
     */
    virtual ~ConfigApproximation();


    /**
     * @brief Set the Fem Approximation object
     * 
     * @param stream 
     */
    void SetFemApproximation(std::ostream &stream) const;


    void SetFpmApproximation(const Parser &parser, const IMP::Voronoi<DIM> &voro, CLOUDEA::Fpm<DIM> &fpm_approx, std::ostream &stream) const;


    void SetMcmApproximation(const Parser &parser, const IMP::Grid<DIM, CELL_NODES> &grid, std::unique_ptr<CLOUDEA::Mfree<DIM>> &mcm_approx,
                             IMP::NodeSet &neumann_nset, std::ostream &stream) const;


};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_APPROXIMATION_HPP_

#include "config_approximation.tpp"