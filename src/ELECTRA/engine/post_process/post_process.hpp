/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file post_process.hpp
   \brief PostProcess class header file.
   \author Konstantinos A. Mountris
   \date 06/05/2019
*/

#ifndef ELECTRA_POST_PROCESS_POST_PROCESS_HPP_
#define ELECTRA_POST_PROCESS_POST_PROCESS_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"

#include <Eigen/Dense>

#include <vector>
#include <iostream>
#include <string>
#include <iterator>
#include <unordered_map>

namespace ELECTRA {

/** \addtogroup Post-Process \{ */

/**
 * \class PostProcess
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting post processing tools for ELECTRA simulations.
 */
class PostProcess
{   
private:

    std::vector<Eigen::MatrixXd> action_potentials_;

    std::unordered_map<int, Eigen::VectorXd> apds_;

    Eigen::VectorXd local_activation_times_;    

public:

    PostProcess();


    virtual ~PostProcess();


    void ComputeActionPotential(const std::vector<Eigen::VectorXd> &potentials, const std::vector<int> &ap_node_set, double record_dt);


    /**
     * \brief Compute the local activation time (lat) for each node in the given node set for a given cycle. 
     * \param [in] nodal_potential The potential value at each node at each time step. 
     * \param [in] dt The size of the time step. 
     * \param [in] cycle_tstart The starting time of the cycle.
     * \param [in] cycle_tend The ending time of the cycle.
     */
    void ComputeLocalActivationTimes(const std::vector<Eigen::VectorXd> &nodal_potential, double dt, double cycle_tstart, double cycle_tend);


    void ComputeAPD(const std::vector<Eigen::VectorXd> &nodal_potential, double dt, double cycle_tstart, double cycle_tend, int apd_lvl);


    const Eigen::VectorXd ComputeTimeDerivative(const Eigen::MatrixXd &action_potential) const;


    inline const std::vector<Eigen::MatrixXd> & ActionPotentials() const { return this->action_potentials_; }


    inline const Eigen::VectorXd & LocalActivationTimes() const { return this->local_activation_times_; }

    
    inline const std::unordered_map<int, Eigen::VectorXd> & APDs() const { return this->apds_; }


};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_POST_PROCESS_POST_PROCESS_HPP_
