/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/post_process/post_process.hpp"


namespace ELECTRA
{

PostProcess::PostProcess() : action_potentials_(), apds_(), local_activation_times_()
{}


PostProcess::~PostProcess()
{}


void PostProcess::ComputeActionPotential(const std::vector<Eigen::VectorXd> &potentials, const std::vector<int> &ap_node_set, double record_dt)
{
    // Reset the action potentials container.
    this->action_potentials_.clear();

    // Iterate over the nodes of the action potential measurement nodeset.
    for (const auto &ap_node : ap_node_set) {

        // Initialize action potential matrix.
        Eigen::MatrixXd action_potential = Eigen::MatrixXd::Zero(potentials.size(), 2);

        // Iterate over the time states of the potential
        for (const auto &current_v : potentials) {
            // Get the index of the saved states.
            auto state_id = &current_v - &potentials[0];

            // Store the action potential at the given state for the requested ap node.
            action_potential.coeffRef(state_id,0) = static_cast<double>(state_id)*record_dt;
            action_potential.coeffRef(state_id,1) = current_v.coeff(ap_node);
        }

        // Store the action potential for the current ap node.
        this->action_potentials_.emplace_back(action_potential);
    }

}


void PostProcess::ComputeLocalActivationTimes(const std::vector<Eigen::VectorXd> &nodal_potential, double dt, double cycle_tstart, double cycle_tend)
{
    auto nodes_num = nodal_potential[0].rows();

    // Reset the local activation times container.
    this->local_activation_times_ = Eigen::VectorXd::Zero(nodes_num);
    this->local_activation_times_.setConstant(std::numeric_limits<double>::quiet_NaN());

    // Estimate first and last time bin.
    int available_tbins = static_cast<int>(nodal_potential.size());
    int tbins_start = static_cast<int>(std::ceil(cycle_tstart/dt));
    int tbins_end = static_cast<int>(std::ceil(cycle_tend/dt));
    if (tbins_start > available_tbins) { tbins_start = 0; }
    if (tbins_end > available_tbins) { tbins_end = available_tbins; }

    // Estimate number of time bins for the given time period.
    int tbins = tbins_end - tbins_start;

    // Iterate over the nodes of the local activation times measurement nodeset.
    double lat = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd ap = Eigen::MatrixXd::Zero(tbins, 2);
    Eigen::VectorXd ap_derivative = Eigen::VectorXd::Zero(tbins-1);
    for (int nid = 0; nid != nodes_num; ++nid) {

        // Iterate over the time bins to compute action potential for this node at the given cycle length.
        for (int i = tbins_start; i != tbins_end; ++i) {

            // Store the action potential values of the requested ap node during the cycle length.
            ap.coeffRef(i-tbins_start,0) = i*dt;
            ap.coeffRef(i-tbins_start,1) = nodal_potential[i].coeff(nid);
        }

        // Check if node is stimulated.
        Eigen::VectorXd::Index max_deriv_id;
        if (ap.col(1).maxCoeff() > 0.) {
            // Compute the time derivative of the action potential.
            ap_derivative = this->ComputeTimeDerivative(ap);
            ap_derivative.maxCoeff(&max_deriv_id);

            // Set local activation time.
            lat = static_cast<double>(max_deriv_id+1) * dt;

            // Store the action potential for the current ap node.
            this->local_activation_times_.coeffRef(nid) = lat;
        }
    }

}


void PostProcess::ComputeAPD(const std::vector<Eigen::VectorXd> &nodal_potential, double dt, double cycle_tstart, double cycle_tend, int apd_lvl)
{
    auto nodes_num = nodal_potential[0].rows();

    // Reset the local activation times container.
    Eigen::VectorXd apd(nodes_num);
    apd.setConstant(std::numeric_limits<double>::quiet_NaN());

    // Estimate first and last time bin.
    int available_tbins = static_cast<int>(nodal_potential.size());
    int tbins_start = static_cast<int>(std::ceil(cycle_tstart/dt));
    int tbins_end = static_cast<int>(std::ceil(cycle_tend/dt));
    if (tbins_start > available_tbins) { tbins_start = 0; }
    if (tbins_end > available_tbins) { tbins_end = available_tbins; }

    // Estimate number of time bins for the given time period.
    int tbins_num = tbins_end - tbins_start;

    // Iterate over the nodes of the nodeset nodes.
    Eigen::MatrixXd ap = Eigen::MatrixXd::Zero(tbins_num, 2);
    Eigen::VectorXd ap_deriv = Eigen::VectorXd::Zero(tbins_num-1);
    for (auto nid = 0; nid != nodes_num; ++nid) {

        // Iterate over the time bins to compute action potential for this node at the given cycle length.
        for (int i = tbins_start; i != tbins_end; ++i) {

            // Store the action potential values of the requested ap node during the cycle length.
            ap(i,0) = i*dt;
            ap(i,1) = nodal_potential[i].coeff(nid);
        }

        double ap_baseline = ap.coeff(0,1);
        Eigen::MatrixXd::Index ap_peak_id;
        double ap_peak = ap.col(1).maxCoeff(&ap_peak_id);

        // Compute the time derivative of the action potential.
        ap_deriv = this->ComputeTimeDerivative(ap);
        Eigen::VectorXd::Index max_deriv_id;
        ap_deriv.maxCoeff(&max_deriv_id);

        // Action potential value after the target APD.
        double voi = ap_baseline + ((100.0-static_cast<double>(apd_lvl))/100.0) * (ap_peak - ap_baseline);

        // Find index of the first value after the target APD.
        long int voi_id = -1;
        for (auto i = ap_peak_id; i != ap.rows(); ++i) {
            if (ap.coeff(i,1) < voi) {
                voi_id = i;
                break;
            }
        }

        if (voi_id > -1) {
            apd.coeffRef(nid) = ap.coeff(voi_id,0) - ap.coeff(max_deriv_id,0);
        }
    
    } // End of Iterate over the nodes of the nodeset nodes.

    // Add this apd measurement in the container.
    this->apds_[apd_lvl] = apd;

}


const Eigen::VectorXd PostProcess::ComputeTimeDerivative(const Eigen::MatrixXd &action_potential) const
{
    // Initialize the time derivative vector.
    Eigen::VectorXd dv_dt = Eigen::VectorXd::Zero(action_potential.rows()-1);

    // Compute the time derivative with finite differences.
    for (int i = 1; i != action_potential.rows(); ++i) {
        dv_dt.coeffRef(i-1) = (action_potential.coeff(i,1) - action_potential.coeff(i-1,1)) / (action_potential.coeff(i,0) - action_potential.coeff(i-1,0));
    }

    return dv_dt;

}


} // End of namespace ELECTRA
