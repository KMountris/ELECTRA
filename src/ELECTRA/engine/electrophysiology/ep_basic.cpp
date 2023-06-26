/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/electrophysiology/ep_basic.hpp"


namespace ELECTRA
{

void EpBasic::UpdateTimeVaryingPrms(const EpVaryingParams &var_params, double dt)
{
    // Get the index for the value of the load curve at the current time.
    auto iter_geq = std::lower_bound(std::begin(var_params.Time()), std::end(var_params.Time()), this->integ_time_);
    std::size_t id = iter_geq - std::begin(var_params.Time());
    if (id >= var_params.Time().size()) { id = var_params.Time().size()-1; }

    // The index should correspond to the load curve entry at time less or equal to the integration time.
    // We correct the index for larger integration time. For equal, the result is as expected.
    double eps = 1e-8;
    if (eps > dt) { eps = dt; }
    if (id != 0 && var_params.Time()[id] - this->integ_time_ > eps) { id--; }

    // Update all the time-varying parameters.
    for (const auto &param_load : var_params.Loads()) {
        this->SetPrm(param_load.first, param_load.second.Data(id));

        // if (this->integ_time_ < 0.05) {
        // std::cout << "int time: " << this->integ_time_ << " id: " << id << " prm first: " << param_load.first <<
        // " prm second: " << param_load.second.Data(id) << "\n";}
    }

    // Increase the integration time.
    // this->integ_time_ = ALGORITHM::KahanSum({this->integ_time_, dt});
    this->integ_time_ += dt;
}

} // End of namespace ELECTRA