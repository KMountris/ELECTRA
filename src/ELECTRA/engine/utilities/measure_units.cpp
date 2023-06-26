/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/utilities/measure_units.hpp"


namespace ELECTRA {

MeasureUnits::MeasureUnits() : ref_time_si_value_(1.), ref_length_si_value_(1.),
                               ref_conductance_si_value_(1.), ref_capacitance_si_value_(1.), ref_current_si_value_(1.), units_()
{
    // Set the default values of the units.
    this->SetTimeUnits();
    this->SetLengthUnits();
    this->SetConductanceUnits();
    this->SetCapacitanceUnits();
    this->SetCurrentUnits();
}


MeasureUnits::~MeasureUnits()
{}


void MeasureUnits::SetRefTimeSIValue(double val) noexcept
{
    this->ref_time_si_value_ = val;

    // Reset the values of the time units to account for the new reference unit.
    this->SetTimeUnits();
}


void MeasureUnits::SetRefLengthSIValue(double val) noexcept
{
    this->ref_length_si_value_ = val;

    // Reset the values of the length units to account for the new reference unit.
    this->SetLengthUnits();
}


void MeasureUnits::SetRefConductanceSIValue(double val) noexcept
{
    this->ref_conductance_si_value_ = val;

    // Reset the values of the conductance units to account for the new reference unit.
    this->SetConductanceUnits();
}


void MeasureUnits::SetRefCapacitanceSIValue(double val) noexcept
{
    this->ref_capacitance_si_value_ = val;

    // Reset the values of the capacitance units to account for the new reference unit.
    this->SetCapacitanceUnits();
}


void MeasureUnits::SetRefCurrentSIValue(double val) noexcept
{
    this->ref_current_si_value_ = val;

    // Reset the values of the capacitance units to account for the new reference unit.
    this->SetCurrentUnits();
}


const double & MeasureUnits::operator[] (const std::string &key) const
{
    return this->units_.at(key);
}


void MeasureUnits::SetTimeUnits() noexcept
{

    this->units_["h"]   = 3600. / this->ref_time_si_value_;
    this->units_["min"] = 60. / this->ref_time_si_value_;
    this->units_["s"]   = 1. / this->ref_time_si_value_;
    this->units_["ds"]  = 1.E-1 / this->ref_time_si_value_;
    this->units_["cs"]  = 1.E-2 / this->ref_time_si_value_;
    this->units_["ms"]  = 1.E-3 / this->ref_time_si_value_;
    this->units_["us"]  = 1.E-6 / this->ref_time_si_value_;
    this->units_["ns"]  = 1.E-9 / this->ref_time_si_value_;
    this->units_["ps"]  = 1.E-12 / this->ref_time_si_value_;
    this->units_["fs"]  = 1.E-15 / this->ref_time_si_value_;
}


void MeasureUnits::SetLengthUnits() noexcept
{
    // Lengh units.
    this->units_["Mm"]  = 1.E6 / this->ref_length_si_value_;
    this->units_["km"]  = 1.E3 / this->ref_length_si_value_;
    this->units_["hm"]  = 1.E2 / this->ref_length_si_value_;
    this->units_["dam"] = 1.E1 / this->ref_length_si_value_;
    this->units_["m"]   = 1. / this->ref_length_si_value_;
    this->units_["dm"]  = 1.E-1 / this->ref_length_si_value_;
    this->units_["cm"]  = 1.E-2 / this->ref_length_si_value_;
    this->units_["mm"]  = 1.E-3 / this->ref_length_si_value_;
    this->units_["um"]  = 1.E-6 / this->ref_length_si_value_;
    this->units_["nm"]  = 1.E-9 / this->ref_length_si_value_;
    this->units_["pm"]  = 1.E-12 / this->ref_length_si_value_;
    this->units_["fm"]  = 1.E-15 / this->ref_length_si_value_;

    // Surface units.
    this->units_["Mm2"]  = this->units_.at("Mm")*this->units_.at("Mm");
    this->units_["km2"]  = this->units_.at("km")*this->units_.at("km");
    this->units_["hm2"]  = this->units_.at("hm")*this->units_.at("hm");
    this->units_["dam2"] = this->units_.at("dam")*this->units_.at("dam");
    this->units_["m2"]   = this->units_.at("m")*this->units_.at("m");
    this->units_["dm2"]  = this->units_.at("dm")*this->units_.at("dm");
    this->units_["cm2"]  = this->units_.at("cm")*this->units_.at("cm");
    this->units_["mm2"]  = this->units_.at("mm")*this->units_.at("mm");
    this->units_["um2"]  = this->units_.at("um")*this->units_.at("um");
    this->units_["nm2"]  = this->units_.at("nm")*this->units_.at("nm");
    this->units_["pm2"]  = this->units_.at("pm")*this->units_.at("pm");
    this->units_["fm2"]  = this->units_.at("fm")*this->units_.at("fm");

    // Volume units.
    this->units_["Mm3"]  = this->units_.at("Mm")*this->units_.at("Mm")*this->units_.at("Mm");
    this->units_["km3"]  = this->units_.at("km")*this->units_.at("km")*this->units_.at("km");
    this->units_["hm3"]  = this->units_.at("hm")*this->units_.at("hm")*this->units_.at("hm");
    this->units_["dam3"] = this->units_.at("dam")*this->units_.at("dam")*this->units_.at("dam");
    this->units_["m3"]   = this->units_.at("m")*this->units_.at("m")*this->units_.at("m");
    this->units_["dm3"]  = this->units_.at("dm")*this->units_.at("dm")*this->units_.at("dm");
    this->units_["cm3"]  = this->units_.at("cm")*this->units_.at("cm")*this->units_.at("cm");
    this->units_["mm3"]  = this->units_.at("mm")*this->units_.at("mm")*this->units_.at("mm");
    this->units_["um3"]  = this->units_.at("um")*this->units_.at("um")*this->units_.at("um");
    this->units_["nm3"]  = this->units_.at("nm")*this->units_.at("nm")*this->units_.at("nm");
    this->units_["pm3"]  = this->units_.at("pm")*this->units_.at("pm")*this->units_.at("pm");
    this->units_["fm3"]  = this->units_.at("fm")*this->units_.at("fm")*this->units_.at("fm");

}


void MeasureUnits::SetConductanceUnits() noexcept
{
    this->units_["MS"] = 1.E6 / this->ref_conductance_si_value_;
    this->units_["kS"] = 1.E3 / this->ref_conductance_si_value_;
    this->units_["S"]  = 1. / this->ref_conductance_si_value_;
    this->units_["dS"] = 1.E-1 / this->ref_conductance_si_value_;
    this->units_["cS"] = 1.E-2 / this->ref_conductance_si_value_;
    this->units_["mS"] = 1.E-3 / this->ref_conductance_si_value_;
    this->units_["uS"] = 1.E-6 / this->ref_conductance_si_value_;
    this->units_["nS"] = 1.E-9 / this->ref_conductance_si_value_;
    this->units_["pS"] = 1.E-12 / this->ref_conductance_si_value_;
    this->units_["fS"] = 1.E-15 / this->ref_conductance_si_value_;
}


void MeasureUnits::SetCapacitanceUnits() noexcept
{
    this->units_["MF"] = 1.E6 / this->ref_capacitance_si_value_;
    this->units_["kF"] = 1.E3 / this->ref_capacitance_si_value_;
    this->units_["F"]  = 1. / this->ref_capacitance_si_value_;
    this->units_["dF"] = 1.E-1 / this->ref_capacitance_si_value_;
    this->units_["cF"] = 1.E-2 / this->ref_capacitance_si_value_;
    this->units_["mF"] = 1.E-3 / this->ref_capacitance_si_value_;
    this->units_["uF"] = 1.E-6 / this->ref_capacitance_si_value_;
    this->units_["nF"] = 1.E-9 / this->ref_capacitance_si_value_;
    this->units_["pF"] = 1.E-12 / this->ref_capacitance_si_value_;
    this->units_["fF"] = 1.E-15 / this->ref_capacitance_si_value_;
}


void MeasureUnits::SetCurrentUnits() noexcept
{
    this->units_["MA"] = 1.E6 / this->ref_current_si_value_;
    this->units_["kA"] = 1.E3 / this->ref_current_si_value_;
    this->units_["A"]  = 1. / this->ref_current_si_value_;
    this->units_["dA"] = 1.E-1 / this->ref_current_si_value_;
    this->units_["cA"] = 1.E-2 / this->ref_current_si_value_;
    this->units_["mA"] = 1.E-3 / this->ref_current_si_value_;
    this->units_["uA"] = 1.E-6 / this->ref_current_si_value_;
    this->units_["nA"] = 1.E-9 / this->ref_current_si_value_;
    this->units_["pA"] = 1.E-12 / this->ref_current_si_value_;
    this->units_["fA"] = 1.E-15 / this->ref_current_si_value_;
}


} // End of namespace ELECTRA