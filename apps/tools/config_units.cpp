/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#include "config_units.hpp"

namespace APP_ELECTRA
{


ConfigUnits::ConfigUnits()
{}


ConfigUnits::~ConfigUnits()
{}


void ConfigUnits::SetReferenceScale(const Parser &parser, ELECTRA::MeasureUnits &units, std::ostream& stream) const
{
    // Create the SI units.
    ELECTRA::MeasureUnits si_units;

    // Set the reference values to the scaled measure units.
    std::string ref_time_unit = parser.GetValue<std::string>("reference units.time");
    std::string ref_length_unit = parser.GetValue<std::string>("tissue.geometry.unit");
    std::string ref_capacitance_unit = parser.GetValue<std::string>("reference units.capacitance");
    std::string ref_current_unit = parser.GetValue<std::string>("reference units.current");

    units.SetRefTimeSIValue(si_units[ref_time_unit]);
    units.SetRefLengthSIValue(si_units[ref_length_unit]);
    units.SetRefCapacitanceSIValue(si_units[ref_capacitance_unit]);
    units.SetRefCurrentSIValue(si_units[ref_current_unit]);

    stream << ELECTRA::Logger::Message("Reference units\n");
    stream << "                  - time: " << ref_time_unit << "\n";
    stream << "                  - length: " << ref_length_unit << "\n";
    stream << "                  - capacitance: " << ref_capacitance_unit << "\n";
    stream << "                  - current: " << ref_current_unit << "\n";
}


} // end of namespace APP_ELECTRA