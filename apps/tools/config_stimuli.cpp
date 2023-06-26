/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#include "config_stimuli.hpp"

namespace APP_ELECTRA
{

ConfigStimuli::ConfigStimuli()
{}


ConfigStimuli::~ConfigStimuli()
{}


void ConfigStimuli::SetStimuli(const Parser &parser, const std::string &body_type, const std::unordered_map<std::string, IMP::NodeSet> &node_sets,
                               int nodes_num, const ELECTRA::MeasureUnits &units, std::vector<ELECTRA::Stimulus> &stimuli, std::ostream &stream) const
{
    // Check that body type is valid.
    if (body_type != "tissue" && body_type != "conduction system") {
        std::string error_msg = "Could not set stimuli. Body type must be: tissue or conduction system.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
    }

    unsigned short stimuli_num = parser.GetValue<unsigned short>(body_type+".stimuli.stimuli number");
    stimuli.clear();
    stimuli.resize(stimuli_num);

    ELECTRA::Stimulus stimulus;
    for (unsigned short i = 1; i <= stimuli_num; ++i) {
        std::string stimulus_path = body_type+".stimuli.stimulus-" + std::to_string(i);
        std::string stim_time_unit = parser.GetValue<std::string>(stimulus_path+".time unit");
        std::string stim_current_unit = parser.GetValue<std::string>(stimulus_path+".current unit");
        std::string stim_nodeset_name = parser.GetValue<std::string>(stimulus_path+".nodeset");

        stimulus.SetAttributes(parser.GetValue<int>(stimulus_path+".id"), stim_nodeset_name);
        stimulus.SetStart(parser.GetValue<double>(stimulus_path+".start") * units[stim_time_unit]);
        if (parser.HasAttribute(stimulus_path+".end")) {
            stimulus.SetEnd(parser.GetValue<double>(stimulus_path+".end") * units[stim_time_unit]);
        }
        stimulus.SetDuration(parser.GetValue<double>(stimulus_path+".duration") * units[stim_time_unit]);
        stimulus.SetAmplitude(parser.GetValue<double>(stimulus_path+".amplitude") * units[stim_current_unit]);

        // Set the stimulated nodes set.
        if (node_sets.find(stim_nodeset_name) == node_sets.end()) {
            std::string error_msg = "Could not set stimulus-" + std::to_string(i) + ". The nodeset [" + stim_nodeset_name + "] does not exist.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_msg));
        }
        stimulus.SetStimulatedNodes(node_sets, nodes_num);
        if (!stimulus.HasStimulatedNodeSet()) {
            stream << termcolor::yellow << ELECTRA::Logger::Warning("Selected node set to apply "+body_type+" stimulus "+
                    std::to_string(i)+" was not found in the "+body_type+" mesh file.") << termcolor::reset << "\n";
        }

        // Set cycle length of the stimulus
        if (parser.IsSingleValue(stimulus_path+".cycle length")) {
            stimulus.SetCycleLengths(parser.GetValue<double>(stimulus_path+".cycle length") * units[stim_time_unit]);

        } else if (parser.IsMultiArray(stimulus_path+".cycle length")) {

            // Get cycles multiarray
            auto object = parser.GetObject(stimulus_path+".cycle length");
            auto entries = object.get<std::vector<std::pair<int,double>>>();

            // Obtain the time length of each cycle.
            std::vector<double> varying_cycle_lengths;
            for (const auto &entry : entries) {
                int cycles_num = entry.first;
                double cycle_length = entry.second;

                if (cycles_num < 1) {
                    std::string error_msg = "Could not set cycle length for stimulus-" + std::to_string(i) + ". The number of cycles must be greater than zero.";
                    throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
                }

                // Store the cycle lengths as many times indicate the given number.
                for (int cl = 0; cl != cycles_num; ++cl) { varying_cycle_lengths.emplace_back(cycle_length * units[stim_time_unit]); }
            }
            stimulus.SetCycleLengths(varying_cycle_lengths);

        } else {
            std::string error_msg = "Could not set "+body_type+" stimulus-"+std::to_string(i)+". Cycle length is given in unknown format.";
            throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
        }// End of Set cycle length of the stimulus.

        // Assign the current stimulus.
        stimuli[i-1] = stimulus;

        // Print information about the current stimulus.
        stream << ELECTRA::Logger::Message("Set up stimulus-") << i << "\n";
        stream << ELECTRA::Logger::Message("       - id: ") << stimulus.Id() << "\n";
        stream << ELECTRA::Logger::Message("       - nodeset:   ") << stimulus.NodeSetName() << "\n";
        stream << ELECTRA::Logger::Message("       - start:     ") << stimulus.Start() << " " << stim_time_unit << "\n";

        if (stimulus.End() < std::numeric_limits<double>::max())
            stream << ELECTRA::Logger::Message("       - end:       ") << stimulus.End() << " " << stim_time_unit << "\n";

        stream << ELECTRA::Logger::Message("       - duration:  ") << stimulus.Duration() << " " << stim_time_unit << "\n";

        // Print cycle length with appropriate format for constant or varying.
        if (stimulus.CycleLengths().size() == 1) {
            stream << ELECTRA::Logger::Message("       - cycle length: ") << stimulus.CycleLength(0) << " " << stim_time_unit << "\n";
        }
        else if (stimulus.CycleLengths().size() > 1) {
            stream << ELECTRA::Logger::Message("       - cycle length:\n");
            for (const auto &cl : stimulus.CycleLengths())
                stream << ELECTRA::Logger::Message("                    ") << cl << " " << stim_time_unit << "\n";
        }
        stream << ELECTRA::Logger::Message("       - amplitude: ") << stimulus.Amplitude() << " " << stim_current_unit << "\n";
    }
}


} // end of namespace APP_ELECTRA