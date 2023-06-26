/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/conditions/body_loads/stimulus.hpp"


namespace ELECTRA
{

Stimulus::Stimulus()  : id_(-1), node_set_name_(""), stimulated_nodes_(), start_(0.), end_(std::numeric_limits<double>::max()), duration_(0.), cycle_lengths_(), amplitude_(0.), has_stimulated_node_set_(false)
{}


Stimulus::~Stimulus()
{}


void Stimulus::SetStart(double start)
{
    // Check if the given starting time is positive.
    if (start < 0.) {
        throw std::invalid_argument(Logger::Error("The start time of a stimulus current is expected greater or equal to zero."));
    }

    this->start_ = start;

}


void Stimulus::SetEnd(double end)
{
    // Check if the given ending time is positive.
    if (end < 0.) {
        throw std::invalid_argument(Logger::Error("The ending time of a stimulus current is expected greater or equal to zero."));
    }

    this->end_ = end;

}


void Stimulus::SetDuration(double duration)
{
    // Check if the given duration time is positive.
    if (duration < 0.) {
        throw std::invalid_argument(Logger::Error("The duration time of a stimulus current is expected greater or equal to zero."));
    }

    this->duration_ = duration;

}


void Stimulus::SetCycleLengths(const std::vector<double> &cycle_lengths)
{
    // It is required that the duration of the stimulus is already set.
    if (this->duration_ < 2*std::numeric_limits<double>::epsilon()) {
        std::string error_str = "Could not set the stimulus cycle lengths. The stimulus time duration must be set first.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    // Set the cycle lengths of the stimulus.
    this->cycle_lengths_ = cycle_lengths;

    // Compute the total activation time and check cycle lengths integrity.
    double active_time = 0.;
    for (const auto &cl : this->cycle_lengths_) {
        // Check if the current cycle time length is positive as expected.
        if (cl < 0.) {
            throw std::invalid_argument(Logger::Error("Could not set the stimulus cycle lengths. The cycle lengths of a stimulus current are expected greater or equal to zero."));
        }

        // Sum in total activation time.
        active_time += cl;
    }

    // Set the end time of the stimulus. It equals to the active time plus its duration.
    // Adds a small increase to ensure full stimulation for last cycle.
    this->end_ = active_time + this->duration_ + 2*std::numeric_limits<double>::epsilon(); 
}


void Stimulus::SetCycleLengths(double cycle_length)
{
    // Check if the given cycle time length is positive.
    if (cycle_length < 0.) {
        throw std::invalid_argument(Logger::Error("The cycle length of a stimulus current is expected greater or equal to zero."));
    }

    this->cycle_lengths_.clear();
    this->cycle_lengths_.emplace_back(cycle_length);

}


bool Stimulus::IsActive(double current_time, double total_time) const
{
    // Compute the number of stimulus current activations.
    int activations_num = 1;
    if (this->cycle_lengths_.size() == 1) {
        // In this case we must compute the activations number based on the cycle length and ending time of the stimulus.
        if (this->cycle_lengths_[0] > 2.*std::numeric_limits<double>::epsilon() && total_time > this->cycle_lengths_[0]) {
            if (total_time < this->end_)
                activations_num = total_time / this->cycle_lengths_[0];
            else
                activations_num = this->end_ / this->cycle_lengths_[0];
        }
    } else {
        // In this case the activations number is equal to the number of cycle lengts.
        activations_num = static_cast<int>(this->cycle_lengths_.size());
    }

    // Check if the given time lies in any of the activation intervals of the stimulus.
    if (this->cycle_lengths_.size() == 1) {
        for (int i = 0; i != activations_num; ++i) {
            // True if the stimulus current is active for the given time.
            if (current_time >= (this->start_ + i*this->cycle_lengths_[0]) &&
                current_time <= (this->start_ + this->duration_ + i*this->cycle_lengths_[0]) &&
                current_time <= this->end_) { return true; }
        }
    } else {
        double pass_cycle_time = 0.;
        for (int i = 0; i != activations_num; ++i) {

            // True if the stimulus current is active for the given time.
            if (current_time >= (this->start_ + pass_cycle_time) &&
                current_time <= (this->start_ + this->duration_ + pass_cycle_time) &&
                current_time <= this->end_) { return true; }

            // Increase the passed cycle time.
            pass_cycle_time += this->cycle_lengths_[i];
        }
    }

    // Return false otherwise.
    return false;

}


void Stimulus::SetStimulatedNodes(const std::unordered_map<std::string,IMP::NodeSet> &model_node_sets, int model_nodes_num)
{
    // Check if a node set name attribute has been assigned to the stimulus.
    if (this->node_set_name_.empty()) {
        throw std::runtime_error(Logger::Error("Could not set stimulated nodes for stimulus. A name attribute has not been assigned."));
    }

    // Initialize the stimulated nodes flags to zero.
    this->stimulated_nodes_.resize(model_nodes_num);

    // Search the node set associated to the stimulus by name.
    bool stim_nset_found = false;
    for (const auto &node_set : model_node_sets) {

        if (this->node_set_name_ == node_set.second.Name()) {
            // Set flags to 1 for nodes belonging in the given nodeset.
            stim_nset_found = true;
            for (const auto &node_id : node_set.second.NodeIds()) {
                this->stimulated_nodes_.coeffRef(node_id) = 1;
            }

            // Stop the node sets search.
            break;
        }
    }

    if (stim_nset_found) {
        this->has_stimulated_node_set_ = true;
    }

}


void Stimulus::AddStimulatedNodesPadding(int padding)
{
    // Get the indices of the non zero stimulated nodes.
    std::vector<int> stim_ids;
    stim_ids.reserve(this->stimulated_nodes_.nonZeros());
    for (Eigen::SparseVector<int>::InnerIterator it(this->stimulated_nodes_); it; ++it) {
        stim_ids.emplace_back(it.index());
    }

    // Apply the padding to the stimulated nodes indices.
    this->stimulated_nodes_.setZero();
    for (const auto &id : stim_ids) {
        this->stimulated_nodes_.coeffRef(id+padding) = 1;
    }

}


} // End of namespace ELECTRA.
