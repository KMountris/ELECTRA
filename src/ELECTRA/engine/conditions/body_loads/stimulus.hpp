/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


/**
   \file stimulus.hpp
   \brief Header file of a stimulus electric current body loading condition.
   \author Konstantinos A. Mountris
   \date 18/03/2019
*/

#ifndef ELECTRA_CONDITIONS_BODY_LOADS_STIMULUS_HPP_
#define ELECTRA_CONDITIONS_BODY_LOADS_STIMULUS_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"

#include <IMP/IMP>

#include <Eigen/Sparse>

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <exception>
#include <stdexcept>


namespace ELECTRA {

/** \addtogroup Conditions \{ */
/** \addtogroup Body-Loads \{ */


/**
 * \class Stimulus
 * \author Konstantinos A. Mountris
 * \brief An electric current used to stimulate a cellular model during cell / tissue / organ simulation.
 */
class Stimulus {

private:

    int id_;                                        /**< The associated index to the stimulus current */

    std::string node_set_name_;                     /**< The name label of the node set where the stimulus current is applied */

    Eigen::SparseVector<int> stimulated_nodes_;     /**< A vector of flags signaling if a node of a given topology (mesh or grid) will be stimulated by the given stimulus */

    double start_;                                  /**< The starting time of the stimulus current */

    double end_;                                    /**< The ending time of the stimulus current */

    double duration_;                               /**< The duration of the stimulus current */

    std::vector<double> cycle_lengths_;             /**< The time length of the stimulus cycles. Describe the periodicity of the stimulus current reactivation. Constant or varying periodicity is supported */

    double amplitude_;                              /**< The amplitude of the stimulus current */

    bool has_stimulated_node_set_;                  /**< Boolean flag to check if the stimulated nodes flag have been set */



public:

    /**
     * \brief The default constructor of the Stimulus.
     */
    Stimulus();


    /**
     * \brief The destructor of the Stimulus.
     */
    virtual ~Stimulus();


    /**
     * \brief Set the attributes of the stimulus.
     * \param [in] id The index attribute of the stimulus.
     * \param [in] node_set_name The name of the corresponding node set where the stimulus is applied.
     * \return [void]
     */
    inline void SetAttributes(int id, std::string node_set_name) { this->id_ = id; this->node_set_name_ = node_set_name; }


    /**
     * \brief Set the stimulated nodes set where the stimulus is applied.
     * \param [in] model_node_sets All the node sets of the geometry model where the stimulus will be applied.
     * \param [in] model_nodes_num The total number of nodes in the geometry model where the stimulus will be applied.
     */
    void SetStimulatedNodes(const std::unordered_map<std::string,IMP::NodeSet> &model_node_sets, int model_nodes_num);


    /**
     * \brief Set the starting time of the stimulus current.
     * \param [in] start The starting time of the stimulus current.
     * \return [void]
     */
    void SetStart(double start);


    /**
     * \brief Set the ending time of the stimulus current. By default is considered infinite.
     * \param [in] end The ending time of the stimulus current.
     * \return [void]
     */
    void SetEnd(double end);


    /**
     * \brief Set the time duration of the stimulus current.
     * \param [in] duration The time duration of the stimulus current.
     * \return [void]
     */
    void SetDuration(double duration);


    /**
     *\brief Set the time length of the stimulus current reactivation cycles.
     * The time length of all the cycles of the stimulus are set for the stimulus. Here varying time length
     * per cycle can be assigned.
     * \warning The stimulus activation terminates after time equal to the sum of the cycles' time lengths.
     * \param [in] cycle_lengths The time length of the stimulus current reactivation cycles.
     * \return [void]
     */
    void SetCycleLengths(const std::vector<double> &cycle_lengths);


    /**
     * \brief Set the time length of the stimulus current reactivation cycles.
     * The time length of all the cycles of the stimulus is considered constant and correspond to the given
     * cycle_length value.
     * \warning The stimulus activation remains during the total time of the simulation except if the end of the stimulus
     *          is set by the user.
     * \param [in] cycle_length The constant time length of all the stimulus current reactivation cycles.
     * \return [void]
     */
    void SetCycleLengths(double cycle_length);


    /**
     * \brief Set the amplitude of the stimulus current.
     * \param [in] amplitude The amplitude of the stimulus current.
     * \return [void]
     */
    inline void SetAmplitude(double amplitude) { this->amplitude_ = amplitude; }


    /**
     * \brief Add a padding to the stimulated nodes indices.
     * \param [in] padding The padding number to be applied at the stimulated nodes indices.
     * \return [void]
     */
    void AddStimulatedNodesPadding(int padding);


    /**
     * \brief Check if the stimulus current is active at the given time moment in the given total time.
     * \param [in] current_time The time moment for which the stimulus current activation status is checked.
     * \param [in] total_time The total time during which the stimulus current activations occur.
     * \return [true] The stimulus current is active at the given time.
     * \return [false] The stimulus current is not active at the given time.
     */
    bool IsActive(double current_time, double total_time) const;


    /**
     * \brief Get the index attribute of the stimulus.
     * \return [int] The index attribute of the stimulus.
     */
    inline int Id() const { return this->id_; }


    /**
     * \brief Get the name attribute of the corresponding node set where the stimulus is applied.
     * \return [const std::string&] The name attribute of the corresponding node set where the stimulus is applied.
     */
    inline const std::string & NodeSetName() const { return this->node_set_name_; }


    inline const Eigen::SparseVector<int> & StimulatedNodes() const { return this->stimulated_nodes_; }


    /**
     * \brief Get the starting time of the stimulus current. 
     * \return [double] The starting time of the stimulus current. 
     * \return [void] 
     */
    inline double Start() const { return this->start_; }


    /**
     * \brief Get the ending time of the stimulus current. 
     * \return [double] The ending time of the stimulus current. 
     * \return [void] 
     */
    inline double End() const { return this->end_; }


    /**
     * \brief Get the duration of the stimulus current. 
     * \return [double] The duration of the stimulus current. 
     * \return [void] 
     */
    inline double Duration() const { return this->duration_; }


    /**
     * \brief Get the time length for all the stimulus current reactivation cycles. 
     * \return [const std::vector<double> &] The time length of all the stimulus current reactivation cycles.
     */
    inline const std::vector<double> & CycleLengths() const { return this->cycle_lengths_; }


    /**
     * \brief Get the time length of a specific reactivation cycle of the stimulus current.
     * Fast access with no range check. 
     * \return [double] The time length of a specific reactivation cycle of the stimulus current.
     */
    inline double CycleLength(std::size_t id) const { return this->cycle_lengths_[id]; }


    /**
     * \brief Get the time length of a specific reactivation cycle of the stimulus current.
     * Slower access with range check.
     * \return [double] The time length of a specific reactivation cycle of the stimulus current.
     */
    inline double CycleLengthAt(std::size_t id) const { return this->cycle_lengths_.at(id); }


    /**
     * \brief Get the amplitude of the stimulus current. 
     * \return [double] The amplitude of the stimulus current. 
     * \return [void] 
     */
    inline double Amplitude() const { return this->amplitude_; }


    /**
     * \brief Gives the status of the nodeset of the stimulus.
     * \return [true] The nodeset for the stimulus has been set.
     * \return [false] The nodeset for the stimulus has not been set.
     */
    inline bool HasStimulatedNodeSet() const { return this->has_stimulated_node_set_; }

};


/** \} End of Doxygen Groups */
/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_CONDITIONS_BODY_LOADS_STIMULUS_HPP_