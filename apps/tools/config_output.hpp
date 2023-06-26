/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_output.hpp
   \brief ConfigOutput class header file.
   \author Konstantinos A. Mountris
   \date 27/06/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/Physics"
#include "ELECTRA/Fibers"
#include "ELECTRA/PostProcess"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <memory>
#include <set>

using namespace ELECTRA;

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigOutput
 * \brief Class to configure the output of simulation results.
 */
template<short DIM, short CELL_NODES>
class ConfigOutput {

protected:

    /**
     * @brief 
     * 
     * @param parser 
     * @param nodes 
     * @param cells 
     * @param monodomain 
     * @param stream 
     */
    void OutputToParaview(const Parser &parser, const std::vector<IMP::Vec<int(DIM), double>> &nodes, const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells,
            const std::shared_ptr<ReactionDiffusion<DIM, CELL_NODES>> &react_diff, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param nodes 
     * @param cells 
     * @param monodomain 
     * @param post_process 
     * @param stream 
     */
    void OutputToEnsight(const Parser &parser, const std::vector<IMP::Vec<int(DIM), double>> &nodes, const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells,
            const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff, const PostProcess &post_process, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param post_process 
     * @param stream 
     */
    void OutputToAscii(const Parser &parser, const PostProcess &post_process, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param monodomain 
     * @param stream 
     */
    void OutputToBinary(const Parser &parser, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff,
            std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param fibers 
     * @param output_filename 
     */
    void SaveFibers(const Eigen::MatrixXd &fibers, const std::string &output_filename) const;


    void SaveDistanceField(const Eigen::VectorXd &distance_field, const std::string &output_filename) const;

    void SaveDirectionField(const Eigen::MatrixXd &direction_field, const std::string &output_filename) const;


public:

    /**
     * \brief ConfigOutput object constructor.
     */
    ConfigOutput();


    /**
     * \brief ConfigOutput object destructor.
     */
    virtual ~ConfigOutput();


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void OutputGeneration(const Parser &parser, const std::vector<IMP::Vec<DIM,double>> &nodes,
            const std::vector<IMP::Cell<DIM,CELL_NODES>> &cells, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM,CELL_NODES>> &react_diff,
            const PostProcess &post_process, std::ostream &stream) const;


    /**
     * @brief 
     * 
     * @param parser 
     * @param mesh 
     * @param voro 
     * @param ldrbm 
     * @param stream 
     */
    void OutputFibers(const Parser &parser, IMP::Mesh<DIM,CELL_NODES> mesh, IMP::Voronoi<DIM> voro,
            const Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const;
};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_HPP_

#include "config_output.tpp"