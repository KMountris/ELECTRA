/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file config_electric.hpp
   \brief ConfigElectric class header file.
   \author Konstantinos A. Mountris
   \date 17/05/2021
*/

#ifndef ELECTRA_APPS_TOOLS_CONFIG_ELECTRIC_HPP_
#define ELECTRA_APPS_TOOLS_CONFIG_ELECTRIC_HPP_

#include "parser.hpp"

#include "ELECTRA/Utilities"
#include "ELECTRA/Materials"

#include <IMP/IMP>
#include <termcolor/termcolor.hpp>

#include <Eigen/Dense>

#include <string>
#include <filesystem>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <memory>

namespace APP_ELECTRA {

/** \addtogroup Application-Tools \{ */

/**
 * \class ConfigElectric
 * \brief Class to configure the electric material properties of a model.
 * \tparam DIM The dimensions of the model.
 * \tparam CELL_NODES The number of nodes of the model's cells.
 */
template<short DIM, short CELL_NODES>
class ConfigElectric {

private:
    std::unordered_map<std::string, ELECTRA::ElectricType> electric_mat_map_;        /**< Map of the electric material types */

protected:

    /**
     * @brief 
     * @param parser 
     * @param model_nodesets 
     * @param values_num 
     * @param diff_values 
     */
    void AssignDiffusivity(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
            const std::string &diff_type, const ELECTRA::MeasureUnits &units, int values_num, std::vector<double> &diff_values) const;


    /**
     * @brief 
     * @param parser 
     * @param model_nodesets 
     * @param values_num 
     * @param diff_values 
     */
    void AssignTransversalRatio(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
                                int values_num, std::vector<double> &ratio_values) const;


    /**
     * @brief 
     * @param parser 
     * @param model_nodesets 
     * @param values_num 
     * @param diff_values 
     */
    void AssignCapacitance(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
                           const ELECTRA::MeasureUnits &units, int values_num, std::vector<double> &cap_values) const;

    
    /**
     * @brief 
     * @param parser 
     * @param model_nodesets 
     * @param values_num 
     * @param diff_values 
     */
    void AssignFibers(const Parser &parser, int values_num, Eigen::MatrixXd &fibers) const;


    /**
     * @brief 
     * @param parser 
     * @param attribute 
     * @param nodesets 
     * @param values_num 
     * @param values 
     */
    void ObtainValuesFromNodesets(const Parser &parser, const std::string &attribute, const std::unordered_map<std::string, IMP::NodeSet> &nodesets,
                                  int values_num, std::vector<double> &values) const;


public:

    /**
     * \brief ConfigElectric object constructor.
     */
    ConfigElectric();


    /**
     * \brief ConfigElectric object destructor.
     */
    virtual ~ConfigElectric();


    /**
     * @brief 
     * 
     * @param parser 
     * @param react_diff 
     * @param stream 
     */
    void InitializeMaterial(const Parser &parser, std::shared_ptr<ELECTRA::ElectricBasic<DIM>> &electric_mat) const;


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void SetMaterialProperties(const Parser &parser, const std::unordered_map<std::string, IMP::NodeSet> &model_nodesets,
            const ELECTRA::MeasureUnits &units, std::shared_ptr<ELECTRA::ElectricBasic<DIM>> &electric_mat, std::ostream &stream) const;


    /**
     * \brief Set up the reference scale of the units for proper unit conversion.
     * \param [in] parser The parser of the simulation file.
     * \param [out] units The units after setting up the reference scale.
     * \param [out] stream The output logging stream.
     * \return [void]
     */
    void ComputeDiffusionTensors(std::shared_ptr<ELECTRA::ElectricBasic<DIM>> &electric_mat, std::ostream &stream) const;

};

/** \} End of Doxygen Groups*/

} //end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_ELECTRIC_HPP_

#include "config_electric.tpp"