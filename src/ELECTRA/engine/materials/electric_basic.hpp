/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
   \file electric_basic.hpp
   \brief ElectricBasic class header file.
   \author Konstantinos A. Mountris
   \date 08/03/2019
*/

#ifndef ELECTRA_MATERIALS_ELECTRIC_BASIC_HPP_
#define ELECTRA_MATERIALS_ELECTRIC_BASIC_HPP_


#include "ELECTRA/engine/utilities/logger.hpp"

#include <Eigen/Dense>

#include <iostream>
#include <initializer_list>
#include <iterator>
#include <vector>
#include <numeric>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <exception>


namespace ELECTRA {

/** \addtogroup Materials \{ */


/**
 * \enum Type of finite element aproximants.
 * \author Konstantinos A. Mountris
 */
enum class ElectricType {
    unknown,          /**< Unknown electric material */
    transversal,      /**< Electric material with transversal anisotropy */
    orthotropic       /**< Electric material with full anisotropy */
};


/**
 * \class ElectricBasic
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting a material with electric properties for electric propagation simulations.
 * \tparam DIM The dimensions of the material's domain.
 */

template<short DIM>
class ElectricBasic
{

public:

    /**
     * \brief ElectricBasic default constructor.
     */
    explicit ElectricBasic() {}


    /**
     * \brief ElectricBasic default destructor.
     */
    virtual ~ElectricBasic() {}


    /**
     * \brief Set the number of the electric material's nodes.
     * \param [in] nodes_num The number of the nodes that belong to the electric material.
     * \return [void]
     */
    virtual void SetNodesNum(int nodes_num) = 0;


    /**
     * \brief Set the internal diffusivity of the material.
     * \param [in] diffusivity The nodal internal diffusivity values.
     * \return [void]
     */
    virtual void SetInterDiffusivity(const std::vector<double> &diffusivity) = 0;


    /**
     * \brief Set the external diffusivity of the material.
     * \param [in] diffusivity The nodal external diffusivity values.
     * \return [void]
     */
    virtual void SetExterDiffusivity(const std::vector<double> &diffusivity) = 0;


    /**
     * \brief Set the transmembrane diffusivity of the material.
     * \param [in] diffusivity The nodal transmembrane diffusivity values.
     * \return [void]
     */
    virtual void SetTransDiffusivity(const std::vector<double> &diffusivity) = 0;


    /**
     * \brief Set the transversal to longitudinal diffusivity ratio of the material.
     * \param [in] trans_ratio The value of the material's transversal-longitudinal diffusivity ratio.
     * \return [void]
     */
    virtual void SetTransversalRatio(const std::vector<double> &trans_ratio) = 0;


    /**
     * \brief Set the capacitance of the material.
     * \param [in] capacitance The value of the material's capacitance.
     * \return [void]
     */
    virtual void SetCapacitance(const std::vector<double> &capacitance) = 0;


    /**
     * \brief Set the direction of the material's fibers.
     * The fibers direction is expressed by a initializer list collecting the axial components
     * of the fiber direction. A value of 0.0 on a specific axis means no fiber direction
     * on this axis and a value of 1.0 means that the fibers direction is parallel to this
     * axis. The sum of the axial components should be equal to 1.0.
     * \param [in] direction The initializer list of the material's fibers direction.
     * \return [void]
     */
    virtual void SetLongFibers(std::initializer_list<double> direction) = 0;


    /**
     * \brief Set the direction of the material's fibers.
     * The fibers direction is expressed by a matrix collecting the nodal fiber directions.
     * \param [in] fiber_direction The vector of the material's fibers direction.
     * \return [void]
     */
    virtual void SetLongFibers(const Eigen::MatrixXd &fiber_direction) = 0;


    /**
     * \brief Compute the material's electric diffusion tensor.
     * \todo Implement diffusion tensor computation for different types of tissue anisotropy.
     * \return [void]
     */
    virtual void ComputeDiffusionTensors() = 0;


    /**
     * \brief Compute the material's electric diffusion tensor for one-dimensional material in higher dimensional space.
     * The diffusion tensor is simply the diffusion coefficient of the material.
     * \return [void]
     */
    virtual void ComputeDiffusionTensors1D() = 0;


    /**
     * \brief Get the number of nodes belonging to the electric material.
     * \return [int] The number of nodes belonging to the electric material.
     */
    virtual int NodesNum() const = 0;


    /**
     * \brief Get the internal diffusivity of the material.
     * \return [const std::vector<double>&] The internal diffusivity values.
     */
    virtual const std::vector<double> & InterDiffusivity() const = 0;


    /**
     * \brief Get the internal diffusivity of the material with permission to modify it.
     * \return [std::vector<double>&] The internal diffusivity values.
     */
    virtual std::vector<double> & EditInterDiffusivity() = 0;


    /**
     * \brief Get the external diffusivity of the material.
     * \return [const std::vector<double>&] The external diffusivity values.
     */
    virtual const std::vector<double> & ExterDiffusivity() const = 0;


    /**
     * \brief Get the external diffusivity of the material with permission to modify it.
     * \return [std::vector<double>&] The external diffusivity values.
     */
    virtual std::vector<double> & EditExterDiffusivity() = 0;


    /**
     * \brief Get the transmembrane diffusivity of the material.
     * \return [const std::vector<double>&] The transmembrane diffusivity values.
     */
    virtual const std::vector<double> & TransDiffusivity() const = 0;


    /**
     * \brief Get the transmembrane diffusivity of the material with permission to modify it.
     * \return [std::vector<double>&] The transmembrane diffusivity values.
     */
    virtual std::vector<double> & EditTransDiffusivity() = 0;


    /**
     * \brief Get the transversal to longitudinal diffusivity ratio of the material.
     * \return [const std::vector<double>&] The transversal to longitudinal diffusivity ratio value.
     */
    virtual const std::vector<double> & TransversalRatio() const = 0;


    /**
     * \brief Get the capacitance of the material.
     * \return [const std::vector<double>&] The capcacitance value.
     */
    virtual const std::vector<double> & Capacitance() const = 0;


    /**
     * \brief Get the capacitance of the material for the node with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the capacitance value.
     * \return [double] The capcacitance value at the node with index id.
     */
    virtual double Capacitance(std::size_t id) const = 0;


    /**
     * \brief Get the fiber direction vector of the material.
     * \return [const Eigen::MatrixXd&] The fiber direction vector of the material.
     */
    virtual const Eigen::MatrixXd & LongFibers() const = 0;


    /**
     * \brief Get the internal diffusion tensor for each node of the material.
     * \return [const std::vector<Eigen::MatrixXd>&] The internal diffusion tensor of each node of the material.
     */
    virtual const std::vector<Eigen::MatrixXd> & InterDiffusionTensors() const = 0;


    /**
     * \brief Get the internal diffusion tensor fot the node of the material with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the internal diffusion tensor.
     * \return [const Eigen::MatrixXd&] The internal diffusion tensor of the node of the material with index id.
     */
    virtual const Eigen::MatrixXd & InterDiffusionTensors(std::size_t id) const = 0;


    /**
     * \brief Get the external diffusion tensor for each node of the material.
     * \return [const std::vector<Eigen::MatrixXd>&] The external diffusion tensor of each node of the material.
     */
    virtual const std::vector<Eigen::MatrixXd> & ExterDiffusionTensors() const = 0;


    /**
     * \brief Get the external diffusion tensor fot the node of the material with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the external diffusion tensor.
     * \return [const Eigen::MatrixXd&] The external diffusion tensor of the node of the material with index id.
     */
    virtual const Eigen::MatrixXd & ExterDiffusionTensors(std::size_t id) const = 0;


    /**
     * \brief Get the transmembrane diffusion tensor for each node of the material.
     * \return [const std::vector<Eigen::MatrixXd>&] The transmembrane diffusion tensor of each node of the material.
     */
    virtual const std::vector<Eigen::MatrixXd> & TransDiffusionTensors() const = 0;


    /**
     * \brief Get the transmembrane diffusion tensor fot the node of the material with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the transmembrane diffusion tensor.
     * \return [const Eigen::MatrixXd&] The transmembrane diffusion tensor of the node of the material with index id.
     */
    virtual const Eigen::MatrixXd & TransDiffusionTensors(std::size_t id) const = 0;

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_MATERIALS_ELECTRIC_BASIC_HPP_

// #include "ELECTRA/engine/materials/electric_basic.tpp"