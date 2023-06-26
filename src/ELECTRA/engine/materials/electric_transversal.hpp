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
   \file electric_transversal.hpp
   \brief ElectricTransversal class header file.
   \author Konstantinos A. Mountris
   \date 02/11/2021
*/

#ifndef ELECTRA_MATERIALS_ELECTRIC_TRANSVERSAL_HPP_
#define ELECTRA_MATERIALS_ELECTRIC_TRANSVERSAL_HPP_


#include "ELECTRA/engine/materials/electric_basic.hpp"
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
 * \class ElectricTransversal
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting an anisotropic electric material with transversal anisotropy.
 * \tparam DIM The dimensions of the material's domain.
 */

template<short DIM>
class ElectricTransversal : public ElectricBasic<DIM>
{

private:

    int nodes_num_;                                     /**< The number of the material's nodes */

    std::vector<double> inter_diffusivity_;          /**< The material's internal diffusivity per domain unit (m, m2, m3, ...). If the material is anisotropic, it represents the longitudinal diffusivity. SI units: m/s */

    std::vector<double> exter_diffusivity_;          /**< The material's external diffusivity per domain unit (m, m2, m3, ...). If the material is anisotropic, it represents the longitudinal diffusivity. SI units: m/s */

    std::vector<double> trans_diffusivity_;     /**< The material's transmembrane diffusivity per domain unit (m, m2, m3, ...). If the material is anisotropic, it represents the longitudinal diffusivity. SI units: m/s */

    std::vector<double> transversal_ratio_;             /**< Transversal to longitudinal diffusivity ratio. It is expressed with a decimal value in the range [0.0 - 1.0] */

    std::vector<double> capacitance_;                   /**< The material's capacitance per domain unit (m, m2, m3, ...). SI units: F/domain_unit */

    Eigen::MatrixXd long_fibers_;                   /**< The material's fiber direction vector. Each axis component is expressed with a decimal value in the range [0.0 - 1.0] */

    std::vector<Eigen::MatrixXd> inter_diffusion_tensors_;    /**< The material's internal diffusion tensors */

    std::vector<Eigen::MatrixXd> exter_diffusion_tensors_;    /**< The material's external diffusion tensors */

    std::vector<Eigen::MatrixXd> trans_diffusion_tensors_;    /**< The material's transmembrane diffusion tensors */


public:

    /**
     * \brief ElectricTransversal default constructor.
     */
    inline ElectricTransversal();


    /**
     * \brief ElectricTransversal default destructor.
     */
    inline virtual ~ElectricTransversal();


    /**
     * \brief Set the number of the electric material's nodes.
     * \param [in] nodes_num The number of the nodes that belong to the electric material.
     * \return [void]
     */
    inline void SetNodesNum(int nodes_num);


    /**
     * \brief Set the internal diffusivity of the material.
     * \param [in] diffusivity The nodal internal diffusivity values.
     * \return [void]
     */
    inline void SetInterDiffusivity(const std::vector<double> &diffusivity);


    /**
     * \brief Set the external diffusivity of the material.
     * \param [in] diffusivity The nodal external diffusivity values.
     * \return [void]
     */
    inline void SetExterDiffusivity(const std::vector<double> &diffusivity);


    /**
     * \brief Set the transmembrane diffusivity of the material.
     * \param [in] diffusivity The nodal transmembrane diffusivity values.
     * \return [void]
     */
    inline void SetTransDiffusivity(const std::vector<double> &diffusivity);


    /**
     * \brief Set the transversal to longitudinal diffusivity ratio of the material.
     * \param [in] trans_ratio The value of the material's transversal-longitudinal diffusivity ratio.
     * \return [void]
     */
    inline void SetTransversalRatio(const std::vector<double> &trans_ratio);


    /**
     * \brief Set the capacitance of the material.
     * \param [in] capacitance The value of the material's capacitance.
     * \return [void]
     */
    inline void SetCapacitance(const std::vector<double> &capacitance);


    /**
     * \brief Set the direction of the material's fibers.
     *
     * The fibers direction is expressed by a initializer list collecting the axial components
     * of the fiber direction. A value of 0.0 on a specific axis means no fiber direction
     * on this axis and a value of 1.0 means that the fibers direction is parallel to this
     * axis. The sum of the axial components should be equal to 1.0.
     *
     * \param [in] direction The initializer list of the material's fibers direction.
     * \return [void]
     */
    inline void SetLongFibers(std::initializer_list<double> direction);


    /**
     * \brief Set the direction of the material's fibers.
     *
     * The fibers direction is expressed by a matrix collecting the nodal fiber directions.
     *
     * \param [in] fiber_direction The vector of the material's fibers direction.
     * \return [void]
     */
    inline void SetLongFibers(const Eigen::MatrixXd &fiber_direction);


    /**
     * \brief Compute the material's electric diffusion tensor.
     * \todo Implement diffusion tensor computation for different types of tissue anisotropy.
     * \return [void]
     */
    inline void ComputeDiffusionTensors();


    /**
     * \brief Compute the material's electric diffusion tensor for one-dimensional material in higher dimensional space.
     * The diffusion tensor is simply the diffusion coefficient of the material.
     * \return [void]
     */
    inline void ComputeDiffusionTensors1D();


    /**
     * \brief Get the number of nodes belonging to the electric material.
     * \return [int] The number of nodes belonging to the electric material.
     */
    inline int NodesNum() const { return this->nodes_num_; }


    /**
     * \brief Get the internal diffusivity of the material.
     * \return [const std::vector<double>&] The internal diffusivity values.
     */
    inline const std::vector<double> & InterDiffusivity() const { return this->inter_diffusivity_; }


    /**
     * \brief Get the internal diffusivity of the material with permission to modify it.
     * \return [std::vector<double>&] The internal diffusivity values.
     */
    inline std::vector<double> & EditInterDiffusivity() { return this->inter_diffusivity_; }


    /**
     * \brief Get the external diffusivity of the material.
     * \return [const std::vector<double>&] The external diffusivity values.
     */
    inline const std::vector<double> & ExterDiffusivity() const { return this->exter_diffusivity_; }


    /**
     * \brief Get the external diffusivity of the material with permission to modify it.
     * \return [std::vector<double>&] The external diffusivity values.
     */
    inline std::vector<double> & EditExterDiffusivity() { return this->exter_diffusivity_; }


    /**
     * \brief Get the transmembrane diffusivity of the material.
     * \return [const std::vector<double>&] The transmembrane diffusivity values.
     */
    inline const std::vector<double> & TransDiffusivity() const { return this->trans_diffusivity_; }


    /**
     * \brief Get the transmembrane diffusivity of the material with permission to modify it.
     * \return [std::vector<double>&] The transmembrane diffusivity values.
     */
    inline std::vector<double> & EditTransDiffusivity() { return this->trans_diffusivity_; }


    /**
     * \brief Get the transversal to longitudinal diffusivity ratio of the material.
     * \return [const std::vector<double>&] The transversal to longitudinal diffusivity ratio value.
     */
    inline const std::vector<double> & TransversalRatio() const { return this->transversal_ratio_; }


    /**
     * \brief Get the capacitance of the material.
     * \return [const std::vector<double>&] The capcacitance value.
     */
    inline const std::vector<double> & Capacitance() const { return this->capacitance_; }


    /**
     * \brief Get the capacitance of the material for the node with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the capacitance value.
     * \return [double] The capcacitance value at the node with index id.
     */
    inline double Capacitance(std::size_t id) const { return this->capacitance_[id]; }


    /**
     * \brief Get the fiber direction vector of the material.
     * \return [const Eigen::MatrixXd&] The fiber direction vector of the material.
     */
    inline const Eigen::MatrixXd & LongFibers() const { return this->long_fibers_; }


    /**
     * \brief Get the internal diffusion tensor for each node of the material.
     * \return [const std::vector<Eigen::MatrixXd>&] The internal diffusion tensor of each node of the material.
     */
    inline const std::vector<Eigen::MatrixXd> & InterDiffusionTensors() const { return this->inter_diffusion_tensors_; }


    /**
     * \brief Get the internal diffusion tensor fot the node of the material with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the internal diffusion tensor.
     * \return [const Eigen::MatrixXd&] The internal diffusion tensor of the node of the material with index id.
     */
    inline const Eigen::MatrixXd & InterDiffusionTensors(std::size_t id) const { return this->inter_diffusion_tensors_[id]; }


    /**
     * \brief Get the external diffusion tensor for each node of the material.
     * \return [const std::vector<Eigen::MatrixXd>&] The external diffusion tensor of each node of the material.
     */
    inline const std::vector<Eigen::MatrixXd> & ExterDiffusionTensors() const { return this->exter_diffusion_tensors_; }


    /**
     * \brief Get the external diffusion tensor fot the node of the material with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the external diffusion tensor.
     * \return [const Eigen::MatrixXd&] The external diffusion tensor of the node of the material with index id.
     */
    inline const Eigen::MatrixXd & ExterDiffusionTensors(std::size_t id) const { return this->exter_diffusion_tensors_[id]; }


    /**
     * \brief Get the transmembrane diffusion tensor for each node of the material.
     * \return [const std::vector<Eigen::MatrixXd>&] The transmembrane diffusion tensor of each node of the material.
     */
    inline const std::vector<Eigen::MatrixXd> & TransDiffusionTensors() const { return this->trans_diffusion_tensors_; }


    /**
     * \brief Get the transmembrane diffusion tensor fot the node of the material with index id.
     * Fast access, no range check.
     * \param [in] id The index of the node for which we will get the transmembrane diffusion tensor.
     * \return [const Eigen::MatrixXd&] The transmembrane diffusion tensor of the node of the material with index id.
     */
    inline const Eigen::MatrixXd & TransDiffusionTensors(std::size_t id) const { return this->trans_diffusion_tensors_[id]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif //ELECTRA_MATERIALS_ELECTRIC_TRANSVERSAL_HPP_

#include "ELECTRA/engine/materials/electric_transversal.tpp"