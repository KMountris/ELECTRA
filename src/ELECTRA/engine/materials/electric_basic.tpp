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


#ifndef ELECTRA_MATERIALS_ELECTRIC_MATERIAL_TPP_
#define ELECTRA_MATERIALS_ELECTRIC_MATERIAL_TPP_


#include "ELECTRA/engine/materials/electric_basic.hpp"


namespace ELECTRA {

template<short DIM>
ElectricBasic<DIM>::ElectricBasic() : nodes_num_(0), diffusivity_(), transversal_ratio_(), capacitance_(),
        fiber_direction_(), diffusion_tensors_()
{}


template<short DIM>
ElectricBasic<DIM>::~ElectricBasic()
{}


template<short DIM>
void ElectricBasic<DIM>::SetNodesNum(int nodes_num)
{
    // Check if given number of material nodes is positive.
    if (nodes_num < 0) {
        throw std::invalid_argument(Logger::Error("Could not set the number of nodes belonging to the electric material. A negative number of material nodes was given."));
    }

    this->nodes_num_ = nodes_num;
}


template<short DIM>
void ElectricBasic<DIM>::SetDiffusivity(const std::vector<double> &diffusivity)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set diffusivity for electric material. The material's nodes number has not been set."));
    }

    // Set diffusivity. Either different for each node or the same.
    if (static_cast<int>(diffusivity.size()) == this->nodes_num_) {
        this->diffusivity_ = diffusivity;
    }
    else if (diffusivity.size() == 1) {
        this->diffusivity_.resize(this->nodes_num_, 0.);

        for (auto &nodal_val : this->diffusivity_)
            nodal_val = diffusivity[0];
    } else {
        throw std::runtime_error(Logger::Error("Could not set diffusivity for electric material.\n"
                                               "The given diffusivity vector should contain either a single value or a value for each node of the material."));
    }


}


template<short DIM>
void ElectricBasic<DIM>::SetTransversalRatio(const std::vector<double> &trans_ratio)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set transversal ratio for electric material. The material's nodes number has not been set."));
    }

    // Set transversal ratio. Either different for each node or the same.
    if (static_cast<int>(trans_ratio.size()) == this->nodes_num_) {
        this->transversal_ratio_ = trans_ratio;
    }
    else if (trans_ratio.size() == 1) {
        this->transversal_ratio_.resize(this->nodes_num_, 0.);

        for (auto &nodal_val : this->transversal_ratio_)
            nodal_val = trans_ratio[0];   
    } else {
        throw std::runtime_error(Logger::Error("Could not set transversal ratio for electric material.\n"
                                               "The given transversal ratio vector should contain either a single value or a value for each node of the material."));
    }

}


template<short DIM>
void ElectricBasic<DIM>::SetCapacitance(const std::vector<double> &capacitance)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set capacitance for electric material. The material's nodes number has not been set."));
    }

    // Set transversal ratio. Either different for each node or the same.
    if (static_cast<int>(capacitance.size()) == this->nodes_num_) {
        this->capacitance_ = capacitance;
    }
    else if (capacitance.size() == 1) {
        this->capacitance_.resize(this->nodes_num_, 0.);

        for (auto &nodal_val : this->capacitance_)
            nodal_val = capacitance[0];   
    } else {
        throw std::runtime_error(Logger::Error("Could not set capacitance for electric material.\n"
                                               "The given capacitance vector should contain either a single value or a value for each node of the material."));
    }

}


template<short DIM>
void ElectricBasic<DIM>::SetFiberDirection(std::initializer_list<double> direction)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set fiber direction for the electric material. The material's nodes number has not been set."));
    }

    // Check the size of the initializer list.
    if (direction.size() != DIM) {
        throw std::invalid_argument(Logger::Error("Could not set fiber direction for the electric material. Fibers direction initializer list should have the same dimensions with the ElectriMaterial."));
    }

    // Check that initializer list sums to unity.
    double dir_sum = 0.;
    for(std::size_t i = 0; i != direction.size(); ++i) {
        dir_sum += std::abs(direction.begin()[i]);
    }
    if (dir_sum < 0.9999 || dir_sum > 1.0001) {
        throw std::invalid_argument(Logger::Error("Could not set fiber direction for the electric material. The fibers direction components should sum to 1."));
    }

    // Initialize the fiber direction matrix.
    this->fiber_direction_ = Eigen::MatrixXd::Zero(this->nodes_num_, DIM);

    // Assign the fiber direction components.
    for (int i = 0; i != this->nodes_num_; ++i) {
        for (const auto &val : direction) {
            auto j = &val - &direction.begin()[0];
            this->fiber_direction_.coeffRef(i, j) = val;
        }
    }

}


template<short DIM>
void ElectricBasic<DIM>::SetFiberDirection(const Eigen::MatrixXd &fiber_direction)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set fiber direction for the electric material. The material's nodes number has not been set."));
    }

    // Check the size of the initializer list.
    if (fiber_direction.rows() != this->nodes_num_ && fiber_direction_.cols() != DIM) {
        throw std::invalid_argument(Logger::Error("Could not set fiber direction for the electric material. Fibers direction matrix does not have entries for all the material's nodes."));
    }

    // Check that initializer list sums to unity.
    double dir_sum = 0.;
    for(int i = 0; i != fiber_direction_.rows(); ++i) {
        dir_sum = 0.;
        for (int j = 0; j != fiber_direction_.cols(); ++j) {
            dir_sum += std::abs(fiber_direction(i,j));
        }
        if (dir_sum < 0.9999 || dir_sum > 1.0001) {
            throw std::invalid_argument(Logger::Error("Could not set fiber direction for the electric material. The fibers direction components for each node should sum to 1."));
        }
    }

    // Set the fiber direction of the material.
    this->fiber_direction_ = fiber_direction;
}


template<short DIM>
void ElectricBasic<DIM>::ComputeDiffusionTensors()
{
    // Initialize the diffusion tensors for the nodes of the material.
    this->diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(DIM, DIM));

    // Compute orthotropic diffusion tensor for each material node.
    for (int i = 0; i != this->nodes_num_; ++i) {
        this->diffusion_tensors_[i] = this->diffusivity_[i]*(1. - this->transversal_ratio_[i])*(this->fiber_direction_.row(i).transpose()*this->fiber_direction_.row(i)) +
                                      this->diffusivity_[i]*this->transversal_ratio_[i]*Eigen::MatrixXd::Identity(DIM, DIM);
    }
}


template<short DIM>
void ElectricBasic<DIM>::ComputeDiffusionTensors1D()
{
    // Initialize the 1D diffusion tensors for the nodes of the material.
    this->diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(1, 1));

    // Diffusion tensor is simply the diffusion coeeficient.
    for (int i = 0; i != this->nodes_num_; ++i) {
        this->diffusion_tensors_[i].setConstant(this->diffusivity_[i]);
    }
}



} // End of namespace ELECTRA

#endif //ELECTRA_MATERIALS_ELECTRIC_MATERIAL_TPP_