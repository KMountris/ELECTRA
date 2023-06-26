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


#ifndef ELECTRA_MATERIALS_ELECTRIC_TRANSVERSAL_TPP_
#define ELECTRA_MATERIALS_ELECTRIC_TRANSVERSAL_TPP_


#include "ELECTRA/engine/materials/electric_transversal.hpp"


namespace ELECTRA {

template<short DIM>
ElectricTransversal<DIM>::ElectricTransversal() : nodes_num_(0), inter_diffusivity_(), exter_diffusivity_(), trans_diffusivity_(),
        transversal_ratio_(), capacitance_(), long_fibers_(), inter_diffusion_tensors_(), exter_diffusion_tensors_(), trans_diffusion_tensors_()
{}


template<short DIM>
ElectricTransversal<DIM>::~ElectricTransversal()
{}


template<short DIM>
void ElectricTransversal<DIM>::SetNodesNum(int nodes_num)
{
    // Check if given number of material nodes is positive.
    if (nodes_num < 0) {
        throw std::invalid_argument(Logger::Error("Could not set the number of nodes belonging to the electric material. A negative number of material nodes was given."));
    }

    this->nodes_num_ = nodes_num;
}


template<short DIM>
void ElectricTransversal<DIM>::SetInterDiffusivity(const std::vector<double> &diffusivity)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set internal diffusivity for electric material. The material's nodes number has not been set."));
    }

    // Set diffusivity. Either different for each node or the same.
    if (static_cast<int>(diffusivity.size()) == this->nodes_num_) {
        this->inter_diffusivity_ = diffusivity;
    }
    else if (diffusivity.size() == 1) {
        this->inter_diffusivity_.resize(this->nodes_num_, 0.);

        for (auto &nodal_val : this->inter_diffusivity_)
            nodal_val = diffusivity[0];
    } else {
        throw std::runtime_error(Logger::Error("Could not set internal diffusivity for electric material.\n"
                "The given diffusivity vector should contain either a single value or a value for each node of the material."));
    }

}


template<short DIM>
void ElectricTransversal<DIM>::SetExterDiffusivity(const std::vector<double> &diffusivity)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set external diffusivity for electric material. The material's nodes number has not been set."));
    }

    // Set diffusivity. Either different for each node or the same.
    if (static_cast<int>(diffusivity.size()) == this->nodes_num_) {
        this->exter_diffusivity_ = diffusivity;
    }
    else if (diffusivity.size() == 1) {
        this->exter_diffusivity_.resize(this->nodes_num_, 0.);

        for (auto &nodal_val : this->exter_diffusivity_)
            nodal_val = diffusivity[0];
    } else {
        throw std::runtime_error(Logger::Error("Could not set external diffusivity for electric material.\n"
                "The given diffusivity vector should contain either a single value or a value for each node of the material."));
    }

}


template<short DIM>
void ElectricTransversal<DIM>::SetTransDiffusivity(const std::vector<double> &diffusivity)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set transmembrane diffusivity for electric material. The material's nodes number has not been set."));
    }

    // Set diffusivity. Either different for each node or the same.
    if (static_cast<int>(diffusivity.size()) == this->nodes_num_) {
        this->trans_diffusivity_ = diffusivity;
    }
    else if (diffusivity.size() == 1) {
        this->trans_diffusivity_.resize(this->nodes_num_, 0.);

        for (auto &nodal_val : this->trans_diffusivity_)
            nodal_val = diffusivity[0];
    } else {
        throw std::runtime_error(Logger::Error("Could not set transmembrane diffusivity for electric material.\n"
                "The given diffusivity vector should contain either a single value or a value for each node of the material."));
    }

}


template<short DIM>
void ElectricTransversal<DIM>::SetTransversalRatio(const std::vector<double> &trans_ratio)
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
void ElectricTransversal<DIM>::SetCapacitance(const std::vector<double> &capacitance)
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
void ElectricTransversal<DIM>::SetLongFibers(std::initializer_list<double> direction)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set fiber direction for the electric material. The material's nodes number has not been set."));
    }

    // Check the size of the initializer list.
    if (direction.size() != DIM) {
        throw std::invalid_argument(Logger::Error("Could not set fiber direction for the electric material. Fibers direction initializer list should have the same dimensions with the ElectriMaterial."));
    }

    // Initialize the fiber direction matrix.
    this->long_fibers_ = Eigen::MatrixXd::Zero(this->nodes_num_, DIM);

    auto dir_norm = double{0.};
    for (const auto &val : direction) {
        dir_norm += val*val;
    }
    dir_norm = std::sqrt(dir_norm);

    // Assign the fiber direction components.
    for (int i = 0; i != this->nodes_num_; ++i) {
        for (const auto &val : direction) {
            auto j = &val - &direction.begin()[0];
            this->long_fibers_.coeffRef(i, j) = val/dir_norm;
        }
    }

}


template<short DIM>
void ElectricTransversal<DIM>::SetLongFibers(const Eigen::MatrixXd &fiber_direction)
{
    // Check if the nodes of the electric material have been set.
    if (this->nodes_num_ == 0) {
        throw std::runtime_error(Logger::Error("Could not set fiber direction for the electric material. The material's nodes number has not been set."));
    }

    // Check the size of the initializer list.
    if (fiber_direction.rows() != this->nodes_num_ && fiber_direction.cols() != DIM) {
        throw std::invalid_argument(Logger::Error("Could not set fiber direction for the electric material. Fibers direction matrix does not have entries for all the material's nodes."));
    }

    // Set the fiber direction of the material.
    this->long_fibers_ = Eigen::MatrixXd::Zero(fiber_direction.rows(), fiber_direction.cols());
    for (auto i = 0; i != fiber_direction.rows(); ++i) {
        this->long_fibers_.row(i) = fiber_direction.row(i) / fiber_direction.row(i).norm();
    }
}


template<short DIM>
void ElectricTransversal<DIM>::ComputeDiffusionTensors()
{
    // Compute internal diffusion tensors.
    if (this->inter_diffusivity_.size() != 0) {
        this->inter_diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(DIM, DIM));
        for (int i = 0; i != this->nodes_num_; ++i) {
            this->inter_diffusion_tensors_[i] = this->inter_diffusivity_[i]*(1. - this->transversal_ratio_[i])*(this->long_fibers_.row(i).transpose()*this->long_fibers_.row(i)) +
                    this->inter_diffusivity_[i]*this->transversal_ratio_[i]*Eigen::MatrixXd::Identity(DIM, DIM);
        }
    }

    // Compute external diffusion tensors.
    if (this->exter_diffusivity_.size() != 0) {
        this->exter_diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(DIM, DIM));
        for (int i = 0; i != this->nodes_num_; ++i) {
            this->exter_diffusion_tensors_[i] = this->exter_diffusivity_[i]*(1. - this->transversal_ratio_[i])*(this->long_fibers_.row(i).transpose()*this->long_fibers_.row(i)) +
                    this->exter_diffusivity_[i]*this->transversal_ratio_[i]*Eigen::MatrixXd::Identity(DIM, DIM);
        }
    }

    // Compute transmembrane diffusion tensors.
    if (this->trans_diffusivity_.size() != 0) {
        this->trans_diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(DIM, DIM));
        for (int i = 0; i != this->nodes_num_; ++i) {
            this->trans_diffusion_tensors_[i] = this->trans_diffusivity_[i]*(1. - this->transversal_ratio_[i])*(this->long_fibers_.row(i).transpose()*this->long_fibers_.row(i)) +
                    this->trans_diffusivity_[i]*this->transversal_ratio_[i]*Eigen::MatrixXd::Identity(DIM, DIM);
        }
    }
}


template<short DIM>
void ElectricTransversal<DIM>::ComputeDiffusionTensors1D()
{
    // Compute internal diffusion tensors.
    if (this->inter_diffusivity_.size() != 0) {
        this->inter_diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(1, 1));
        for (int i = 0; i != this->nodes_num_; ++i) {
            this->inter_diffusion_tensors_[i].setConstant(this->inter_diffusivity_[i]);
        }
    }

    // Compute external diffusion tensors.
    if (this->exter_diffusivity_.size() != 0) {
        this->exter_diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(1, 1));
        for (int i = 0; i != this->nodes_num_; ++i) {
            this->exter_diffusion_tensors_[i].setConstant(this->exter_diffusivity_[i]);
        }
    }

    // Compute transmembrane diffusion tensors.
    if (this->trans_diffusivity_.size() != 0) {
        this->trans_diffusion_tensors_.resize(this->nodes_num_, Eigen::MatrixXd::Zero(1, 1));
        for (int i = 0; i != this->nodes_num_; ++i) {
            this->trans_diffusion_tensors_[i].setConstant(this->trans_diffusivity_[i]);
        }
    }
}



} // End of namespace ELECTRA

#endif //ELECTRA_MATERIALS_ELECTRIC_TRANSVERSAL_TPP_