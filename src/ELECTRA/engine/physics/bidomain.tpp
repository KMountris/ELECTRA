/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#ifndef ELECTRA_ENGINE_PHYSICS_BIDOMAIN_TPP_
#define ELECTRA_ENGINE_PHYSICS_BIDOMAIN_TPP_


#include "ELECTRA/engine/physics/bidomain.hpp"


namespace ELECTRA {


template<short DIM, short CELL_NODES>
Bidomain<DIM, CELL_NODES>::Bidomain() : vout_(), stiff_int_mat_(), mass_vec_(),
        thread_loop_manager_(), threads_number_(0.), thread_mutex_()
{
    // Get the number of parallel threads.
    const std::size_t available_threads = std::thread::hardware_concurrency()-1;
    this->threads_number_ = std::max(available_threads, 1ul);
}


template<short DIM, short CELL_NODES>
Bidomain<DIM, CELL_NODES>::~Bidomain()
{}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::AssembleMatrices(const IMP::Mesh<DIM, CELL_NODES> &mesh, const std::shared_ptr<ElectricBasic<DIM>> &material)
{
    if (static_cast<int>(material->InterDiffusionTensors().size()) != mesh.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for monodomain model. Material internal diffusion tensors are not available for all mesh nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }
    if (static_cast<int>(material->ExterDiffusionTensors().size()) != mesh.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for monodomain model. Material external diffusion tensors are not available for all mesh nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Clear the stiffness matrix and the lumped mass vector.
    auto n = material->NodesNum() + this->ConductSystem().NodesNum();
    this->stiff_int_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n,n);
    this->stiff_bulk_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n,n);
    this->mass_vec_.setZero(n);

    // Compute the cell dimensions.
    constexpr short CELL_DIM = (DIM > 1 && CELL_NODES == 2) ? 1 : DIM;

    // Create a finite element of the given type.
    auto elem = CLOUDEA::FemFactory<CELL_DIM>::Create(mesh.CellsShape());

    // Set the quadrature of the finite element.
    switch (mesh.CellsShape()) {
        case IMP::CellShape::edge :
            elem->SetQuadrature({2});
        break;
        case IMP::CellShape::tri :
            elem->SetQuadrature({3});
        break;
        case IMP::CellShape::quad :
            elem->SetQuadrature({2,2});
        break;
        case IMP::CellShape::tet :
            elem->SetQuadrature({4});
        break;
        case IMP::CellShape::hex :
            elem->SetQuadrature({2,2,2});
        break;
        default:
            std::string error_msg = "Unknown finite element type requested during assembling the monodomain model matrices.";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

    // Initialize triplets for assembly of internal and bulk stiffness matrices.
    typedef Eigen::Triplet<double> EigTripD;
    auto stiff_int_assemble_triplet = std::vector<EigTripD>{};
    auto stiff_bulk_assemble_triplet = std::vector<EigTripD>{};
    stiff_int_assemble_triplet.reserve(n * CELL_NODES);
    stiff_bulk_assemble_triplet.reserve(n * CELL_NODES);

    // Initialize container to store the cell node coordinates.
    auto nodes = std::vector<IMP::Vec<CELL_DIM, double>>(CELL_NODES, IMP::Vec<CELL_DIM, double>());

    // Iterate over the mesh cells.
    auto stiff_int_cell = Eigen::MatrixXd(CELL_NODES, CELL_NODES);     // Internal stiffness matrix of a cell.
    auto stiff_bulk_cell = Eigen::MatrixXd(CELL_NODES, CELL_NODES);    // Bulk stiffness matrix of a cell (internal+external).
    auto mass_cell = Eigen::MatrixXd(CELL_NODES, CELL_NODES);          // Mass matrix of a cell.
    auto mean_sigma_int = Eigen::MatrixXd(CELL_DIM, CELL_DIM);         // Mean internal diffusion tensor of a cell.
    auto mean_sigma_bulk = Eigen::MatrixXd(CELL_DIM, CELL_DIM);         // Mean external diffusion tensor of a cell.
    auto mean_capacitance = 0.;
    for (const auto &cell : mesh.Cells()) {

        // Reset mean internal, external diffusion tensors and mean capacitance values.
        mean_sigma_int.setZero();
        mean_sigma_bulk.setZero();
        mean_capacitance = 0.;

        // Compute mean cell diffusion tensor and capacitance.
        for (const auto &nid : cell.Connectivity()) {
            mean_sigma_int += material->InterDiffusionTensors(nid);
            mean_sigma_bulk += (material->InterDiffusionTensors(nid)+material->ExterDiffusionTensors(nid));
            mean_capacitance += material->Capacitance()[nid];
        }
        mean_sigma_int /= CELL_NODES;
        mean_sigma_bulk /= CELL_NODES;
        mean_capacitance /= CELL_NODES;

        // Get the cell nodes coordinates.
        if (CELL_DIM != DIM) {
            // Set nodes coordinates from distance for 1D nodes in nD space.
            nodes[0][0] = 0.;
            nodes[1][0] = std::sqrt(mesh.Nodes(cell.N(1)).Distance2( mesh.Nodes(cell.N(0))) );
        } else {
            // Get coordinates from mesh nodes container for nD nodes in nD space.
            for (const auto &nid : cell.Connectivity()) {
                auto i = &nid - &cell.Connectivity()[0];
                for (short d = 0; d != CELL_DIM; ++d) { nodes[i][d] = mesh.Nodes(nid)[d]; }
            }
        }

        // Compute the jacobian, shape function and derivatives for the current cell.
        elem->ComputeDerivsNatural();  elem->ComputeJacobians(nodes);
        elem->ComputeShapeFunctions(); elem->ComputeDerivs();

        // Construct the local stiffness and mass over the element's quadrature.
        stiff_int_cell.setZero();
        stiff_bulk_cell.setZero();
        mass_cell.setZero();
        for (const auto &qweight : elem->Quadrature().Weights()) {
            auto id = std::distance(&elem->Quadrature().Weights()[0], &qweight);
            stiff_int_cell += qweight*elem->DetJacobians()[id] *
                    (elem->Derivs()[id]*mean_sigma_int*elem->Derivs()[id].transpose());
            stiff_bulk_cell += qweight*elem->DetJacobians()[id] *
                    (elem->Derivs()[id]*mean_sigma_bulk*elem->Derivs()[id].transpose());
            mass_cell += mean_capacitance*qweight *
                    elem->DetJacobians()[id]*elem->ShapeFunctions()[id]*elem->ShapeFunctions()[id].transpose();
        }

        // Collect the cell stiffnesses and mass contribution.
        for (short local_i = 0; local_i != CELL_NODES; ++local_i) {
            auto global_i = cell.N(local_i);

            for (short local_j = 0; local_j != CELL_NODES; ++local_j) {
                auto global_j = cell.N(local_j);

                stiff_int_assemble_triplet.emplace_back( EigTripD(global_i, global_j, stiff_int_cell(local_i, local_j)) );
                stiff_bulk_assemble_triplet.emplace_back( EigTripD(global_i, global_j, stiff_bulk_cell(local_i, local_j)) );
                this->mass_vec_[global_i] += mass_cell(local_i, local_j);
            }
        } // End Collect stiffness contribution.

    } // End of Iterate over the mesh cells.

    // Set cardiac stiffness matrix.
    this->stiff_int_mat_.setFromTriplets(std::begin(stiff_int_assemble_triplet), std::end(stiff_int_assemble_triplet));
    this->stiff_bulk_mat_.setFromTriplets(std::begin(stiff_bulk_assemble_triplet), std::end(stiff_bulk_assemble_triplet));

    stiff_int_assemble_triplet.clear(); stiff_int_assemble_triplet.shrink_to_fit();
    stiff_bulk_assemble_triplet.clear(); stiff_bulk_assemble_triplet.shrink_to_fit();

    // Compute the stable time step for the diffusion term of the monodomain model.
    // this->ComputeCriticalTimeStep();

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::AssembleMatrices(const IMP::Grid<DIM, CELL_NODES> &grid, const std::shared_ptr<ElectricBasic<DIM>> &material,
        const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx, const IMP::NodeSet &neumann_set)
{

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::AssembleMatrices(const IMP::Voronoi<DIM> &voro,
        const std::shared_ptr<ElectricBasic<DIM>> &material, const CLOUDEA::Fpm<DIM> &fpm)
{
    if (static_cast<int>(fpm.Phi().size()) != voro.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for bidomain model. Fpm approximation is not available for all voronoi tesselation nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    if (static_cast<int>(material->InterDiffusionTensors().size()) != voro.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for bidomain model. Material internal diffusion tensors are not available for all voronoi tesselation nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    if (static_cast<int>(material->ExterDiffusionTensors().size()) != voro.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for bidomain model. Material external diffusion tensors are not available for all voronoi tesselation nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    typedef Eigen::Triplet<double> EigTripD;

    // Clear the stiffness matrix and the lumped mass vector.
    int n = material->NodesNum() + this->ConductSystem().NodesNum();
    this->stiff_int_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n,n);
    this->stiff_bulk_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n,n);
    this->mass_vec_.setZero(n);

    // Compute mean diffusion value.
    Eigen::VectorXd cell_volumes(n);
    double mean_int_kappa = 0., mean_ext_kappa = 0., total_volume = 0.;
    for (const auto &cell : voro.Cells()) {
        auto id = &cell - &voro.Cells()[0];
        cell_volumes.coeffRef(id) = cell.Measure(voro.Points(), voro.Facets());
        mean_int_kappa += std::pow(material->InterDiffusionTensors(id).diagonal().prod(),1./DIM)*cell_volumes.coeff(id);
        mean_ext_kappa += std::pow(material->ExterDiffusionTensors(id).diagonal().prod(),1./DIM)*cell_volumes.coeff(id);
        total_volume += cell_volumes.coeff(id);
    }
    double int_penal = fpm.Penalty() * (mean_int_kappa/total_volume);
    double bulk_penal = fpm.Penalty() * ((mean_int_kappa+mean_ext_kappa)/total_volume);

    // Initialize triplets for assembly.
    std::vector<EigTripD> stiff_int_triplet, stiff_bulk_triplet, mass_triplet;
    stiff_int_triplet.reserve(n * CELL_NODES);
    stiff_bulk_triplet.reserve(n * CELL_NODES);
    mass_triplet.reserve(n * CELL_NODES);

    // Assemble stiffness and mass matrices.
    Eigen::MatrixXd stiff_int_cell, stiff_bulk_cell, mass_cell;
    for (const auto &cell : voro.Cells()) {
        auto id = &cell - &voro.Cells()[0];

        // Compute stiffness and mass matrices of the cell.
        stiff_int_cell = fpm.PhiGrad(id)*material->InterDiffusionTensors(id)*fpm.PhiGrad(id).transpose()*cell_volumes.coeff(id);
        stiff_bulk_cell = fpm.PhiGrad(id)*(material->InterDiffusionTensors(id)+material->ExterDiffusionTensors(id))*fpm.PhiGrad(id).transpose()*cell_volumes.coeff(id);
        mass_cell = fpm.Phi(id)*fpm.Phi(id).transpose()*cell_volumes.coeff(id)*material->Capacitance(id);

        // Assemble matrices.
        for (const auto &gi : fpm.Support().InfluenceNodeIds(id)) {
            auto li = &gi - &fpm.Support().InfluenceNodeIds(id)[0];

            for (const auto &gj : fpm.Support().InfluenceNodeIds(id)) {
                auto lj = &gj - &fpm.Support().InfluenceNodeIds(id)[0];
                stiff_int_triplet.emplace_back( EigTripD(gi, gj, stiff_int_cell(li, lj)) );
                stiff_bulk_triplet.emplace_back( EigTripD(gi, gj, stiff_bulk_cell(li, lj)) );
                mass_triplet.emplace_back( EigTripD(gi, gj, mass_cell(li, lj)) );
            }
        }
    } // End of Assemble stiffness and mass matrices.

    // Assemble mass matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> mass_mat(n,n);
    mass_mat.setFromTriplets(std::begin(mass_triplet), std::end(mass_triplet));
    mass_triplet.clear(); mass_triplet.shrink_to_fit();

    // Convert to lumped mass vector.
    for (int i = 0; i != mass_mat.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator mass(mass_mat, i); mass; ++mass) {
            this->mass_vec_.coeffRef(i) = this->mass_vec_.coeff(i) + mass.value();
        }
    }
    mass_mat.setZero(); mass_mat.data().squeeze();

    // Set cardiac stiffness matrices.
    this->stiff_int_mat_.setFromTriplets(std::begin(stiff_int_triplet), std::end(stiff_int_triplet));
    stiff_int_triplet.clear(); stiff_int_triplet.shrink_to_fit();

    this->stiff_bulk_mat_.setFromTriplets(std::begin(stiff_bulk_triplet), std::end(stiff_bulk_triplet));
    stiff_bulk_triplet.clear(); stiff_bulk_triplet.shrink_to_fit();

    // Apply numerical correction at internal facets.
    switch (DIM) {
        case 2 : this->ApplyFpmCorrection2D(voro, material, fpm, int_penal, bulk_penal); break;

        case 3 : this->ApplyFpmCorrection3D(voro, material, fpm, int_penal, bulk_penal); break;

        default:
            std::string error_msg = "Could not apply numerical correction for FPM assembly. Domain dimensions are not supported.";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::ApplyFpmCorrection2D(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
        const CLOUDEA::Fpm<DIM> &fpm, double int_penalty, double bulk_penalty)
{
    typedef Eigen::Triplet<double> EigTripD;

    // Estimate facets division number for memory handling.
    int divs = 1;           // facets division.
    int rest = 0;           // resting facets number for not exact division.
    int chunk = 50000;     // division chunk size.
    int facets_num = static_cast<int>(voro.Facets().size());
    if (facets_num > chunk) {
        divs = facets_num / chunk;
        rest = facets_num - divs*chunk;
    } else { chunk = facets_num; }

    // Estimate correction triplet entries number.
    int entries_num = chunk*fpm.Support().MaxInfluenceNodesNum()*fpm.Support().MaxInfluenceNodesNum();

    // Initialize indices of cells sharing the facet.
    int n1 = 0, n2 = 0;

    // Initialize neighbor nodes of the nodes inside the facet-sharing cells.
    std::vector<int> neigh1, neigh2;

    // Initialize normals on facet and shape function.
    Eigen::VectorXd ne1, ne2, N1, N2;

    // Initialize shape function derivatives, diffusion tensors, and correction matrices.
    Eigen::MatrixXd B1, B2, lambda1, lambda2;
    Eigen::MatrixXd Cint11, Cint12, Cint21, Cint22;
    Eigen::MatrixXd Cbulk11, Cbulk12, Cbulk21, Cbulk22;

    // Initialize correction triplets vector and correction matrix.
    std::vector<EigTripD> corr_int_triplets, corr_bulk_triplets;
    corr_int_triplets.reserve(entries_num); corr_bulk_triplets.reserve(entries_num);
    auto corr_int_mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(this->stiff_int_mat_.rows(), this->stiff_int_mat_.cols());
    auto corr_bulk_mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(this->stiff_bulk_mat_.rows(), this->stiff_bulk_mat_.cols());

    // Iterate facets in chuncks.
    for (int i = 0; i != divs; ++i) {
        for (int fid = i*chunk; fid != (i+1)*chunk; ++fid) {
            // Process only internal facets.
            if (!voro.Facets(fid).IsFree()) {

                // Compute the facet length.
                double L = voro.Facets(fid).Length(voro.Points());

                // Indices of cells sharing the facet.
                n1 = voro.Facets(fid).ParentCellId();
                n2 = voro.Facets(fid).NeighCellId();

                // Compute the facet outward normal vector.
                IMP::Vec<DIM, double> normal = voro.Facets(fid).Normal(voro.Points());
                if (normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0)  normal = -normal;

                // Compute opposite vector to normal.
                ne1 = normal.CopyToEigen();
                ne2 = -ne1;

                // Boundary dependent parameter with unit of length
                double he = std::sqrt(voro.Nodes(n2).Distance2(voro.Nodes(n1)));

                // Get fpm gradients for the cells.
                B1 = fpm.PhiGrad(n1);
                B2 = fpm.PhiGrad(n2);

                // Get neighbor nodes to the nodes enclosed by the cells.
                neigh1 = fpm.Support().InfluenceNodeIds(n1);
                neigh2 = fpm.Support().InfluenceNodeIds(n2);

                // Compute integration point at the center of the facet.
                IMP::Vec<DIM, double> xt = voro.Facets(fid).Centroid(voro.Points());

                // Compute shape function of the parent cell.
                N1 = B1*(xt-voro.Nodes(n1)).CopyToEigen();
                N1.coeffRef(0) = N1.coeff(0) + 1.;

                // Compute shape function of the neighbor cell.
                N2 = B2*(xt-voro.Nodes(n2)).CopyToEigen();
                N2.coeffRef(0) = N2.coeff(0) + 1.;

                // Apply internal matrix numerical correction.
                lambda1 = material->InterDiffusionTensors(n1);
                lambda2 = material->InterDiffusionTensors(n2);
                Cint11 = L * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (int_penalty/he)*(N1*N1.transpose()));
                Cint12 = L * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (int_penalty/he)*(N1*N2.transpose()));
                Cint21 = L * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (int_penalty/he)*(N2*N1.transpose()));
                Cint22 = L * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (int_penalty/he)*(N2*N2.transpose()));


                // Apply bulk matrix numerical correction.
                lambda1 = material->InterDiffusionTensors(n1) + material->ExterDiffusionTensors(n1);
                lambda2 = material->InterDiffusionTensors(n2) + material->ExterDiffusionTensors(n2);
                Cbulk11 = L * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (bulk_penalty/he)*(N1*N1.transpose()));
                Cbulk12 = L * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (bulk_penalty/he)*(N1*N2.transpose()));
                Cbulk21 = L * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (bulk_penalty/he)*(N2*N1.transpose()));
                Cbulk22 = L * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (bulk_penalty/he)*(N2*N2.transpose()));

                // Assemble Cint11/Cint12 and Cbulk11/Cbulk12 matrices.
                for (const auto &gi : neigh1) {
                    auto li = &gi - &neigh1[0];

                    // Cint11 and Cbulk11 matrix.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint11.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk11.coeff(li, lj)) );
                    }

                    // Cint12 and Cbulk12 matrix.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint12.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk12.coeff(li, lj)) );
                    }
                } // End of Assemble Cint11/Cint12 and Cbulk11/Cbulk12 matrices.

                // Assemble Cint21/Cint22 and Cbulk21/Cbulk22 matrices.
                for (const auto &gi : neigh2) {
                    auto li = &gi - &neigh2[0];

                    // Cint21 and Cbulk21 matrices.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint21.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk21.coeff(li, lj)) );
                    }

                    // Cint22 and Cbulk22 matrices.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint22.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk22.coeff(li, lj)) );
                    }
                } // End of Assemble Cint21/Cint22 and Cbulk21/Cbulk22 matrices.
            }
        } // End of Process only internal facets.

        // Add cardiac stiffness matrix correction.
        corr_int_mat.setFromTriplets(std::begin(corr_int_triplets), std::end(corr_int_triplets));
        this->stiff_int_mat_ += corr_int_mat;

        corr_bulk_mat.setFromTriplets(std::begin(corr_bulk_triplets), std::end(corr_bulk_triplets));
        this->stiff_bulk_mat_ += corr_bulk_mat;

        // Clean correction triples and matrix.
        corr_int_triplets.clear();
        corr_bulk_triplets.clear();
        corr_int_mat.setZero();
        corr_bulk_mat.setZero();

    } // End iterate facets in chuncks.

    // Iterate over resting facets if any.
    if (rest > 0) {
        // Clean correction triplets and matrix.
        corr_int_triplets.clear(); corr_bulk_triplets.clear();
        corr_int_mat.setZero(); corr_bulk_mat.setZero();

        for (int fid = divs*chunk; fid != divs*chunk+rest; ++fid) {
            // Process only internal facets.
            if (!voro.Facets(fid).IsFree()) {

                // Compute the facet length.
                double L = voro.Facets(fid).Length(voro.Points());

                // Indices of cells sharing the facet.
                int n1 = voro.Facets(fid).ParentCellId();
                int n2 = voro.Facets(fid).NeighCellId();

                // Compute the facet outward normal vector.
                IMP::Vec<DIM, double> normal = voro.Facets(fid).Normal(voro.Points());
                if (normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0)  normal = -normal;

                // Compute opposite vector to normal.
                ne1 = normal.CopyToEigen();
                ne2 = -ne1;

                // Boundary dependent parameter with unit of length
                double he = std::sqrt(voro.Nodes(n2).Distance2(voro.Nodes(n1)));

                // Get fpm gradients for the cells.
                B1 = fpm.PhiGrad(n1);
                B2 = fpm.PhiGrad(n2);

                // Get neighbor nodes to the nodes enclosed by the cells.
                neigh1 = fpm.Support().InfluenceNodeIds(n1);
                neigh2 = fpm.Support().InfluenceNodeIds(n2);

                // Compute integration point at the center of the facet.
                IMP::Vec<DIM, double> xt = voro.Facets(fid).Centroid(voro.Points());

                // Compute shape function of the parent cell.
                N1 = B1*(xt-voro.Nodes(n1)).CopyToEigen();
                N1.coeffRef(0) = N1.coeff(0) + 1.;

                // Compute shape function of the neighbor cell.
                N2 = B2*(xt-voro.Nodes(n2)).CopyToEigen();
                N2.coeffRef(0) = N2.coeff(0) + 1.;

                // Apply internal matrix numerical correction.
                lambda1 = material->InterDiffusionTensors(n1);
                lambda2 = material->InterDiffusionTensors(n2);
                Cint11 = L * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (int_penalty/he)*(N1*N1.transpose()));
                Cint12 = L * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (int_penalty/he)*(N1*N2.transpose()));
                Cint21 = L * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (int_penalty/he)*(N2*N1.transpose()));
                Cint22 = L * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (int_penalty/he)*(N2*N2.transpose()));

                // Apply bulk matrix numerical correction.
                lambda1 = material->InterDiffusionTensors(n1) + material->ExterDiffusionTensors(n1);
                lambda2 = material->InterDiffusionTensors(n2) + material->ExterDiffusionTensors(n2);
                Cbulk11 = L * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (bulk_penalty/he)*(N1*N1.transpose()));
                Cbulk12 = L * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (bulk_penalty/he)*(N1*N2.transpose()));
                Cbulk21 = L * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (bulk_penalty/he)*(N2*N1.transpose()));
                Cbulk22 = L * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (bulk_penalty/he)*(N2*N2.transpose()));

                // Assemble C11 & C12 matrices.
                for (const auto &gi : neigh1) {
                    auto li = &gi - &neigh1[0];

                    // C11 matrix.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint11.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk11.coeff(li, lj)) );
                    }

                    // C12 matrix.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint12.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk12.coeff(li, lj)) );
                    }
                } // End of Assemble C11 & C12 matrices.

                // Assemble C21 & C22 matrices.
                for (const auto &gi : neigh2) {
                    auto li = &gi - &neigh2[0];

                    // C21 matrix.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint21.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk21.coeff(li, lj)) );
                    }

                    // C22 matrix.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint22.coeff(li, lj)) );
                        corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk22.coeff(li, lj)) );
                    }
                } // End of Assemble C21 & C22 matrices.
            }
        } // End of Process only internal facets.

        // Add cardiac stiffness matrix correction.
        corr_int_mat.setFromTriplets(std::begin(corr_int_triplets), std::end(corr_int_triplets));
        this->stiff_int_mat_ += corr_int_mat;

        corr_bulk_mat.setFromTriplets(std::begin(corr_bulk_triplets), std::end(corr_bulk_triplets));
        this->stiff_bulk_mat_ += corr_bulk_mat;
    }

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::ApplyFpmCorrection3D(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
        const CLOUDEA::Fpm<DIM> &fpm, double int_penalty, double bulk_penalty)
{
    typedef Eigen::Triplet<double> EigTripD;

    // Estimate facets division number for memory handling.
    auto divs = int{1};           // facets division.
    auto rest = int{0};           // resting facets number for not exact division.
    auto chunk = int{50000};     // division chunk size.
    auto facets_num = static_cast<int>(voro.Facets().size());
    if (facets_num > chunk) {
        divs = facets_num / chunk;
        rest = facets_num - divs*chunk;
    } else { chunk = facets_num; }

    // Estimate correction triplet entries number.
    auto entries_num = chunk*fpm.Support().MaxInfluenceNodesNum()*fpm.Support().MaxInfluenceNodesNum();

    // Initialize indices of cells sharing the facet.
    auto n1 = int{0};
    auto n2 = int{0};

    // Initialize neighbor nodes of the nodes inside the facet-sharing cells.
    auto neigh1 = std::vector<int>{};
    auto neigh2 = std::vector<int>{};

    // Initialize normals on facet and shape function.
    auto ne1 = Eigen::VectorXd{};
    auto ne2 = Eigen::VectorXd{};
    auto N1 = Eigen::VectorXd{};
    auto N2 = Eigen::VectorXd{};

    // Initialize shape function derivatives, diffusion tensors, and correction matrices.
    auto B1 = Eigen::MatrixXd{};
    auto B2 = Eigen::MatrixXd{};
    auto lambda1 = Eigen::MatrixXd{};
    auto lambda2 = Eigen::MatrixXd{};
    auto Cint11 = Eigen::MatrixXd{};
    auto Cint12 = Eigen::MatrixXd{};
    auto Cint21 = Eigen::MatrixXd{};
    auto Cint22 = Eigen::MatrixXd{};
    auto Cbulk11 = Eigen::MatrixXd{};
    auto Cbulk12 = Eigen::MatrixXd{};
    auto Cbulk21 = Eigen::MatrixXd{};
    auto Cbulk22 = Eigen::MatrixXd{};

    // Initialize correction triplets vectors.
    auto corr_int_triplets = std::vector<EigTripD>{};
    auto corr_bulk_triplets = std::vector<EigTripD>{};
    corr_int_triplets.reserve(entries_num);
    corr_bulk_triplets.reserve(entries_num);

    // Initialize correction matrices.
    auto corr_int_mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(this->stiff_int_mat_.rows(), this->stiff_int_mat_.cols());
    auto corr_bulk_mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(this->stiff_bulk_mat_.rows(), this->stiff_bulk_mat_.cols());

    // Iterate facets in chuncks.
    for (int i = 0; i != divs; ++i) {
        for (int fid = i*chunk; fid != (i+1)*chunk; ++fid) {
            // Process only internal facets.
            if (!voro.Facets(fid).IsFree()) {

                // Indices of cells sharing the facet.
                n1 = voro.Facets(fid).ParentCellId();
                n2 = voro.Facets(fid).NeighCellId();

                // Get fpm gradients for the cells.
                B1 = fpm.PhiGrad(n1);
                B2 = fpm.PhiGrad(n2);

                // Get neighbor nodes to the nodes enclosed by the cells.
                neigh1 = fpm.Support().InfluenceNodeIds(n1);
                neigh2 = fpm.Support().InfluenceNodeIds(n2);

                // Triangulate the polygon facet.
                auto facet_pointsnum = static_cast<int>(voro.Facets(fid).Connectivity().size());
                auto triangles = std::vector<IMP::Cell<DIM,3>>(facet_pointsnum-2);

                // Construct triangles.
                const auto v0 = voro.Facets(fid).C(0);
                auto tid = int{0}, v1 = int{0}, v2 = int{0};
                for (int i = 1; i != facet_pointsnum-1; ++i) {
                    v1 = voro.Facets(fid).C(i);
                    v2 = voro.Facets(fid).C(i+1);
                    triangles[tid++].SetConnectivity({v0, v1, v2});
                }

                // Apply correction for each triangle of the facet.
                auto tri_edge1 = IMP::Vec<DIM, double>{};
                auto tri_edge2 = IMP::Vec<DIM, double>{};
                auto tri_normal = IMP::Vec<DIM, double>{};
                for (const auto &tri : triangles) {

                    // Compute triangle normal.
                    tri_edge1 = voro.Points(tri.N(1)) - voro.Points(tri.N(0));
                    tri_edge2 = voro.Points(tri.N(2)) - voro.Points(tri.N(0));
                    tri_normal.Set({tri_edge1[1]*tri_edge2[2] - tri_edge1[2]*tri_edge2[1],
                        tri_edge1[2]*tri_edge2[0] - tri_edge1[0]*tri_edge2[2],
                        tri_edge1[0]*tri_edge2[1] - tri_edge1[1]*tri_edge2[0]}
                    );

                    // Compute the triangle area.
                    auto area = 0.5*std::abs(tri_normal.Norm());

                    // Make normal vector unit.
                    tri_normal /= tri_normal.Norm();

                    // Make sure that normal unit vector is outward.
                    if (tri_normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0)  tri_normal = -tri_normal;

                    // Set triangle normal pointing toward the parent and neighbor cell.
                    ne1 = tri_normal.CopyToEigen();
                    ne2 = -ne1;

                    // Boundary dependent parameter with unit of length
                    auto he = std::sqrt(voro.Nodes(n2).Distance2(voro.Nodes(n1)));

                    // Compute integration point at the center of the triangle.
                    auto xt = (voro.Points(tri.N(0)) + voro.Points(tri.N(1)) + voro.Points(tri.N(2))) / 3.;

                    // Compute shape function of the parent cell.
                    N1 = B1*(xt-voro.Nodes(n1)).CopyToEigen();
                    N1.coeffRef(0) = N1.coeff(0) + 1.;

                    // Compute shape function of the neighbor cell.
                    N2 = B2*(xt-voro.Nodes(n2)).CopyToEigen();
                    N2.coeffRef(0) = N2.coeff(0) + 1.;

                    // Apply internal matrix numerical correction.
                    lambda1 = material->InterDiffusionTensors(n1);
                    lambda2 = material->InterDiffusionTensors(n2);
                    Cint11 = area * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (int_penalty/he)*(N1*N1.transpose()));
                    Cint12 = area * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (int_penalty/he)*(N1*N2.transpose()));
                    Cint21 = area * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (int_penalty/he)*(N2*N1.transpose()));
                    Cint22 = area * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (int_penalty/he)*(N2*N2.transpose()));

                    // Apply bulk matrix numerical correction.
                    lambda1 = material->InterDiffusionTensors(n1) + material->ExterDiffusionTensors(n1);
                    lambda2 = material->InterDiffusionTensors(n2) + material->ExterDiffusionTensors(n2);
                    Cbulk11 = area * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (bulk_penalty/he)*(N1*N1.transpose()));
                    Cbulk12 = area * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (bulk_penalty/he)*(N1*N2.transpose()));
                    Cbulk21 = area * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (bulk_penalty/he)*(N2*N1.transpose()));
                    Cbulk22 = area * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (bulk_penalty/he)*(N2*N2.transpose()));

                    // Assemble C11 & C12 matrices.
                    for (const auto &gi : neigh1) {
                        auto li = &gi - &neigh1[0];

                        // C11 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint11.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk11.coeff(li, lj)) );
                        }

                        // C12 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint12.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk12.coeff(li, lj)) );
                        }
                    } // End of Assemble C11 & C12 matrices.

                    // Assemble C21 & C22 matrices.
                    for (const auto &gi : neigh2) {
                        auto li = &gi - &neigh2[0];

                        // C21 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint21.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk21.coeff(li, lj)) );
                        }

                        // C22 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint22.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk22.coeff(li, lj)) );
                        }
                    } // End of Assemble C21 & C22 matrices.
                } // End of Apply correction for each triangle of the facet.
            }
        } // End of Apply numerical correction at internal facets.

        // Add cardiac stiffness matrix correction.
        corr_int_mat.setFromTriplets(std::begin(corr_int_triplets), std::end(corr_int_triplets));
        this->stiff_int_mat_ += corr_int_mat;

        corr_bulk_mat.setFromTriplets(std::begin(corr_bulk_triplets), std::end(corr_bulk_triplets));
        this->stiff_bulk_mat_ += corr_bulk_mat;

        // Clean correction triples and matrix.
        corr_int_triplets.clear();
        corr_bulk_triplets.clear();
        corr_int_mat.setZero();
        corr_bulk_mat.setZero();

    } // End iterate facets in chuncks.

    // Iterate over resting facets if any.
    if (rest > 0) {
        for (int fid = divs*chunk; fid != divs*chunk+rest; ++fid) {
            // Process only internal facets.
            if (!voro.Facets(fid).IsFree()) {

                // Indices of cells sharing the facet.
                n1 = voro.Facets(fid).ParentCellId();
                n2 = voro.Facets(fid).NeighCellId();

                // Get fpm gradients for the cells.
                B1 = fpm.PhiGrad(n1);
                B2 = fpm.PhiGrad(n2);

                // Get neighbor nodes to the nodes enclosed by the cells.
                neigh1 = fpm.Support().InfluenceNodeIds(n1);
                neigh2 = fpm.Support().InfluenceNodeIds(n2);

                // Triangulate the polygon facet.
                auto facet_pointsnum = static_cast<int>(voro.Facets(fid).Connectivity().size());
                auto triangles = std::vector<IMP::Cell<DIM,3>>(facet_pointsnum-2);

                // Construct triangles.
                const auto v0 = voro.Facets(fid).C(0);
                auto tid = int{0}, v1 = int{0}, v2 = int{0};
                for (int i = 1; i != facet_pointsnum-1; ++i) {
                    v1 = voro.Facets(fid).C(i);
                    v2 = voro.Facets(fid).C(i+1);
                    triangles[tid++].SetConnectivity({v0, v1, v2});
                }

                // Apply correction for each triangle of the facet.
                auto tri_edge1 = IMP::Vec<DIM, double>{};
                auto tri_edge2 = IMP::Vec<DIM, double>{};
                auto tri_normal = IMP::Vec<DIM, double>{};
                for (const auto &tri : triangles) {

                    // Compute triangle normal.
                    tri_edge1 = voro.Points(tri.N(1)) - voro.Points(tri.N(0));
                    tri_edge2 = voro.Points(tri.N(2)) - voro.Points(tri.N(0));
                    tri_normal.Set({tri_edge1[1]*tri_edge2[2] - tri_edge1[2]*tri_edge2[1],
                                    tri_edge1[2]*tri_edge2[0] - tri_edge1[0]*tri_edge2[2],
                                    tri_edge1[0]*tri_edge2[1] - tri_edge1[1]*tri_edge2[0]});

                    // Compute the triangle area.
                    auto area = 0.5*std::abs(tri_normal.Norm());

                    // Make normal vector unit.
                    tri_normal /= tri_normal.Norm();

                    // Make sure that normal unit vector is outward.
                    if (tri_normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0)  tri_normal = -tri_normal;

                    // Set triangle normal pointing toward the parent and neighbor cell.
                    ne1 = tri_normal.CopyToEigen();
                    ne2 = -ne1;

                    // Boundary dependent parameter with unit of length
                    auto he = std::sqrt(voro.Nodes(n2).Distance2(voro.Nodes(n1)));

                    // Compute integration point at the center of the triangle.
                    auto xt = (voro.Points(tri.N(0)) + voro.Points(tri.N(1)) + voro.Points(tri.N(2))) / 3.;

                    // Compute shape function of the parent cell.
                    N1 = B1*(xt-voro.Nodes(n1)).CopyToEigen();
                    N1.coeffRef(0) = N1.coeff(0) + 1.;

                    // Compute shape function of the neighbor cell.
                    N2 = B2*(xt-voro.Nodes(n2)).CopyToEigen();
                    N2.coeffRef(0) = N2.coeff(0) + 1.;

                    // Apply internal matrix numerical correction.
                    lambda1 = material->InterDiffusionTensors(n1);
                    lambda2 = material->InterDiffusionTensors(n2);
                    Cint11 = area * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (int_penalty/he)*(N1*N1.transpose()));
                    Cint12 = area * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (int_penalty/he)*(N1*N2.transpose()));
                    Cint21 = area * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (int_penalty/he)*(N2*N1.transpose()));
                    Cint22 = area * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (int_penalty/he)*(N2*N2.transpose()));

                    // Apply bulk matrix numerical correction.
                    lambda1 = material->InterDiffusionTensors(n1) + material->ExterDiffusionTensors(n1);
                    lambda2 = material->InterDiffusionTensors(n2) + material->ExterDiffusionTensors(n2);
                    Cbulk11 = area * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (bulk_penalty/he)*(N1*N1.transpose()));
                    Cbulk12 = area * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (bulk_penalty/he)*(N1*N2.transpose()));
                    Cbulk21 = area * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (bulk_penalty/he)*(N2*N1.transpose()));
                    Cbulk22 = area * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (bulk_penalty/he)*(N2*N2.transpose()));

                    // Assemble C11 & C12 matrices.
                    for (const auto &gi : neigh1) {
                        auto li = &gi - &neigh1[0];

                        // C11 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint11.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk11.coeff(li, lj)) );
                        }

                        // C12 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint12.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk12.coeff(li, lj)) );
                        }
                    } // End of Assemble C11 & C12 matrices.

                    // Assemble C21 & C22 matrices.
                    for (const auto &gi : neigh2) {
                        auto li = &gi - &neigh2[0];

                        // C21 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint21.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk21.coeff(li, lj)) );
                        }

                        // C22 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_int_triplets.emplace_back( EigTripD(gi, gj, Cint22.coeff(li, lj)) );
                            corr_bulk_triplets.emplace_back( EigTripD(gi, gj, Cbulk22.coeff(li, lj)) );
                        }
                    } // End of Assemble C21 & C22 matrices.
                } // End of Apply correction for each triangle of the facet.
            }
        } // End of Apply numerical correction at internal facets.

        // Add cardiac stiffness matrix correction.
        corr_int_mat.setFromTriplets(std::begin(corr_int_triplets), std::end(corr_int_triplets));
        this->stiff_int_mat_ += corr_int_mat;

        corr_bulk_mat.setFromTriplets(std::begin(corr_bulk_triplets), std::end(corr_bulk_triplets));
        this->stiff_bulk_mat_ += corr_bulk_mat;
    }

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::ComputeCriticalTimeStep()
{

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::Compute(const std::shared_ptr<ElectricBasic<DIM>> &material, const std::vector<Stimulus> &stimuli)
{
    // Check if stiffness matrices are available.
    if (this->stiff_int_mat_.rows() != static_cast<int>(this->Cells().size()) ||
        this->stiff_bulk_mat_.rows() != static_cast<int>(this->Cells().size())) {
        throw std::invalid_argument(Logger::Error("Could not compute bidodomain solver solution. Stiffness matrices have not been set correctly."));
    }
    // Check if minimum time step was set.
    if (this->DtMin() < 0.) {
        throw std::runtime_error(Logger::Error("Could not compute bidodomain solver solution. The minimum time step was not set."));
    }
    // Check if simulation time was set.
    if (this->Dt() < 0.) {
        throw std::runtime_error(Logger::Error("Could not compute bidodomain solver solution. The time step was not set."));
    }
    // Check if simulation time was set.
    if (this->SimulationTime() < 0.) {
        throw std::runtime_error(Logger::Error("Could not compute bidodomain solver solution. The simulation time was not set."));
    }
    // Check if output steps number was set.
    if (this->OutputSteps() < 0) {
        throw std::runtime_error(Logger::Error("Could not compute bidodomain solver solution. The output steps number was not set."));
    }
    if (this->Cells().size() != static_cast<std::size_t>(material->NodesNum()+this->ConductSystem().NodesNum())) {
        throw std::runtime_error(Logger::Error("Could not compute bidodomain solver solution. Cells have not been initialized."));
    }

    for (const auto &stim : stimuli) {
        if (!stim.HasStimulatedNodeSet() && this->ConductSystem().NodesNum() == 0) {
            throw std::runtime_error(Logger::Error("Could not compute bidodomain solver solution. Simulated nodes for applied stimulus not set."));
        }
    }

    // Initialize potential output container.
    std::size_t cells_num = this->Cells().size();
    this->vout_.clear();
    this->vout_.resize(std::ceil(this->SimulationSteps()/this->OutputSteps())+1, Eigen::VectorXd::Zero(cells_num));

    // Construct lumped mass matrix.
    typedef Eigen::Triplet<double> T;
    std::vector<T> mass_triplets;
    mass_triplets.reserve(cells_num);
    for (Eigen::Index i = 0; i != this->mass_vec_.size(); ++i)
        mass_triplets.emplace_back(T(i,i,this->mass_vec_.coeffRef(i)));
    Eigen::SparseMatrix<double, Eigen::RowMajor> mass_mat(cells_num,cells_num);
    mass_mat.setFromTriplets(std::begin(mass_triplets), std::end(mass_triplets));
    mass_triplets.clear(); mass_triplets.shrink_to_fit();

    // Initialize solution vector.
    Eigen::VectorXd sol = Eigen::VectorXd::Zero(2*cells_num);

    // Initialize the potential of the nodal cells.
    for (auto &n_cell : this->Cells()) {
        auto id = &n_cell - &this->Cells()[0];

        // Set initial potential values.
        sol.coeffRef(id) = n_cell->V();
    }

    // Store the initial nodal potential to the output container.
    this->vout_[0] = sol.head(cells_num);

    // Choose theta parameter
    auto theta = double{0.};
    if (this->DiffusionSolver() == DiffusionSolverType::be) { theta = 1.; }
    else if (this->DiffusionSolver() == DiffusionSolverType::cn) { theta = 0.5;}
    else {
        auto error_msg = "Could not compute bidomain model. Unknown solver type. Supported: [be | cn]";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Implicit solution.
    auto steps_counter = int{0};
    auto pos = int{1};
    for (int step = 1; step <= this->SimulationSteps(); ++step) {
        // Compute Reaction diffusion with Godunov Operator Splitting and Backward Euler.
        this->ComputeReaction(stimuli, step, this->Dt(), sol, this->EditCells());
        this->ComputeDiffusionTheta(mass_mat, theta, sol);
        if (step % this->OutputSteps() == 0) {
            this->vout_[pos++] = sol.head(cells_num);
            std::cout << Logger::Message("Completed simulation steps: ") << step << "/" << this->SimulationSteps() << "\r" << std::flush;
        }
        steps_counter++;
    }
    std::cout << "\n";

    // Store final state of current potential if the last step was not saved in the loop.
    if (this->OutputSteps() != 0) {
        if (steps_counter % this->OutputSteps() != 0) {
            this->vout_.emplace_back(sol.head(cells_num));
        }
    } else {
        // Store final state if output_steps_ is zero.
        this->vout_.emplace_back(sol.head(cells_num));
    }

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::FictitiousValuesToReal(const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx)
{

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::ComputeReaction(const std::vector<Stimulus> &stimuli, int step, double dt,
        Eigen::VectorXd &solution, std::vector<std::unique_ptr<EpBasic>> &cells)
{
    // Set thread ranges.
    this->thread_loop_manager_.SetLoopRanges(cells.size(), this->threads_number_);

    // Multithreaded reaction term computation.
    std::vector<std::thread> threads;
    if (this->AdaptiveReaction()) {
        for (std::size_t t = 0; t != this->threads_number_; ++t) {
            threads.emplace_back(std::thread(&Bidomain::ReactionCallbackAdaptive, this, t,
                                 std::cref(stimuli), step, dt, std::ref(solution), std::ref(cells)));
        }
    }
    else {
        for (std::size_t t = 0; t != this->threads_number_; ++t) {
            threads.emplace_back(std::thread(&Bidomain::ReactionCallbackStandard, this, t,
                                 std::cref(stimuli), step, dt, std::ref(solution), std::ref(cells)));
        }
    }
    // Join threads.
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::ReactionCallbackStandard(std::size_t thread_id, const std::vector<Stimulus> &stimuli,
        int step, double dt, Eigen::VectorXd &solution, std::vector<std::unique_ptr<EpBasic>> &cells)
{
    // Iterate over nodes (nodal cells) for reaction term.
    for (auto i = this->thread_loop_manager_.LoopStartId(thread_id);
         i != this->thread_loop_manager_.LoopEndId(thread_id); ++i) {

        // Compute the current applied to the nodal cell from the active stimuli.
        double stimulus_current = 0.;
        for (const auto &stimulus : stimuli) {
            if (stimulus.IsActive(step*this->Dt(), this->SimulationTime()) && stimulus.StimulatedNodes().coeff(i) == 1) {
                stimulus_current += stimulus.Amplitude();
            }
        }

        // Update the time varying parameters of the cell's ep model if any.
        if (this->CellVarParamFlags()[i] != -1)
            cells[i]->UpdateTimeVaryingPrms(this->CellVarParamGroups()[this->CellVarParamFlags()[i]], dt);

        // Compute new cell model state  and update the cell membrane potential.
        cells[i]->Compute(solution.coeff(i), dt, stimulus_current);
        cells[i]->SetV(ALGORITHM::ForwardEuler(solution.coeff(i), dt, cells[i]->dVdt()));

        // Update the solution value that correspond to the cell.
        solution.coeffRef(i) = cells[i]->V();

    } // End of Iterate over nodes (nodal cells) for reaction term.

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::ReactionCallbackAdaptive(std::size_t thread_id, const std::vector<Stimulus> &stimuli,
        int step, double dt, Eigen::VectorXd &solution, std::vector<std::unique_ptr<EpBasic>> &cells)
{

    // Compute maximum adaptive time step multiplication factor.
    int kappa_max = static_cast<int>(std::ceil(dt / this->DtMin()));

    // Iterate over nodes (nodal cells) for reaction term.
    for (auto i = this->thread_loop_manager_.LoopStartId(thread_id);
         i != this->thread_loop_manager_.LoopEndId(thread_id); ++i) {

        // Initialize the current applied to the nodal cell from the active stimuli.
        double stimulus_current = 0.;

        // Compute the current applied to the nodal cell from the active stimuli.
        for (const auto &stimulus : stimuli) {
            if (stimulus.IsActive(step*this->Dt(), this->SimulationTime()) && stimulus.StimulatedNodes().coeff(i) == 1) {
                stimulus_current += stimulus.Amplitude();
            }
        }

        // Initialize adaptive time step multiplication factor to maximum value.
        int kappa = kappa_max;
        double dt_adapt = this->DtMin();
        double t_cell = 0.;

        // Initialize the adaptive potential value.
        double v_adapt = solution.coeff(i);

        int safe_guard = 100000*kappa_max, steps_done = 0;
        while (t_cell < dt) {

            // Ensure that integration time is exactly equal to the total time step.
            if (t_cell + dt_adapt > dt) { dt_adapt = dt - t_cell; }

            // Update the time varying parameters of the cell's ap model if any.
            if (this->CellVarParamFlags()[i] != -1) {
               cells[i]->UpdateTimeVaryingPrms(this->CellVarParamGroups()[this->CellVarParamFlags()[i]], dt_adapt);
            }

            // Compute the new transmembrane potential for the cell.
            cells[i]->Compute(v_adapt, dt_adapt, stimulus_current);

            // Compute the transmembrane potential time derivative.
            double dvdt = cells[i]->dVdt();

            // Update the adaptive potential value.
            v_adapt += dt_adapt*dvdt;

            // Advance time in cell
            t_cell += dt_adapt;

            // Update the kappa factor value.
            if (dvdt > 0.) { kappa = 5 + static_cast<int>(std::floor(std::abs(dvdt))); }
            else { kappa = 1 + static_cast<int>(std::floor(std::abs(dvdt))); }

            // Check kappa factor range.
            if (kappa > kappa_max) { kappa = kappa_max; }

            // Update adaptive dt.
            dt_adapt = dt / kappa;

            steps_done++;
            if (steps_done > safe_guard) {
                std::string error_msg = "Error occurred during adaptive integration of the reaction term due to NaN in the potential derivative computation.";
                throw std::runtime_error(Logger::Error(error_msg));
            }
        }

        // Set the final transmembrane potential value to the cell.
        cells[i]->SetV(v_adapt);

        // Update the nodal potential value.
        solution.coeffRef(i) = v_adapt;

    } // End of Iterate over nodes (nodal cells) for reaction term.

}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::ComputeDiffusionTheta(const Eigen::SparseMatrix<double, Eigen::RowMajor> &mass_mat,
        double theta, Eigen::VectorXd &solution)
{
    Eigen::initParallel();
    Eigen::setNbThreads(std::thread::hardware_concurrency()-1);

    // Sparse matrices type definitions.
    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpRowMat;
    typedef Eigen::SparseMatrix<double> SpColMat;

    // Construct new matrices for Backward Euler solution.
    Eigen::Index n = mass_mat.rows();
    SpRowMat dtAi = this->Dt()*this->stiff_int_mat_;

    auto A12 = SpColMat(n,2*n);
    A12.leftCols(n) = mass_mat + theta*dtAi;
    A12.rightCols(n) = dtAi;

    auto A34 = SpColMat(n,2*n);
    A34.leftCols(n) = dtAi.transpose();
    A34.rightCols(n) = (this->Dt()/theta)*this->stiff_bulk_mat_;

    auto A = SpRowMat(2*n,2*n);
    A.topRows(n) = SpRowMat(A12);
    A.bottomRows(n) = SpRowMat(A34);

    A12.setZero(); A12.data().squeeze();
    A34.setZero(); A34.data().squeeze();

    // Construct right-hand side vector.
    Eigen::VectorXd b(2*n);
    b.head(n) = (mass_mat - (1.0-theta)*dtAi) * solution.head(n);
    b.tail(n) = ((-(1.0-theta)/theta)*dtAi) * solution.tail(n);
    dtAi.setZero(); dtAi.data().squeeze();

    // Solve bidomain algebraic system.
    auto solver = Eigen::ConjugateGradient<SpRowMat>{};
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        auto error_msg = "Could not compute bidomain model. System matrix is not invertible.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
    solution = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
}


template<short DIM, short CELL_NODES>
void Bidomain<DIM, CELL_NODES>::DiffusionCallback(std::size_t thread_id, const Eigen::SparseMatrix<double, Eigen::RowMajor> &stiff_mat,
        const Eigen::VectorXd &mass_vec, double dt, Eigen::VectorXd &v_nodal)
{
}


} // End of namespace ELECTRA


#endif //ELECTRA_ENGINE_PHYSICS_BIDOMAIN_TPP_