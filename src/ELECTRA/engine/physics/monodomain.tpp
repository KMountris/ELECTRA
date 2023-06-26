/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#ifndef ELECTRA_ENGINE_PHYSICS_MONODOMAIN_TPP_
#define ELECTRA_ENGINE_PHYSICS_MONODOMAIN_TPP_


#include "ELECTRA/engine/physics/monodomain.hpp"


namespace ELECTRA {


template<short DIM, short CELL_NODES>
Monodomain<DIM, CELL_NODES>::Monodomain() : vout_(), stiff_mat_(), mass_vec_(), thread_loop_manager_(), threads_number_(0.), thread_mutex_()
{
    // Get the number of parallel threads.
    const std::size_t available_threads = std::thread::hardware_concurrency()-1;
    this->threads_number_ = std::max(available_threads, 1ul);
}


template<short DIM, short CELL_NODES>
Monodomain<DIM, CELL_NODES>::~Monodomain()
{}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::AssembleMatrices(const IMP::Mesh<DIM, CELL_NODES> &mesh, const std::shared_ptr<ElectricBasic<DIM>> &material)
{
    if (static_cast<int>(material->TransDiffusionTensors().size()) != mesh.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for monodomain model. Material transmembrane diffusion tensors are not available for all mesh nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Clear the stiffness matrix and the lumped mass vector.
    int n = material->NodesNum() + this->ConductSystem().NodesNum();
    this->stiff_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n, n);
    this->mass_vec_.setZero(n);

    // Compute the cell dimensions.
    constexpr short CELL_DIM = (DIM > 1 && CELL_NODES == 2) ? 1 : DIM;

    // Create a finite element of the given type.
    std::unique_ptr<CLOUDEA::BaseFem<CELL_DIM>> elem = CLOUDEA::FemFactory<CELL_DIM>::Create(mesh.CellsShape());

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

    // Initialize triplets for assembly.
    typedef Eigen::Triplet<double> T;
    std::vector<T> stiff_assemble_triplet;
    stiff_assemble_triplet.reserve(n * CELL_NODES);

    // Initialize container to store the cell node coordinates.
    std::vector<IMP::Vec<CELL_DIM, double>> nodes(CELL_NODES, IMP::Vec<CELL_DIM, double>());

    // Iterate over the mesh cells.
    Eigen::MatrixXd mean_diffusion_tensor(CELL_DIM, CELL_DIM);
    double mean_capacitance;
    for (const auto &cell : mesh.Cells()) {

        // Reset mean diffusion tensor and mean capacitance values for the cell.
        mean_diffusion_tensor.setZero();
        mean_capacitance = 0.;

        // Compute mean cell diffusion tensor and capacitance.
        for (const auto &nid : cell.Connectivity()) {
            mean_diffusion_tensor += material->TransDiffusionTensors()[nid];
            mean_capacitance += material->Capacitance()[nid];
        }
        mean_diffusion_tensor /= CELL_NODES;
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
        Eigen::MatrixXd stiff_cell = Eigen::MatrixXd::Zero(CELL_NODES, CELL_NODES);
        Eigen::MatrixXd mass_cell = Eigen::MatrixXd::Zero(CELL_NODES, CELL_NODES);
        for (const auto &qweight : elem->Quadrature().Weights()) {
            auto id = std::distance(&elem->Quadrature().Weights()[0], &qweight);
            stiff_cell += qweight*elem->DetJacobians()[id] *
                         (elem->Derivs()[id]*mean_diffusion_tensor*elem->Derivs()[id].transpose());
            mass_cell += mean_capacitance*qweight *
                         elem->DetJacobians()[id]*elem->ShapeFunctions()[id]*elem->ShapeFunctions()[id].transpose();
        }

        // Collect the cell stiffness and mass contribution.
        for (short local_i = 0; local_i != CELL_NODES; ++local_i) {
            auto global_i = cell.N(local_i);

            for (short local_j = 0; local_j != CELL_NODES; ++local_j) {
                auto global_j = cell.N(local_j);

                stiff_assemble_triplet.emplace_back(T(global_i, global_j, stiff_cell(local_i, local_j)));
                this->mass_vec_[global_i] += mass_cell(local_i, local_j);
            }
        } // End Collect stiffness contribution.

    } // End of Iterate over the mesh cells.

    // Add conduction system contribution.
    if (this->ConductSystem().NodesNum() > 0) {
        // Create 1D finite element for the conduction system segments.
        CLOUDEA::FemD1v2<1> seg_elem;
        seg_elem.SetQuadrature({2});

        // Initialize container to store node coordinates for the conduction system segments.
        std::vector<IMP::Vec<1, double>> seg_nodes(2, IMP::Vec<1, double>());

        // Initialize the diffusivity value of the conduction system segment.
        double seg_diff = 0.;

        // Iterate the conduction system segments.
        for (const auto &seg : this->ConductSystem().Segments()) {

            // Set the segment diffusivity value.
            seg_diff = 0.5 * (this->ConductSystem().DiffusionCoeffs(seg.N(0)) + this->ConductSystem().DiffusionCoeffs(seg.N(1)));

            // Get the segment nodes coordinates.
            seg_nodes[0][0] = 0.;
            seg_nodes[1][0] = std::sqrt( this->ConductSystem().Nodes(seg.N(1)).Distance2(this->ConductSystem().Nodes(seg.N(0))) );

            // Compute the jacobian, shape function and derivatives for the current segment.
            seg_elem.ComputeDerivsNatural();  seg_elem.ComputeJacobians(seg_nodes);
            seg_elem.ComputeShapeFunctions(); seg_elem.ComputeDerivs();

            // Construct the local stiffness and mass over the segment's quadrature.
            Eigen::Matrix2d stiff_seg = Eigen::Matrix2d::Zero(2, 2);
            Eigen::Matrix2d mass_seg = Eigen::Matrix2d::Zero(2, 2);
            for (const auto &qweight : seg_elem.Quadrature().Weights()) {
                auto id = std::distance(&seg_elem.Quadrature().Weights()[0], &qweight);
                stiff_seg += qweight * seg_elem.DetJacobians()[id] * (seg_elem.Derivs()[id]*seg_diff*seg_elem.Derivs()[id].transpose());
                mass_seg += qweight * seg_elem.DetJacobians()[id] * (seg_elem.ShapeFunctions()[id]*seg_elem.ShapeFunctions()[id].transpose());
            }

            // Collect the seg stiffness and mass contribution.
            for (short local_i = 0; local_i != 2; ++local_i) {
                // Add the padding to account that the conduction system nodes follow the material nodes.
                auto global_i = material->NodesNum() + seg.N(local_i);

                for (short local_j = 0; local_j != 2; ++local_j) {
                    auto global_j = material->NodesNum() + seg.N(local_j);

                    stiff_assemble_triplet.emplace_back(T(global_i, global_j, stiff_seg(local_i, local_j)));
                    this->mass_vec_.coeffRef(global_i) += mass_seg(local_i, local_j);

                }
            } // End Collect stiffness contribution.
        } // End Iterate the conduction system segments.


        // Add connector elements at pmj sites.
        double pmj_diff_coeff = 0.;
        for (const auto &pmj : this->ConductSystem().Pmj()) {
            int pmj_id = &pmj - &this->ConductSystem().Pmj()[0];

            // Index of the purkinje (left) corner of the connector element.
            int n1 = material->NodesNum() + this->ConductSystem().TerminalNodeIds()[pmj_id];

            // Get the pmj diffisovity coefficient
            pmj_diff_coeff = this->ConductSystem().PmjDiffusionCoeffs(pmj_id);

            // Iterate over the indices of the connected myocardial nodes (right corner of connector element).
            for (const auto &n2 : pmj) {
                // Add contribution in the stiffness matrix assemblying triplet.
                stiff_assemble_triplet.emplace_back(T(n1, n1, pmj_diff_coeff));
                stiff_assemble_triplet.emplace_back(T(n1, n2, -pmj_diff_coeff));
                stiff_assemble_triplet.emplace_back(T(n2, n1, -pmj_diff_coeff));
                stiff_assemble_triplet.emplace_back(T(n2, n2, pmj_diff_coeff));
            }
        }

    } // End Add conduction system contribution.

    stiff_assemble_triplet.shrink_to_fit();

    // Set cardiac stiffness matrix.
    this->stiff_mat_.setFromTriplets(std::begin(stiff_assemble_triplet), std::end(stiff_assemble_triplet));

    // Compute the stable time step for the diffusion term of the monodomain model.
    this->ComputeCriticalTimeStep();
}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::AssembleMatrices(const IMP::Grid<DIM, CELL_NODES> &grid, const std::shared_ptr<ElectricBasic<DIM>> &material,
        const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx, const IMP::NodeSet &neumann_set)
{
    if (static_cast<int>(mfree_approx->Phi().size()) != grid.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for monodomain model. Meshfree approximation is not available for all grid nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    if (static_cast<int>(material->TransDiffusionTensors().size()) != grid.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for monodomain model. Material transmembrane diffusion tensors are not available for all grid nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Clear the conductivity matrix.
    int n = grid.NodesNum() + this->ConductSystem().NodesNum();
    this->mass_vec_.resize(n);

    // Set flag "0" for traction free and "1" for neumann nodes and degrees of freedom.
    int dof_num = DIM*n;
    int neu_nodes_num = 0;
    int neu_dof_num = 0;
    Eigen::SparseVector<int> neumann_node_flags(grid.NodesNum());
    Eigen::SparseVector<int> dof_flags(dof_num);
    for (const auto nid : neumann_set.NodeIds()) {
        // Set flag to 1 for Neumann node.
        neumann_node_flags.coeffRef(nid) = 1;
        neu_nodes_num++;

        // Set flag to 1 for Neumann dof.
        for (short d = 0; d != DIM; ++d) {
            dof_flags.coeffRef(DIM*nid+d) = 1;
            neu_dof_num++;
        }

    }

    // Initialize triplets for sparse mats assembly.
    typedef Eigen::Triplet<double> T;
    std::vector<T> Q_triplet, Ks_triplet, Ka_triplet;
    Q_triplet.reserve(neu_dof_num*mfree_approx->Support().MaxInfluenceNodesNum());
    Ks_triplet.reserve(dof_num*mfree_approx->Support().MaxInfluenceNodesNum());
    Ka_triplet.reserve(dof_num*mfree_approx->Support().MaxInfluenceNodesNum());

    // Iterate over the grid nodes.
    int rep_n = 0;                // Neumann entries counter.
    Eigen::VectorXd temp;      // Temporary vector for Ka contribution.
    for (const auto &node : grid.Nodes()) {

        // Get the index of the current node of the grid.
        auto n_id = &node - &grid.Nodes()[0];

        // Iterate over support nodes for the current node.
        int pos = 0;
        for (const auto &neigh_id : mfree_approx->Support().InfluenceNodeIds(n_id)) {

            // Get the position of the current neighbor node in the container.
            pos = &neigh_id - &mfree_approx->Support().InfluenceNodeIds(n_id)[0];

            // Compute temporal vector to be assigned in the Ka matrix.
            temp = - (material->TransDiffusionTensors(n_id) * mfree_approx->PhiGrad(n_id).row(pos).transpose());

            // Populate Ks, Ka triplets for neumann and free nodes.
            for (short d = 0; d != DIM; ++d) {
                Ks_triplet.emplace_back( T(n_id, DIM*neigh_id+d, mfree_approx->PhiGrad(n_id)(pos, d)) );
                Ka_triplet.emplace_back( T(DIM*n_id+d, neigh_id, temp(d)) );
            }
        } // End of Iterate over support nodes for the current node.

        // Compute inverse normals matrix.
        if (neumann_node_flags.coeff(n_id) == 1) {
            double norm = 0.;
            for (short d = 0; d != DIM; ++d) { norm += grid.NodeNormals().coeff(n_id,d)*grid.NodeNormals().coeff(n_id,d); }

            for (short d = 0; d != DIM; ++d) {
                // Populate normals matrix.
                Q_triplet.emplace_back( T(rep_n, rep_n, 1. / (1. + 1e6*grid.NodeNormals().coeff(n_id,d)*grid.NodeNormals().coeff(n_id,d)/norm)) );
                rep_n++;
            }
        }

        // Populate capacitance lumped vector.
        this->mass_vec_.coeffRef(n_id) = material->Capacitance(n_id)*mfree_approx->Phi(n_id).sum();

    } // End of Iterate over the grid nodes.
    neumann_node_flags.setZero(); neumann_node_flags.data().squeeze();

    // Set Q matrix from triplets.
    Eigen::SparseMatrix<double, Eigen::RowMajor> Q(neu_dof_num, neu_dof_num);
    Q.setFromTriplets(std::begin(Q_triplet), std::end(Q_triplet));
    Q_triplet.clear(); Q_triplet.shrink_to_fit();

    // Set Ks matrix from triplets.
    Eigen::SparseMatrix<double> Ks(n, dof_num);
    Ks.setFromTriplets(std::begin(Ks_triplet), std::end(Ks_triplet));
    Ks_triplet.clear(); Ks_triplet.shrink_to_fit();

    // Set Ka matrix from triplets.
    Eigen::SparseMatrix<double, Eigen::RowMajor> Ka(dof_num, n);
    Ka.setFromTriplets(std::begin(Ka_triplet), std::end(Ka_triplet));
    Ka_triplet.clear(); Ka_triplet.shrink_to_fit();

    // Initialize triplets to separate K for neumann and traction-free dof.
    std::vector<T> Kn_triplet, Kf_triplet;
    Kn_triplet.reserve(dof_num*mfree_approx->Support().MaxInfluenceNodesNum());
    Kf_triplet.reserve(dof_num*mfree_approx->Support().MaxInfluenceNodesNum());

    // Populate Kn_triplet and Kf_triplet from Ks.
    int kn_pad = 0; int kf_pad = 0;
    for (int i = 0; i != Ks.outerSize(); ++i) {
        if (dof_flags.coeff(i) == 1) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Ks,i); it; ++it) {
                Kn_triplet.emplace_back( T(it.row(), kn_pad, it.value()) );
            }
            kn_pad++;
        } else {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Ks,i); it; ++it) {
                Kf_triplet.emplace_back( T(it.row(), kf_pad, it.value()) );
            }
            kf_pad++;
        }
    }

        std::cout << "********* pass\n";

    // Set Ks_n matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> Ks_n(n, neu_dof_num);
    Ks_n.setFromTriplets(std::begin(Kn_triplet), std::end(Kn_triplet));

    // Set Ks_f matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> Ks_f(n, dof_num-neu_dof_num);
    Ks_f.setFromTriplets(std::begin(Kf_triplet), std::end(Kf_triplet));

    // Delete original Ks matrix to save memory.
    Ks.setZero(); Ks.data().squeeze();

    // Populate Kn_triplet and Kf_triplet from Ka.
    Kn_triplet.clear(); Kf_triplet.clear();
    kn_pad = 0; kf_pad = 0;
    for (int i = 0; i != Ka.outerSize(); ++i) {
        if (dof_flags.coeff(i) == 1) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Ka,i); it; ++it) {
                Kn_triplet.emplace_back(T(kn_pad, it.col(), it.value()) );
            }
            kn_pad++;
        } else {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Ka,i); it; ++it) {
                Kf_triplet.emplace_back(T(kf_pad, it.col(), it.value()) );
            }
            kf_pad++;
        }
    }
    dof_flags.setZero(); dof_flags.data().squeeze();

    // Delete original Ka matrix to save memory.
    Ka.setZero(); Ka.data().squeeze();

    // Set Ka_n matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> Ka_n(neu_dof_num, n);
    Ka_n.setFromTriplets(std::begin(Kn_triplet), std::end(Kn_triplet));
    Kn_triplet.clear(); Kn_triplet.shrink_to_fit();

    // Set Ka_f matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> Ka_f(dof_num-neu_dof_num, n);
    Ka_f.setFromTriplets(std::begin(Kf_triplet), std::end(Kf_triplet));
    Kf_triplet.clear(); Kf_triplet.shrink_to_fit();

    // Set the A part of the final stiffness matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> A = Ks_n*Q*Ka_n;
    Ks_n.setZero(); Ks_n.data().squeeze();
    Q.setZero(); Q.data().squeeze();
    Ka_n.setZero(); Ka_n.data().squeeze();

    // Set the B part of the final stiffness matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> B = Ks_f*Ka_f;
    Ks_f.setZero(); Ks_f.data().squeeze();
    Ka_f.setZero(); Ka_f.data().squeeze();

    // Set the final stiffness matrix.
    this->stiff_mat_ = A+B;
    A.setZero(); A.data().squeeze();
    B.setZero(); B.data().squeeze();

    // Add conduction system contribution.
    if (this->ConductSystem().NodesNum() > 0) {

        // Initialize triplet to store conduction system stiffness contribution.
        std::vector<T> C_triplet;

        // Create 1D finite element for the conduction system segments.
        CLOUDEA::FemD1v2<1> seg_elem;
        seg_elem.SetQuadrature({2});

        // Initialize container to store node coordinates for the conduction system segments.
        std::vector<IMP::Vec<1, double>> seg_nodes(2, IMP::Vec<1, double>());

        // Initialize the diffusivity value of the conduction system segment.
        double seg_diff = 0.;

        // Iterate the conduction system segments.
        for (const auto &seg : this->ConductSystem().Segments()) {

            // Set the segment diffusivity value.
            seg_diff = 0.5 * (this->ConductSystem().DiffusionCoeffs(seg.N(0)) + this->ConductSystem().DiffusionCoeffs(seg.N(1)));

            // Get the segment nodes coordinates.
            seg_nodes[0][0] = 0.;
            seg_nodes[1][0] = std::sqrt( this->ConductSystem().Nodes(seg.N(1)).Distance2(this->ConductSystem().Nodes(seg.N(0))) );

            // Compute the jacobian, shape function and derivatives for the current segment.
            seg_elem.ComputeDerivsNatural();  seg_elem.ComputeJacobians(seg_nodes);
            seg_elem.ComputeShapeFunctions(); seg_elem.ComputeDerivs();

            // Construct the local stiffness and mass over the segment's quadrature.
            Eigen::Matrix2d stiff_seg = Eigen::Matrix2d::Zero(2, 2);
            Eigen::Matrix2d mass_seg = Eigen::Matrix2d::Zero(2, 2);
            for (const auto &qweight : seg_elem.Quadrature().Weights()) {
                auto id = std::distance(&seg_elem.Quadrature().Weights()[0], &qweight);
                stiff_seg += qweight * seg_elem.DetJacobians()[id] * (seg_elem.Derivs()[id]*seg_diff*seg_elem.Derivs()[id].transpose());
                mass_seg += qweight * seg_elem.DetJacobians()[id] * (seg_elem.ShapeFunctions()[id]*seg_elem.ShapeFunctions()[id].transpose());
            }

            // Collect the seg stiffness and mass contribution.
            for (short local_i = 0; local_i != 2; ++local_i) {
                // Add the padding to account that the conduction system nodes follow the grid nodes.
                auto global_i = grid.NodesNum() + seg.N(local_i);

                for (short local_j = 0; local_j != 2; ++local_j) {
                    auto global_j = grid.NodesNum() + seg.N(local_j);

                    C_triplet.emplace_back(T(global_i, global_j, stiff_seg(local_i, local_j)));
                    this->mass_vec_.coeffRef(global_i) += mass_seg(local_i, local_j);

                }
            } // End Collect stiffness contribution.
        } // End Iterate the conduction system segments.


        // Add connector elements at pmj sites.
        double pmj_diff_coeff = 0.;
        for (const auto &pmj : this->ConductSystem().Pmj()) {
            int pmj_id = &pmj - &this->ConductSystem().Pmj()[0];

            // Index of the purkinje (left) corner of the connector element.
            int n1 = grid.NodesNum() + this->ConductSystem().TerminalNodeIds()[pmj_id];

            // Get the pmj diffisuvity coefficient
            pmj_diff_coeff = this->ConductSystem().PmjDiffusionCoeffs(pmj_id);

            // Iterate over the indices of the connected myocardial nodes (right corner of connector element).
            for (const auto &n2 : pmj) {
                // Add contribution in the stiffness matrix assemblying triplet.
                C_triplet.emplace_back(T(n1, n1, pmj_diff_coeff));
                C_triplet.emplace_back(T(n1, n2, -pmj_diff_coeff));
                C_triplet.emplace_back(T(n2, n1, -pmj_diff_coeff));
                C_triplet.emplace_back(T(n2, n2, pmj_diff_coeff));
            }
        }


        // Create conduction system stiffness matrix contribution.
        Eigen::SparseMatrix<double, Eigen::RowMajor> C(n, n);
        C.setFromTriplets(std::begin(C_triplet), std::end(C_triplet));
        C_triplet.clear(); C_triplet.shrink_to_fit();

        // Add to the stiffness matrix of the model.
        this->stiff_mat_ += C;

    } // End Add conduction system contribution.

    // std::cout << std::setprecision(16);
    // for (int k=0; k != this->stiff_mat_.outerSize(); ++k) {
    //     for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(this->stiff_mat_,k); it; ++it) {
    //         std::cout << "(" << it.row()+1 << "," << it.col()+1 << ") " << it.value() << std::endl;
    //     }
    // }

    // Compute the stable time step for the diffusion term of the monodomain model.
    this->ComputeCriticalTimeStep();
}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::AssembleMatrices(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
        const CLOUDEA::Fpm<DIM> &fpm)
{
    if (static_cast<int>(fpm.Phi().size()) != voro.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for monodomain model. Fpm approximation is not available for all voronoi tesselation nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    if (static_cast<int>(material->TransDiffusionTensors().size()) != voro.NodesNum()) {
        std::string error_msg = "Could not assemble matrices for monodomain model. Material transmembrane diffusion tensors are not available for all voronoi tesselation nodes.";
        throw std::invalid_argument(Logger::Error(error_msg));
    }

    // Clear the stiffness matrix and the lumped mass vector.
    int n = material->NodesNum() + this->ConductSystem().NodesNum();

    // Compute mean diffusion value.
    Eigen::VectorXd cell_volumes(n);
    double mean_kappa = 0., total_volume = 0.;
    for (const auto &cell : voro.Cells()) {
        auto id = &cell - &voro.Cells()[0];
        cell_volumes.coeffRef(id) = cell.Measure(voro.Points(), voro.Facets());
        mean_kappa += material->TransDiffusionTensors(id).diagonal().mean()*cell_volumes.coeff(id);
        total_volume += cell_volumes.coeff(id);
    }
    mean_kappa /= total_volume;
    double penal = mean_kappa * fpm.Penalty();

    // Initialize triplets for assembly.
    std::vector<Eigen::Triplet<double>> stiff_triplet, mass_triplet;
    stiff_triplet.reserve(n);  mass_triplet.reserve(n);

    // Assemble stiffness and mass matrices.
    Eigen::MatrixXd stiff_cell, mass_cell;
    for (const auto &cell : voro.Cells()) {
        auto id = &cell - &voro.Cells()[0];

        // Compute stiffness and mass matrices of the cell.
        stiff_cell = fpm.PhiGrad(id)*material->TransDiffusionTensors(id)*fpm.PhiGrad(id).transpose()*cell_volumes.coeff(id);
        mass_cell = fpm.Phi(id)*fpm.Phi(id).transpose()*cell_volumes.coeff(id)*material->Capacitance(id);

        // Assemble matrices.
        for (const auto &gi : fpm.Support().InfluenceNodeIds(id)) {
            auto li = &gi - &fpm.Support().InfluenceNodeIds(id)[0];

            for (const auto &gj : fpm.Support().InfluenceNodeIds(id)) {
                auto lj = &gj - &fpm.Support().InfluenceNodeIds(id)[0];
                stiff_triplet.emplace_back( Eigen::Triplet<double>(gi, gj, stiff_cell(li, lj)) );
                mass_triplet.emplace_back( Eigen::Triplet<double>(gi, gj, mass_cell(li, lj)) );
            }
        }
    } // End of Assemble stiffness and mass matrices.

    // Assemble mass matrix.
    Eigen::SparseMatrix<double, Eigen::RowMajor> mass_mat(n,n);
    mass_mat.setFromTriplets(std::begin(mass_triplet), std::end(mass_triplet));
    mass_triplet.clear(); mass_triplet.shrink_to_fit();

    // Convert to lumped mass vector.
    this->mass_vec_.setZero(n);
    for (int i = 0; i != mass_mat.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator mass(mass_mat, i); mass; ++mass) {
            this->mass_vec_.coeffRef(i) = this->mass_vec_.coeff(i) + mass.value();
        }
    }
    mass_mat.setZero(); mass_mat.data().squeeze();

    // Set cardiac stiffness matrix.
    this->stiff_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n, n);
    this->stiff_mat_.setFromTriplets(std::begin(stiff_triplet), std::end(stiff_triplet));
    stiff_triplet.clear(); stiff_triplet.shrink_to_fit();

    // Apply numerical correction at internal facets.
    switch (DIM) {
        case 2 : this->ApplyFpmCorrection2D(voro, material, fpm, penal); break;

        case 3 : this->ApplyFpmCorrection3D(voro, material, fpm, penal); break;

        default:
            std::string error_msg = "Could not apply numerical correction for FPM assembly. Domain dimensions are not supported.";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

    // Compute the stable time step for the diffusion term of the monodomain model.
    this->ComputeCriticalTimeStep();

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ApplyFpmCorrection2D(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
                                                       const CLOUDEA::Fpm<DIM> &fpm, double penalty)
{
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
    Eigen::MatrixXd B1, B2, lambda1, lambda2, C11, C12, C21, C22;

    // Initialize correction triplets vector and correction matrix.
    std::vector<Eigen::Triplet<double>> corr_triplets;
    corr_triplets.reserve(entries_num);
    auto corr_mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(this->stiff_mat_.rows(), this->stiff_mat_.cols());

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

                // Apply numerical correction.
                lambda1 = material->TransDiffusionTensors(n1);
                lambda2 = material->TransDiffusionTensors(n2);
                C11 = L * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (penalty/he)*(N1*N1.transpose()));
                C12 = L * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (penalty/he)*(N1*N2.transpose()));
                C21 = L * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (penalty/he)*(N2*N1.transpose()));
                C22 = L * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (penalty/he)*(N2*N2.transpose()));

                // Assemble C11 & C12 matrices.
                for (const auto &gi : neigh1) {
                    auto li = &gi - &neigh1[0];

                    // C11 matrix.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C11.coeff(li, lj)) );
                    }

                    // C12 matrix.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C12.coeff(li, lj)) );
                    }
                } // End of Assemble C11 & C12 matrices.

                // Assemble C21 & C22 matrices.
                for (const auto &gi : neigh2) {
                    auto li = &gi - &neigh2[0];

                    // C21 matrix.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C21.coeff(li, lj)) );
                    }

                    // C22 matrix.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C22.coeff(li, lj)) );
                    }
                } // End of Assemble C21 & C22 matrices.
            }
        } // End of Process only internal facets.

        // Add cardiac stiffness matrix correction.
        corr_mat.setFromTriplets(std::begin(corr_triplets), std::end(corr_triplets));
        this->stiff_mat_ += corr_mat;

        // Clean correction triples and matrix.
        corr_triplets.clear();
        corr_mat.setZero();

    } // End iterate facets in chuncks.

    // Iterate over resting facets if any.
    if (rest > 0) {
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

                // Apply numerical correction.
                lambda1 = material->TransDiffusionTensors(n1);
                lambda2 = material->TransDiffusionTensors(n2);
                C11 = L * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (penalty/he)*(N1*N1.transpose()));
                C12 = L * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (penalty/he)*(N1*N2.transpose()));
                C21 = L * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (penalty/he)*(N2*N1.transpose()));
                C22 = L * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (penalty/he)*(N2*N2.transpose()));

                // Assemble C11 & C12 matrices.
                for (const auto &gi : neigh1) {
                    auto li = &gi - &neigh1[0];

                    // C11 matrix.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C11.coeff(li, lj)) );
                    }

                    // C12 matrix.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C12.coeff(li, lj)) );
                    }
                } // End of Assemble C11 & C12 matrices.

                // Assemble C21 & C22 matrices.
                for (const auto &gi : neigh2) {
                    auto li = &gi - &neigh2[0];

                    // C21 matrix.
                    for (const auto &gj : neigh1) {
                        auto lj = &gj - &neigh1[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C21.coeff(li, lj)) );
                    }

                    // C22 matrix.
                    for (const auto &gj : neigh2) {
                        auto lj = &gj - &neigh2[0];
                        corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C22.coeff(li, lj)) );
                    }
                } // End of Assemble C21 & C22 matrices.
            }
        } // End of Process only internal facets.

        // Add cardiac stiffness matrix correction.
        corr_mat.setFromTriplets(std::begin(corr_triplets), std::end(corr_triplets));
        this->stiff_mat_ += corr_mat;
    }
}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ApplyFpmCorrection3D(const IMP::Voronoi<DIM> &voro, const std::shared_ptr<ElectricBasic<DIM>> &material,
                                                       const CLOUDEA::Fpm<DIM> &fpm, double penalty)
{
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
    Eigen::MatrixXd B1, B2, lambda1, lambda2, C11, C12, C21, C22;

    // Initialize correction triplets vector and correction matrix.
    std::vector<Eigen::Triplet<double>> corr_triplets;
    corr_triplets.reserve(entries_num);
    auto corr_mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(this->stiff_mat_.rows(), this->stiff_mat_.cols());

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
                int facet_pointsnum = static_cast<int>(voro.Facets(fid).Connectivity().size());
                std::vector<IMP::Cell<DIM,3>> triangles(facet_pointsnum-2);

                // Construct triangles.
                const int v0 = voro.Facets(fid).C(0);
                int tid = 0, v1 = 0, v2 = 0;
                for (int i = 1; i != facet_pointsnum-1; ++i) {
                    v1 = voro.Facets(fid).C(i);
                    v2 = voro.Facets(fid).C(i+1);
                    triangles[tid++].SetConnectivity({v0, v1, v2});
                }

                // Apply correction for each triangle of the facet.
                IMP::Vec<DIM, double> tri_edge1, tri_edge2, tri_normal;
                for (const auto &tri : triangles) {

                    // Compute triangle normal.
                    tri_edge1 = voro.Points(tri.N(1)) - voro.Points(tri.N(0));
                    tri_edge2 = voro.Points(tri.N(2)) - voro.Points(tri.N(0));
                    tri_normal.Set({tri_edge1[1]*tri_edge2[2] - tri_edge1[2]*tri_edge2[1],
                                    tri_edge1[2]*tri_edge2[0] - tri_edge1[0]*tri_edge2[2],
                                    tri_edge1[0]*tri_edge2[1] - tri_edge1[1]*tri_edge2[0]});

                    // Compute the triangle area.
                    double area = 0.5*std::abs(tri_normal.Norm());

                    // Make normal vector unit.
                    tri_normal /= tri_normal.Norm();

                    // Make sure that normal unit vector is outward.
                    if (tri_normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0)  tri_normal = -tri_normal;

                    // Set triangle normal pointing toward the parent and neighbor cell.
                    ne1 = tri_normal.CopyToEigen();
                    ne2 = -ne1;

                    // Boundary dependent parameter with unit of length
                    double he = std::sqrt(voro.Nodes(n2).Distance2(voro.Nodes(n1)));

                    // Compute integration point at the center of the triangle.
                    IMP::Vec<DIM, double> xt = (voro.Points(tri.N(0)) + voro.Points(tri.N(1)) + voro.Points(tri.N(2))) / 3;

                    // Compute shape function of the parent cell.
                    N1 = B1*(xt-voro.Nodes(n1)).CopyToEigen();
                    N1.coeffRef(0) = N1.coeff(0) + 1.;

                    // Compute shape function of the neighbor cell.
                    N2 = B2*(xt-voro.Nodes(n2)).CopyToEigen();
                    N2.coeffRef(0) = N2.coeff(0) + 1.;

                    // Apply numerical correction.
                    lambda1 = material->TransDiffusionTensors(n1);
                    lambda2 = material->TransDiffusionTensors(n2);
                    C11 = area * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (penalty/he)*(N1*N1.transpose()));
                    C12 = area * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (penalty/he)*(N1*N2.transpose()));
                    C21 = area * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (penalty/he)*(N2*N1.transpose()));
                    C22 = area * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (penalty/he)*(N2*N2.transpose()));

                    // Assemble C11 & C12 matrices.
                    for (const auto &gi : neigh1) {
                        auto li = &gi - &neigh1[0];

                        // C11 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C11.coeff(li, lj)) );
                        }

                        // C12 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C12.coeff(li, lj)) );
                        }
                    } // End of Assemble C11 & C12 matrices.

                    // Assemble C21 & C22 matrices.
                    for (const auto &gi : neigh2) {
                        auto li = &gi - &neigh2[0];

                        // C21 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C21.coeff(li, lj)) );
                        }

                        // C22 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C22.coeff(li, lj)) );
                        }
                    } // End of Assemble C21 & C22 matrices.
                } // End of Apply correction for each triangle of the facet.
            }
        } // End of Apply numerical correction at internal facets.

        // Add cardiac stiffness matrix correction.
        corr_mat.setFromTriplets(std::begin(corr_triplets), std::end(corr_triplets));
        this->stiff_mat_ += corr_mat;

        // Clean correction triples and matrix.
        corr_triplets.clear();
        corr_mat.setZero();

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
                int facet_pointsnum = static_cast<int>(voro.Facets(fid).Connectivity().size());
                std::vector<IMP::Cell<DIM,3>> triangles(facet_pointsnum-2);

                // Construct triangles.
                const int v0 = voro.Facets(fid).C(0);
                int tid = 0, v1 = 0, v2 = 0;
                for (int i = 1; i != facet_pointsnum-1; ++i) {
                    v1 = voro.Facets(fid).C(i);
                    v2 = voro.Facets(fid).C(i+1);
                    triangles[tid++].SetConnectivity({v0, v1, v2});
                }

                // Apply correction for each triangle of the facet.
                IMP::Vec<DIM, double> tri_edge1, tri_edge2, tri_normal;
                for (const auto &tri : triangles) {

                    // Compute triangle normal.
                    tri_edge1 = voro.Points(tri.N(1)) - voro.Points(tri.N(0));
                    tri_edge2 = voro.Points(tri.N(2)) - voro.Points(tri.N(0));
                    tri_normal.Set({tri_edge1[1]*tri_edge2[2] - tri_edge1[2]*tri_edge2[1],
                                    tri_edge1[2]*tri_edge2[0] - tri_edge1[0]*tri_edge2[2],
                                    tri_edge1[0]*tri_edge2[1] - tri_edge1[1]*tri_edge2[0]});

                    // Compute the triangle area.
                    double area = 0.5*std::abs(tri_normal.Norm());

                    // Make normal vector unit.
                    tri_normal /= tri_normal.Norm();

                    // Make sure that normal unit vector is outward.
                    if (tri_normal.CwiseMul(voro.Nodes(n2)-voro.Nodes(n1)).Sum() < 0)  tri_normal = -tri_normal;

                    // Set triangle normal pointing toward the parent and neighbor cell.
                    ne1 = tri_normal.CopyToEigen();
                    ne2 = -ne1;

                    // Boundary dependent parameter with unit of length
                    double he = std::sqrt(voro.Nodes(n2).Distance2(voro.Nodes(n1)));

                    // Compute integration point at the center of the triangle.
                    IMP::Vec<DIM, double> xt = (voro.Points(tri.N(0)) + voro.Points(tri.N(1)) + voro.Points(tri.N(2))) / 3;

                    // Compute shape function of the parent cell.
                    N1 = B1*(xt-voro.Nodes(n1)).CopyToEigen();
                    N1.coeffRef(0) = N1.coeff(0) + 1.;

                    // Compute shape function of the neighbor cell.
                    N2 = B2*(xt-voro.Nodes(n2)).CopyToEigen();
                    N2.coeffRef(0) = N2.coeff(0) + 1.;

                    // Apply numerical correction.
                    lambda1 = material->TransDiffusionTensors(n1);
                    lambda2 = material->TransDiffusionTensors(n2);
                    C11 = area * (-0.5*(N1*ne1.transpose()*lambda1*B1.transpose() + B1*lambda1*ne1*N1.transpose()) + (penalty/he)*(N1*N1.transpose()));
                    C12 = area * (-0.5*(N1*ne1.transpose()*lambda2*B2.transpose() + B1*lambda1*ne2*N2.transpose()) - (penalty/he)*(N1*N2.transpose()));
                    C21 = area * (-0.5*(N2*ne2.transpose()*lambda1*B1.transpose() + B2*lambda2*ne1*N1.transpose()) - (penalty/he)*(N2*N1.transpose()));
                    C22 = area * (-0.5*(N2*ne2.transpose()*lambda2*B2.transpose() + B2*lambda2*ne2*N2.transpose()) + (penalty/he)*(N2*N2.transpose()));

                    // Assemble C11 & C12 matrices.
                    for (const auto &gi : neigh1) {
                        auto li = &gi - &neigh1[0];

                        // C11 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C11.coeff(li, lj)) );
                        }

                        // C12 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C12.coeff(li, lj)) );
                        }
                    } // End of Assemble C11 & C12 matrices.

                    // Assemble C21 & C22 matrices.
                    for (const auto &gi : neigh2) {
                        auto li = &gi - &neigh2[0];

                        // C21 matrix.
                        for (const auto &gj : neigh1) {
                            auto lj = &gj - &neigh1[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C21.coeff(li, lj)) );
                        }

                        // C22 matrix.
                        for (const auto &gj : neigh2) {
                            auto lj = &gj - &neigh2[0];
                            corr_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, C22.coeff(li, lj)) );
                        }
                    } // End of Assemble C21 & C22 matrices.
                } // End of Apply correction for each triangle of the facet.
            }
        } // End of Apply numerical correction at internal facets.

        // Add cardiac stiffness matrix correction.
        corr_mat.setFromTriplets(std::begin(corr_triplets), std::end(corr_triplets));
        this->stiff_mat_ += corr_mat;
    }

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ComputeCriticalTimeStep()
{

    // Check stiffness matrix and lamped mass vector size consistency.
    if (this->stiff_mat_.rows() != this->mass_vec_.rows() || this->mass_vec_.rows() == 0) {
        throw std::invalid_argument(Logger::Error("Could not compute Monodomain solver's critical time step in Monodomain solver."
                                                  " Stiffness matrix and mass lumped vector are not initialized correctly."));
    }

    // Initialize to the maximum machine value.
    double dt_critical = std::numeric_limits<double>::max();

    // Iterate over the rows of the stiffness matrix.
    double row_sum;
    for (int i = 0; i != this->stiff_mat_.outerSize(); ++i) {
        row_sum = 0.;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(this->stiff_mat_, i); it; ++it) {
            row_sum += std::abs(it.value());
        }
        row_sum -= std::abs(this->stiff_mat_.coeff(i,i));

        // Compute the critical time step for the represented node by the current column.
        double crit_dt = this->mass_vec_.coeff(i) / (row_sum + this->stiff_mat_.coeff(i,i));

        // Update stable time step.
        if (crit_dt < dt_critical) { dt_critical = crit_dt; }
    }

    if (dt_critical < 0.) {
        std::string error_msg = "Unexpected error in the critical time step calculation. Time step is negative (critical dt = " +
                                std::to_string(dt_critical) + "). Check stiffness matrix and mass lumped vector.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Reduce the stable time step by a 10% safety factor and store it.
    this->SetDtCritical(0.9*dt_critical);
}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::Compute(const std::shared_ptr<ElectricBasic<DIM>> &material, const std::vector<Stimulus> &stimuli)
{
    // Check if conductivity matrix is available.
    if (this->stiff_mat_.rows() != static_cast<int>(this->Cells().size())) {
        throw std::invalid_argument(Logger::Error("Could not compute solution in Monodomain solver. Stiffness matrix has not been set correctly."));
    }

    // Check if simulation time was set.
    if (this->SimulationTime() < 0.) {
        throw std::runtime_error(Logger::Error("Could not compute solution in Monodomain solver. The simulation time was not set."));
    }

    // Check if output steps number was set.
    if (this->OutputSteps() < 0) {
        throw std::runtime_error(Logger::Error("Could not compute solution in Monodomain solver. The output steps number was not set."));
    }

    if (this->Cells().size() != static_cast<std::size_t>(material->NodesNum()+this->ConductSystem().NodesNum())) {
        throw std::runtime_error(Logger::Error("Could not compute solution in Monodomain solver. Cells have not been initialized."));
    }

    for (const auto &stim : stimuli) {
        if (!stim.HasStimulatedNodeSet() && this->ConductSystem().NodesNum() == 0) {
            throw std::runtime_error(Logger::Error("Could not compute solution in Monodomain solver. Simulated nodes for applied stimulus not set."));
        }
    }

    // Initialize potential output container.
    this->vout_.clear();
    this->vout_.resize(std::ceil(this->SimulationSteps()/this->OutputSteps())+1, Eigen::VectorXd::Zero(this->Cells().size()));

    // Initialize nodal potential vector and cardiac cell potential vector.
    Eigen::VectorXd v_nodal = Eigen::VectorXd::Zero(this->Cells().size());

    // Initialize the potential of the nodal cells.
    for (auto &n_cell : this->Cells()) {
        auto id = &n_cell - &this->Cells()[0];

        // Set initial potential values.
        v_nodal.coeffRef(id) = n_cell->V();
    }

    // Store the initial nodal potential to the output container.
    this->vout_[0] = v_nodal;

    // Iterate over time.
    int steps_counter = 0;
    int pos = 1;
    for (int step = 1; step <= this->SimulationSteps(); ++step) {
        // Increase the steps counter.
        steps_counter++;

        // Compute Reaction diffusion with Strang Operator Splitting.
        this->ComputeDiffusionExplicit(this->stiff_mat_, this->mass_vec_, 0.5*this->Dt(), v_nodal);
        this->ComputeReaction(stimuli, step, v_nodal, this->EditCells());
        this->ComputeDiffusionExplicit(this->stiff_mat_, this->mass_vec_, 0.5*this->Dt(), v_nodal);

        if (step % this->OutputSteps() == 0) {
            // Store updated nodal potential to output container.
            this->vout_[pos++] = v_nodal;
            std::cout << Logger::Message("Completed simulation steps: ") << step << "/" << this->SimulationSteps() << "\r" << std::flush;
        }

    } // End of Iterate over time.
    std::cout << "\n";

    // Store final state of current potential if the last step was not saved in the loop.
    if (this->OutputSteps() != 0) {
        if (steps_counter % this->OutputSteps() != 0) {
            this->vout_.emplace_back(v_nodal);
        }
    }
    else {  // Store final state if output_steps_ is zero.
        this->vout_.emplace_back(v_nodal);
    }

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::FictitiousValuesToReal(const std::unique_ptr<CLOUDEA::Mfree<DIM>> &mfree_approx)
{
    if ( mfree_approx->Phi().size() != static_cast<std::size_t>(this->vout_[0].rows() - this->ConductSystem().NodesNum()) ) {
        std::string error_str = "Could not convert fictitious nodal values to real in monodomain model. Meshfree approximation is not available for all nodes.";
        throw std::runtime_error(Logger::Error(error_str));
    }

    std::size_t nodes_num = mfree_approx->Phi().size();

    // Iterate over potential values for each time step.
    Eigen::VectorXd fict_vals;
    double real_val = 0.;
    for (auto &v : this->vout_) {
        // Temporary copy of fictitious potential values.
        fict_vals = v;

        // Iterate over potential nodal values.
        for (std::size_t n_id = 0; n_id != nodes_num; ++n_id) {

            // Get the fictitious values for the neighbor nodes.
            for (const auto &neigh_id : mfree_approx->Support().InfluenceNodeIds(n_id)) {
                auto i = &neigh_id - &mfree_approx->Support().InfluenceNodeIds(n_id)[0];

                real_val += mfree_approx->Phi(n_id)(i) * fict_vals.coeff(neigh_id);
            }

            v.coeffRef(n_id) = real_val;
            real_val = 0.;
        } // End of Iterate over potential nodal values.

    } // End of Iterate over potential values for each time step.

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ComputeReaction(const std::vector<Stimulus> &stimuli, int step, Eigen::VectorXd &v_nodal,
                                                  std::vector<std::unique_ptr<EpBasic>> &cells)
{

    // Set thread ranges.
    this->thread_loop_manager_.SetLoopRanges(v_nodal.size(), this->threads_number_);

    // Multithreaded reaction term computation.
    std::vector<std::thread> threads;
    if (this->AdaptiveReaction()) {
        for (std::size_t t = 0; t != this->threads_number_; ++t) {
            threads.emplace_back(std::thread(&Monodomain::ReactionCallbackAdaptive, this, t,
                                 std::cref(stimuli), step, std::ref(v_nodal), std::ref(cells)));
        }
    }
    else {
        for (std::size_t t = 0; t != this->threads_number_; ++t) {
            threads.emplace_back(std::thread(&Monodomain::ReactionCallbackStandard, this, t,
                                 std::cref(stimuli), step, std::ref(v_nodal), std::ref(cells)));
        }
    }
    // Join threads.
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ReactionCallbackStandard(std::size_t thread_id, const std::vector<Stimulus> &stimuli,
                                                           int step, Eigen::VectorXd &v_nodal, std::vector<std::unique_ptr<EpBasic>> &cells)
{

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

        // Update the time varying parameters of the cell's ep model if any.
        if (this->CellVarParamFlags()[i] != -1) {
            cells[i]->UpdateTimeVaryingPrms(this->CellVarParamGroups()[this->CellVarParamFlags()[i]], this->Dt());
        }

        // Compute new cell model state  and update the cell membrane potential.
        cells[i]->Compute(v_nodal.coeff(i), this->Dt(), stimulus_current);
        cells[i]->SetV(ALGORITHM::ForwardEuler(v_nodal.coeff(i), this->Dt(), cells[i]->dVdt()));

        // Update the nodal potential value.
        v_nodal.coeffRef(i) = cells[i]->V();

    } // End of Iterate over nodes (nodal cells) for reaction term.

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ReactionCallbackAdaptive(std::size_t thread_id, const std::vector<Stimulus> &stimuli,
    int step, Eigen::VectorXd &v_nodal, std::vector<std::unique_ptr<EpBasic>> &cells)
{

    // Compute maximum adaptive time step multiplication factor.
    int kappa_max = static_cast<int>(std::ceil(this->Dt() / this->DtMin()));

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
        double v_adapt = v_nodal.coeff(i);

        int safe_guard = 100000*kappa_max, steps_done = 0;
        while (t_cell < this->Dt()) {

            // Ensure that integration time is exactly equal to the total time step.
            if (t_cell + dt_adapt > this->Dt()) { dt_adapt = this->Dt() - t_cell; }

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
            dt_adapt = this->Dt() / kappa;

            steps_done++;
            if (steps_done > safe_guard) {
                std::string error_msg = "Error occured during adaptive integration of the reaction term due to NaN in the potential derivative computation.";
                throw std::runtime_error(Logger::Error(error_msg));
            }
        }

        // Set the final transmembrane potential value to the cell.
        cells[i]->SetV(v_adapt);

        // Update the nodal potential value.
        v_nodal.coeffRef(i) = v_adapt;

    } // End of Iterate over nodes (nodal cells) for reaction term.

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ComputeDiffusionExplicit(const Eigen::SparseMatrix<double, Eigen::RowMajor> &stiff_mat, const Eigen::VectorXd &mass_vec, 
                                                           double dt, Eigen::VectorXd &v_nodal)
{
    // Set thread ranges.
    this->thread_loop_manager_.SetLoopRanges(v_nodal.rows(), this->threads_number_);

    // Multithreaded diffusion term computation.
    std::vector<std::thread> threads;
    if (this->AdaptiveDiffusion() == true) {

        // Set adaptive time step.
        double adaptive_dt = ((dt < this->DtCritical()) ? dt : this->DtCritical());

        // Compute membrane potential diffusion adaptively.
        double total_time = 0.;
        // while (total_time+adaptive_dt <= dt) {
        while (dt - total_time > 2*std::numeric_limits<double>::epsilon()) {

            // Compute the diffusion propagation for the adaptive time step.
            threads.clear();
            for (std::size_t t = 0; t != this->threads_number_; ++t) {
                threads.emplace_back(std::thread(&Monodomain::DiffusionCallback, this, t, std::cref(stiff_mat),
                        std::cref(mass_vec), adaptive_dt, std::ref(v_nodal)));
            }

            // Join threads to ensure that the new potential values are updated by all threads.
            std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

            // Increase the processed total time.
            total_time += adaptive_dt;

            // Ensure that time integration is exact up to the predefined dt.
            if (dt - total_time < adaptive_dt)  { adaptive_dt = (dt - total_time); }

        } // End of Compute membrane potential diffusion adaptively.

        // // Ensure that time integration is exact up to the predefined dt.
        // if (total_time < (dt - 2*std::numeric_limits<double>::epsilon())) {
        //     // Set the adaptive time to the remaining time until the predefined dt.
        //     adaptive_dt = dt - total_time;

        //     // Compute the diffusion propagation for the last part of the time step.
        //     threads.clear();
        //     for (std::size_t t = 0; t != this->threads_number_; ++t) {
        //         threads.emplace_back(std::thread(&Monodomain::DiffusionCallback, this, t, std::cref(stiff_mat),
        //                                          std::cref(mass_vec), adaptive_dt, std::ref(v_nodal)));
        //     }

        //     // Join threads.
        //     std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
        // }
    }
    else {
        threads.clear();
        for (std::size_t t = 0; t != this->threads_number_; ++t) {
            threads.emplace_back(std::thread(&Monodomain::DiffusionCallback, this, t, std::cref(stiff_mat),
                                 std::cref(mass_vec), dt, std::ref(v_nodal)));
        }

        // Join threads.
        std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
    }

}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::ComputeDiffusionImplicit(const Eigen::SparseMatrix<double, Eigen::RowMajor> &stiff_mat, 
                                                           const Eigen::SparseMatrix<double, Eigen::RowMajor> &mass_mat, 
                                                           const Eigen::VectorXd &v_old, Eigen::VectorXd &v_new)
{
    // Solve implicitly the system (M - dt*K)*v_new = M*vold.
    // v_new = arma::spsolve((mass_mat - this->dt_*stiffness_mat), (mass_mat*v_old));
}


template<short DIM, short CELL_NODES>
void Monodomain<DIM, CELL_NODES>::DiffusionCallback(std::size_t thread_id, const Eigen::SparseMatrix<double, Eigen::RowMajor> &stiff_mat, 
                                                    const Eigen::VectorXd &mass_vec, double dt, Eigen::VectorXd &v_nodal)
{

    // Iterate over the nodes in the thread for diffusion computation.
    for (auto i = this->thread_loop_manager_.LoopStartId(thread_id);
            i != this->thread_loop_manager_.LoopEndId(thread_id); ++i) {
        double stiff_x_v_nodal = 0.;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(stiff_mat, i); it; ++it) {
            stiff_x_v_nodal += ( it.value() * v_nodal.coeff(it.col()) );
        }

        // Compute membrane potential difussion.
        v_nodal.coeffRef(i) = ALGORITHM::ForwardEuler(v_nodal.coeff(i), dt, -stiff_x_v_nodal/mass_vec.coeff(i));
    } // End Iterate over nodes (nodal cells) for diffusion term.
}


} // End of namespace ELECTRA


#endif //ELECTRA_ENGINE_PHYSICS_MONODOMAIN_TPP_