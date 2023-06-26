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

#ifndef ELECTRA_FIBERS_LDRBM_TPP_
#define ELECTRA_FIBERS_LDRBM_TPP_

#include "ELECTRA/engine/fibers/ldrbm.hpp"


namespace ELECTRA
{

template<short DIM, short CELL_NODES>
Ldrbm<DIM, CELL_NODES>::Ldrbm() : laplacian_mat_(), transmural_direction_(), appendage_veins_direction_(), inter_veins_direction_(),
    valve_veins_direction_(), atri_tricuspid_direction_(), apicobasal_direction_(), fiber_direction_(), sheet_direction_(), transmural_distance_(),
    appendage_veins_distance_(), inter_veins_distance_(), valve_veins_distance_(), atri_tricuspid_distance_(), septal_distance_(), apicobasal_distance_(),
    intraventricular_function_(), la_node_ids_(), ra_node_ids_(), lv_node_ids_(), rv_node_ids_(), ventricle_interface_node_ids_(), septum_node_ids_()
{}


template<short DIM, short CELL_NODES>
Ldrbm<DIM, CELL_NODES>::~Ldrbm()
{}


template<short DIM, short CELL_NODES>
Eigen::VectorXd Ldrbm<DIM, CELL_NODES>::SliceVector(const Eigen::VectorXd &eigen_vec, const std::vector<int> &row_ids)
{
    Eigen::VectorXd sliced(row_ids.size());

    int id = 0;
    for (const auto &row_id : row_ids) {
        sliced.coeffRef(id++) = eigen_vec.coeff(row_id);
    }

    return sliced;
}


template<short DIM, short CELL_NODES>
Eigen::SparseMatrix<double, Eigen::RowMajor> Ldrbm<DIM, CELL_NODES>::SliceSparseMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor> &eigen_spmat,
                                                                                       const std::vector<int> &row_ids, const std::vector<int> &col_ids)
{
    // Initialize triplets to slice the entries in the rows first.
    std::vector<Eigen::Triplet<double>> sliced_entries;
    sliced_entries.reserve(eigen_spmat.nonZeros());

    // Get non-zero entries in the rows of interest.
    int i = 0;
    for (const auto &row_id : row_ids) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(eigen_spmat, row_id); it; ++it) {
            sliced_entries.emplace_back(Eigen::Triplet<double>(i, it.col(), it.value()) );
        }
        i++;
    }

    // Construct the sparse matrix with sliced rows and column major order.
    Eigen::SparseMatrix<double> row_sliced_mat(row_ids.size(), eigen_spmat.cols());
    row_sliced_mat.setFromTriplets(std::begin(sliced_entries), std::end(sliced_entries));

    // Get non-zero entries in the colums of interest.
    int j = 0;
    sliced_entries.clear();
    for (const auto &col_id : col_ids) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(row_sliced_mat, col_id); it; ++it) {
            sliced_entries.emplace_back(Eigen::Triplet<double>(it.row(), j, it.value()) );
        }
        j++;
    }

    // Construct the final sparse matrix with sliced rows and columns and row major order.
    Eigen::SparseMatrix<double, Eigen::RowMajor> final_sliced_mat(row_ids.size(), col_ids.size());
    final_sliced_mat.setFromTriplets(std::begin(sliced_entries), std::end(sliced_entries));

    return final_sliced_mat;
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::SolveLaplacian(const std::unordered_map<std::string, IMP::NodeSet> &nsets, const std::vector<std::string> &dirichlet_tags, 
                                            const std::vector<double> &dirichlet_values, Eigen::VectorXd &solution)
{
    if (this->laplacian_mat_.rows() == 0) {
        std::string error_msg = "Could not solve Laplacian. Check that the laplacian matrix has been assembled correctly.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    //Initialize solution vector.
    int nodes_num = this->laplacian_mat_.rows();
    solution = Eigen::VectorXd::Zero(nodes_num);

    // Set free node flag indices and boundary conditions in solution vector.
    std::vector<short> free_node_flags(nodes_num, 1);
    int tag_id = 0;
    for (const auto &tag : dirichlet_tags) {
        for (const auto &nid : nsets.at(tag).NodeIds()) {
            free_node_flags[nid] = 0;
            solution.coeffRef(nid) = dirichlet_values[tag_id];
        }
        tag_id++;
    }

    // Collect free node indices.
    std::vector<int> free_node_ids;
    free_node_ids.reserve(nodes_num);
    int id = 0;
    for (const auto &flag : free_node_flags) {
        if (flag == 1) { free_node_ids.emplace_back(id); }
        id++;
    }
    free_node_flags.clear();
    free_node_flags.shrink_to_fit();

    // Compute body load vector.
    Eigen::VectorXd b = - this->laplacian_mat_*solution;

    // Extract laplacian matrix and body load vector of free nodes.
    auto b_free = this->SliceVector(b, free_node_ids);
    auto laplacian_free = this->SliceSparseMatrix(this->laplacian_mat_, free_node_ids, free_node_ids);

    // Solve laplacian for free nodes with BICGSTAB solver.
    Eigen::initParallel();
    Eigen::setNbThreads(std::thread::hardware_concurrency()-1);
    // std::cout << "Using threads: " << Eigen::nbThreads() << "\n";

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>> solver;
    solver.compute(laplacian_free);
    Eigen::VectorXd sol_free = solver.solve(b_free);

    // Update the solution vector for free nodes.
    id = 0;
    for (const auto &fnode_id : free_node_ids) {
        solution.coeffRef(fnode_id) = sol_free.coeff(id++);
    }
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeGradient(const IMP::Mesh<DIM,CELL_NODES> &mesh, const Eigen::VectorXd &nodal_scalars, Eigen::MatrixXd &nodal_gradients)
{
    // Set cell faces connectivity.
    Eigen::MatrixXd faces;
    if (CELL_NODES == 4) {
        faces.resize(4,3);
        faces << 0, 2, 1,
                 0, 3, 2,
                 1, 2, 3,
                 0, 1, 3;
    } else {
        std::string error_msg = "Could not compute gradient. Currently supports only tetrahedral meshes.";
        throw std::runtime_error(Logger::Error(error_msg));
    }
    // } else if (CELL_NODES == 8) {
    //     faces << 0, 3, 2, 1,
    //              4, 5, 6, 7,
    //              0, 1, 5, 4,
    //              2, 3, 7, 6,
    //              1, 2, 6, 5,
    //              0, 4, 7, 3;
    // }

    // Tetrahedron volume lambda.
    auto TetVolume = [](IMP::Vec<DIM,double> v0, IMP::Vec<DIM,double> v1,
                          IMP::Vec<DIM,double> v2, IMP::Vec<DIM,double> v3)
    {
        IMP::Vec<DIM,double> a = v0 - v3;
        IMP::Vec<DIM,double> b = v1 - v3;
        IMP::Vec<DIM,double> c = v2 - v3;

        double det = a[0] * (b[1]*c[2] - c[1]*b[2]) -
                     b[0] * (a[1]*c[2] - c[1]*a[2]) +
                     c[0] * (a[1]*b[2] - b[1]*a[2]);

        return std::abs(det)/6.;
    };

    // // Hexahedron volume lambda.
    // double HexVolume = [](IMP::Vec<DIM,double> v0, IMP::Vec<DIM,double> v1,
    //                       IMP::Vec<DIM,double> v2, IMP::Vec<DIM,double> v3,
    //                       IMP::Vec<DIM,double> v4, IMP::Vec<DIM,double> v5,
    //                       IMP::Vec<DIM,double> v6, IMP::Vec<DIM,double> v7)
    // { };

    // Compute gradient at cell centroids.
    Eigen::MatrixXd cell_gradients = Eigen::MatrixXd::Zero(mesh.CellsNum(),DIM);
    Eigen::VectorXd cell_grad;
    std::vector<IMP::Vec<DIM,double>> cell_centroids(mesh.CellsNum());
    std::vector<std::vector<int>> attached_cells_to_nodes(mesh.NodesNum());
    IMP::Vec<DIM,double> face_centroid, face_normal;
    double cell_vol = 0., face_scalar = 0.;
    int cid = 0;
    for (const auto &cell : mesh.Cells()) {

        // Compute the centroid of the cell.
        cell_centroids[cid].SetZero();
        for (const auto &nid : cell.Connectivity()) {
            // Add the cell to the container of attached cell to each of its nodes.
            attached_cells_to_nodes[nid].emplace_back(cid);
            cell_centroids[cid] += mesh.Nodes(nid);
        }
        cell_centroids[cid] /= 4.;

        // Compute the volume of the cell.
        cell_vol = TetVolume(mesh.Nodes(cell.N(0)), mesh.Nodes(cell.N(1)),
                             mesh.Nodes(cell.N(2)), mesh.Nodes(cell.N(3)));

        // Compute the gradient at the cell's centroid.
        cell_grad = Eigen::VectorXd::Zero(DIM);
        for (int f=0; f!=faces.rows(); ++f) {

            // Compute face normal vector with Newell's method.
            face_normal.SetZero();
            face_centroid.SetZero();
            face_scalar = 0.;
            for (int i=0; i!=faces.cols(); ++i) {
                auto v = mesh.Nodes( cell.N(faces(f,i)) );
                auto w = mesh.Nodes( cell.N(faces(f,(i+1)%faces.cols())) );

                // Compute face normal by cross product of sequential nodes
                face_normal[0] += v[1]*w[2] - v[2]*w[1];
                face_normal[1] += v[2]*w[0] - v[0]*w[2];
                face_normal[2] += v[0]*w[1] - v[1]*w[0];

                // Compute face centroid.
                face_centroid += mesh.Nodes( cell.N(faces(f,i)) );

                // Compute scalar at face centroid.
                face_scalar += nodal_scalars.coeffRef( cell.N(faces(f,i)) );
            }

            // Normalize face centroid and face scalar.
            face_centroid /= static_cast<double>(faces.cols());
            face_scalar /= static_cast<double>(faces.cols());

            // Check face normal direction.
            if (face_normal.CwiseMul(face_centroid-cell_centroids[cid]).Sum() < 0)  face_normal = -face_normal;

            // Add the gradient contribution from the face centroid.
            cell_grad += face_scalar*0.5*face_normal.CopyToEigen();
        }
        cell_gradients.row(cid++) = cell_grad / cell_vol;

    } // End of Iterate over the mesh cells.

    // Compute nodal gradients.
    int nid = 0;
    double weight = 0., weights_sum = 0.;
    Eigen::RowVectorXd nodal_grad = Eigen::RowVectorXd::Zero(DIM);
    nodal_gradients = Eigen::MatrixXd::Zero(nodal_scalars.rows(),DIM);
    for (const auto &attached_cells : attached_cells_to_nodes) {
        weights_sum = 0.;
        nodal_grad = Eigen::RowVectorXd::Zero(DIM);
        for (const auto &cid : attached_cells) {
            weight = 1. / std::sqrt(mesh.Nodes(nid).Distance2(cell_centroids[cid]));
            weights_sum += weight;
            nodal_grad += weight*cell_gradients.row(cid);
        }
        nodal_gradients.row(nid++) = nodal_grad/weights_sum;
    }
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeGradient(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const Eigen::VectorXd &nodal_scalars, Eigen::MatrixXd &nodal_gradients)
{
    // Compute nodal gradients.
    std::vector<int> neighs;
    Eigen::RowVectorXd neighs_scalar;
    nodal_gradients = Eigen::MatrixXd::Zero(nodal_scalars.rows(),DIM);
    for (int nid = 0; nid != voro.NodesNum(); ++nid) {

        // Get neighbor nodes.
        neighs = fpm.Support().InfluenceNodeIds(nid);

        // Collect the nodal scalar values of the neighbor nodes.
        neighs_scalar = Eigen::RowVectorXd::Zero(neighs.size());
        for (std::size_t i = 0; i != neighs.size(); ++i) {
            neighs_scalar.coeffRef(i) = nodal_scalars.coeff(neighs[i]);
        }

        // Compute nodal gradient.
        nodal_gradients.row(nid) = neighs_scalar*fpm.PhiGrad(nid);
    }
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::AssembleLaplacian(const IMP::Mesh<DIM,CELL_NODES> &mesh)
{
    // Create a finite element of the given type.
    std::unique_ptr<CLOUDEA::BaseFem<DIM>> elem = CLOUDEA::FemFactory<DIM>::Create(mesh.CellsShape());

    // Set the quadrature of the finite element.
    switch (mesh.CellsShape()) {
        case IMP::CellShape::tet :
            elem->SetQuadrature({4});
        break;
        case IMP::CellShape::hex :
            elem->SetQuadrature({2,2,2});
        break;
        default:
            std::string error_msg = "Could not assemble laplacian matrix. A 3D mesh with either tetrahedral or hexahedral elements is required.";
            throw std::invalid_argument(Logger::Error(error_msg));
        break;
    }

    // Compute natural derivatives and shape functions of the element.
    elem->ComputeDerivsNatural();
    elem->ComputeShapeFunctions();

    // Initialize triplets for assembly.
    std::vector<Eigen::Triplet<double>> laplacian_triplets;
    laplacian_triplets.reserve(mesh.NodesNum() * CELL_NODES);

    // Initialize container to store the coordinates of the cell nodes.
    std::vector<IMP::Vec<DIM, double>> nodes(CELL_NODES);

    // Iterate over the mesh cells.
    for (const auto &cell : mesh.Cells()) {
        // Get coordinates from mesh nodes container for nD nodes in nD space.
        for (const auto &nid : cell.Connectivity()) {
            auto i = &nid - &cell.Connectivity()[0];
            for (short d = 0; d != DIM; ++d) { nodes[i][d] = mesh.Nodes(nid)[d]; }
        }

        // Compute the jacobian and physical derivatives for the current cell.
        elem->ComputeJacobians(nodes);
        elem->ComputeDerivs();

        // Construct the local laplacian operator over the element's quadrature.
        int qid = 0;
        Eigen::MatrixXd  laplacian_cell = Eigen::MatrixXd::Zero(CELL_NODES, CELL_NODES);
        for (const auto &qweight : elem->Quadrature().Weights()) {
            laplacian_cell += qweight*elem->DetJacobians()[qid] * (elem->Derivs()[qid]*elem->Derivs()[qid].transpose());
            qid++;
        }

        // Assemble the local laplacian operator of the element.
        int gi = 0, gj = 0;
        for (short li = 0; li != CELL_NODES; ++li) {
            gi = cell.N(li);
            for (short lj = 0; lj != CELL_NODES; ++lj) {
                gj = cell.N(lj);
                laplacian_triplets.emplace_back(Eigen::Triplet<double>(gi, gj, laplacian_cell(li, lj)));
            }
        }

    } // End of Iterate over the mesh cells.
    laplacian_triplets.shrink_to_fit();

    // Set cardiac stiffness matrix.
    this->laplacian_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(mesh.NodesNum(), mesh.NodesNum());
    this->laplacian_mat_.setFromTriplets(std::begin(laplacian_triplets), std::end(laplacian_triplets));
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::AssembleLaplacian(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, double penalty)
{
    if (static_cast<int>(fpm.Phi().size()) != voro.NodesNum()) {
        std::string error_str = "Could not assemble laplacian matrix. Fpm approximation is not available for all voronoi tesselation nodes.";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    // Initialize triplets for assembly.
    int n = voro.NodesNum();
    std::vector<Eigen::Triplet<double>> laplacian_triplets;
    laplacian_triplets.reserve(n);

    // Assemble stiffness matrix.
    Eigen::MatrixXd stiff_cell;
    for (const auto &cell : voro.Cells()) {
        auto id = &cell - &voro.Cells()[0];

        // Compute laplacian matrix of the cell.
        stiff_cell = fpm.PhiGrad(id)*fpm.PhiGrad(id).transpose()*cell.Measure(voro.Points(), voro.Facets());

        // Assemble matrices.
        for (const auto &gi : fpm.Support().InfluenceNodeIds(id)) {
            auto li = &gi - &fpm.Support().InfluenceNodeIds(id)[0];

            for (const auto &gj : fpm.Support().InfluenceNodeIds(id)) {
                auto lj = &gj - &fpm.Support().InfluenceNodeIds(id)[0];
                laplacian_triplets.emplace_back( Eigen::Triplet<double>(gi, gj, stiff_cell(li, lj)) );
            }
        }
    } // End of Assemble stiffness matrix.

    // Set cardiac stiffness matrix.
    this->laplacian_mat_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n, n);
    this->laplacian_mat_.setFromTriplets(std::begin(laplacian_triplets), std::end(laplacian_triplets));
    laplacian_triplets.clear(); laplacian_triplets.shrink_to_fit();

    // Apply numerical correction at internal facets.
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

    // Initialize shape function derivatives and correction matrices.
    Eigen::MatrixXd B1, B2, C11, C12, C21, C22;

    // Initialize correction triplets vector and correction matrix.
    std::vector<Eigen::Triplet<double>> corr_triplets;
    corr_triplets.reserve(entries_num);
    auto corr_mat = Eigen::SparseMatrix<double, Eigen::RowMajor>(n, n);

    // Iterate facets in chunks.
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
                    C11 = area * (-0.5*(N1*ne1.transpose()*B1.transpose() + B1*ne1*N1.transpose()) + (penalty/he)*(N1*N1.transpose()));
                    C12 = area * (-0.5*(N1*ne1.transpose()*B2.transpose() + B1*ne2*N2.transpose()) - (penalty/he)*(N1*N2.transpose()));
                    C21 = area * (-0.5*(N2*ne2.transpose()*B1.transpose() + B2*ne1*N1.transpose()) - (penalty/he)*(N2*N1.transpose()));
                    C22 = area * (-0.5*(N2*ne2.transpose()*B2.transpose() + B2*ne2*N2.transpose()) + (penalty/he)*(N2*N2.transpose()));

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

        // Add laplacian matrix matrix correction.
        corr_mat.setFromTriplets(std::begin(corr_triplets), std::end(corr_triplets));
        this->laplacian_mat_ += corr_mat;

        // Clean correction triplets and matrix.
        corr_triplets.clear();
        corr_mat.setZero();

    } // End iterate facets in chunks.

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
                    C11 = area * (-0.5*(N1*ne1.transpose()*B1.transpose() + B1*ne1*N1.transpose()) + (penalty/he)*(N1*N1.transpose()));
                    C12 = area * (-0.5*(N1*ne1.transpose()*B2.transpose() + B1*ne2*N2.transpose()) - (penalty/he)*(N1*N2.transpose()));
                    C21 = area * (-0.5*(N2*ne2.transpose()*B1.transpose() + B2*ne1*N1.transpose()) - (penalty/he)*(N2*N1.transpose()));
                    C22 = area * (-0.5*(N2*ne2.transpose()*B2.transpose() + B2*ne2*N2.transpose()) + (penalty/he)*(N2*N2.transpose()));

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
        this->laplacian_mat_ += corr_mat;
    }

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriTransmural(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags)
{
    // Check that epiwall and at least one endowall tag are defined.
    if ((atri_tags.LaEndoWallTag().empty() && atri_tags.RaEndoWallTag().empty()) || atri_tags.EpiWallTag().empty()) {
        std::string error_msg = "Could not compute fibers transmural distance. Epicardium and endocardium wall tags should be defined.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Set epicardium boundary tags and values.
    std::vector<std::string> dirichlet_tags{atri_tags.EpiWallTag()};
    std::vector<double> dirichlet_values{0.};

    // Add left atrium endocardium boundary tags and values.
    if (!atri_tags.LaEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(atri_tags.LaEndoWallTag());
        dirichlet_values.emplace_back(-1.);
    }

    // Add right atrium endocardium boundary tags and values.
    if (!atri_tags.RaEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(atri_tags.RaEndoWallTag());
        dirichlet_values.emplace_back(1.);
    }

    // Compute transmural distance and gradient.
    this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, this->transmural_distance_);
    this->ComputeGradient(mesh, this->transmural_distance_, this->transmural_direction_);

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriTransmural(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags)
{
    // Check that epiwall and at least one endowall tag are defined.
    if ((atri_tags.LaEndoWallTag().empty() && atri_tags.RaEndoWallTag().empty()) || atri_tags.EpiWallTag().empty()) {
        std::string error_msg = "Could not compute fibers transmural distance. Epicardium and endocardium wall tags should be defined.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Set epicardium boundary tags and values.
    std::vector<std::string> dirichlet_tags{atri_tags.EpiWallTag()};
    std::vector<double> dirichlet_values{0.};

    // Add left atrium endocardium boundary tags and values.
    if (!atri_tags.LaEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(atri_tags.LaEndoWallTag());
        dirichlet_values.emplace_back(-1.);
    }

    // Add right atrium endocardium boundary tags and values.
    if (!atri_tags.RaEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(atri_tags.RaEndoWallTag());
        dirichlet_values.emplace_back(1.);
    }

    // Compute transmural distance and gradient.
    this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, this->transmural_distance_);
    this->ComputeGradient(voro, fpm, this->transmural_distance_, this->transmural_direction_);
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriAppendageVeins(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags)
{
    // Check if at least one of the atria has been extracted.
    if (this->ra_node_ids_.empty() && this->la_node_ids_.empty()) {
        std::string error_msg = "Could not compute atria appendage to veins distance and direction. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize appendage to veins distance and direction.
    this->appendage_veins_distance_ = Eigen::VectorXd::Zero(mesh.NodesNum());
    this->appendage_veins_direction_ = Eigen::MatrixXd::Zero(mesh.NodesNum(),DIM);

    // Compute appendage to veins distance and direction at left atrium.
    if (!this->la_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.RightPulmonVeinTag(), atri_tags.LeftPulmonVeinTag(),
                                                atri_tags.MitralValveTag(), atri_tags.LaAppendageTag()};
        std::vector<double> dirichlet_values{2., 0., 1., -1.};

        // Solve appendage to veins component.
        Eigen::VectorXd psi_ab;
        Eigen::MatrixXd psi_ab_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_ab);
        this->ComputeGradient(mesh, psi_ab, psi_ab_grad);

        for (const auto &la_nid : this->la_node_ids_) {
            this->appendage_veins_distance_.coeffRef(la_nid) = psi_ab.coeff(la_nid);
            this->appendage_veins_direction_.row(la_nid) = psi_ab_grad.row(la_nid);
        }
    } // End of Compute appendage to veins distance and direction at left atrium.


    // Compute appendage to veins distance and direction at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.InferiorCavalVeinTag(), atri_tags.SuperiorCavalVeinTag(),
                                                atri_tags.TricuspidValveSeptumTag(), atri_tags.TricuspidValveFreeTag(), atri_tags.RaAppendageTag()};
        std::vector<double> dirichlet_values{2., 0., 1., 1., -1.};

        // Solve appendage to veins component.
        Eigen::VectorXd psi_ab;
        Eigen::MatrixXd psi_ab_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_ab);
        this->ComputeGradient(mesh, psi_ab, psi_ab_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->appendage_veins_distance_.coeffRef(ra_nid) = psi_ab.coeff(ra_nid);
            this->appendage_veins_direction_.row(ra_nid) = psi_ab_grad.row(ra_nid);
        }

    } // End of Compute appendage to veins distance and direction at right ventricle.

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriAppendageVeins(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags)
{
    // Check if at least one of the atria has been extracted.
    if (this->ra_node_ids_.empty() && this->la_node_ids_.empty()) {
        std::string error_msg = "Could not compute atria appendage to veins distance and direction. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize appendage to veins distance and direction.
    this->appendage_veins_distance_ = Eigen::VectorXd::Zero(voro.NodesNum());
    this->appendage_veins_direction_ = Eigen::MatrixXd::Zero(voro.NodesNum(),DIM);

    // Compute appendage to veins distance and direction at left atrium.
    if (!this->la_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.RightPulmonVeinTag(), atri_tags.LeftPulmonVeinTag(),
                                                atri_tags.MitralValveTag(), atri_tags.LaAppendageTag()};
        std::vector<double> dirichlet_values{2., 0., 1., -1.};

        // Solve appendage to veins component.
        Eigen::VectorXd psi_ab;
        Eigen::MatrixXd psi_ab_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_ab);
        this->ComputeGradient(voro, fpm, psi_ab, psi_ab_grad);

        for (const auto &la_nid : this->la_node_ids_) {
            this->appendage_veins_distance_.coeffRef(la_nid) = psi_ab.coeff(la_nid);
            this->appendage_veins_direction_.row(la_nid) = psi_ab_grad.row(la_nid);
        }
    } // End of Compute appendage to veins distance and direction at left atrium.


    // Compute appendage to veins distance and direction at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.InferiorCavalVeinTag(), atri_tags.SuperiorCavalVeinTag(),
                                                atri_tags.TricuspidValveSeptumTag(), atri_tags.TricuspidValveFreeTag(), atri_tags.RaAppendageTag()};
        std::vector<double> dirichlet_values{2., 0., 1., 1., -1.};

        // Solve appendage to veins component.
        Eigen::VectorXd psi_ab;
        Eigen::MatrixXd psi_ab_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_ab);
        this->ComputeGradient(voro, fpm, psi_ab, psi_ab_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->appendage_veins_distance_.coeffRef(ra_nid) = psi_ab.coeff(ra_nid);
            this->appendage_veins_direction_.row(ra_nid) = psi_ab_grad.row(ra_nid);
        }

    } // End of Compute appendage to veins distance and direction at right ventricle.
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriInterVeins(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags)
{
    // Check if at least one of the atria has been extracted.
    if (this->ra_node_ids_.empty() && this->la_node_ids_.empty()) {
        std::string error_msg = "Could not compute atria interveins distance and direction. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize distance and direction between the atrial veins.
    this->inter_veins_distance_ = Eigen::VectorXd::Zero(mesh.NodesNum());
    this->inter_veins_direction_ = Eigen::MatrixXd::Zero(mesh.NodesNum(),DIM);

    // Compute distance and direction between the atrial veins at left atrium.
    if (!this->la_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.RightPulmonVeinTag(), atri_tags.LeftPulmonVeinTag()};
        std::vector<double> dirichlet_values{1., 0.};

        // Solve distance and direction between the atrial veins component.
        Eigen::VectorXd psi_v;
        Eigen::MatrixXd psi_v_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_v);
        this->ComputeGradient(mesh, psi_v, psi_v_grad);

        for (const auto &la_nid : this->la_node_ids_) {
            this->inter_veins_distance_.coeffRef(la_nid) = psi_v.coeff(la_nid);
            this->inter_veins_direction_.row(la_nid) = psi_v_grad.row(la_nid);
        }
    } // End of Compute distance and direction between the atrial veins at left atrium.


    // Compute distance and direction between the atrial veins at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.InferiorCavalVeinTag(), atri_tags.SuperiorCavalVeinTag(), atri_tags.RaAppendageTag()};
        std::vector<double> dirichlet_values{1., 0., 0.};

        // Solve distance and direction between the atrial veins component.
        Eigen::VectorXd psi_v;
        Eigen::MatrixXd psi_v_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_v);
        this->ComputeGradient(mesh, psi_v, psi_v_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->inter_veins_distance_.coeffRef(ra_nid) = psi_v.coeff(ra_nid);
            this->inter_veins_direction_.row(ra_nid) = psi_v_grad.row(ra_nid);
        }

    } // End of Compute distance and direction between the atrial veins at right ventricle.

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriInterVeins(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags)
{
    // Check if at least one of the atria has been extracted.
    if (this->ra_node_ids_.empty() && this->la_node_ids_.empty()) {
        std::string error_msg = "Could not compute atria interveins distance and direction. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize distance and direction between the atrial veins.
    this->inter_veins_distance_ = Eigen::VectorXd::Zero(voro.NodesNum());
    this->inter_veins_direction_ = Eigen::MatrixXd::Zero(voro.NodesNum(),DIM);

    // Compute distance and direction between the atrial veins at left atrium.
    if (!this->la_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.RightPulmonVeinTag(), atri_tags.LeftPulmonVeinTag()};
        std::vector<double> dirichlet_values{1., 0.};

        // Solve distance and direction between the atrial veins component.
        Eigen::VectorXd psi_v;
        Eigen::MatrixXd psi_v_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_v);
        this->ComputeGradient(voro, fpm, psi_v, psi_v_grad);

        for (const auto &la_nid : this->la_node_ids_) {
            this->inter_veins_distance_.coeffRef(la_nid) = psi_v.coeff(la_nid);
            this->inter_veins_direction_.row(la_nid) = psi_v_grad.row(la_nid);
        }
    } // End of Compute distance and direction between the atrial veins at left atrium.


    // Compute distance and direction between the atrial veins at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.InferiorCavalVeinTag(), atri_tags.SuperiorCavalVeinTag(), atri_tags.RaAppendageTag()};
        std::vector<double> dirichlet_values{1., 0., 0.};

        // Solve distance and direction between the atrial veins component.
        Eigen::VectorXd psi_v;
        Eigen::MatrixXd psi_v_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_v);
        this->ComputeGradient(voro, fpm, psi_v, psi_v_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->inter_veins_distance_.coeffRef(ra_nid) = psi_v.coeff(ra_nid);
            this->inter_veins_direction_.row(ra_nid) = psi_v_grad.row(ra_nid);
        }

    } // End of Compute distance and direction between the atrial veins at right ventricle.
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriValveVeins(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags)
{
    // Check if at least one of the atria has been extracted.
    if (this->ra_node_ids_.empty() && this->la_node_ids_.empty()) {
        std::string error_msg = "Could not compute atria valve to veins distance and direction. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize valve to veins distance and direction.
    this->valve_veins_distance_ = Eigen::VectorXd::Zero(mesh.NodesNum());
    this->valve_veins_direction_ = Eigen::MatrixXd::Zero(mesh.NodesNum(),DIM);

    // Compute valve to veins distance and direction at left atrium.
    if (!this->la_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.MitralValveTag(), atri_tags.RightPulmonVeinTag(),
                                                atri_tags.LeftPulmonVeinTag(), atri_tags.LaAppendageTag()};
        std::vector<double> dirichlet_values{1., 0., 0., 0.};

        // Solve valve to veins component.
        Eigen::VectorXd psi_r;
        Eigen::MatrixXd psi_r_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_r);
        this->ComputeGradient(mesh, psi_r, psi_r_grad);

        for (const auto &la_nid : this->la_node_ids_) {
            this->valve_veins_distance_.coeffRef(la_nid) = psi_r.coeff(la_nid);
            this->valve_veins_direction_.row(la_nid) = psi_r_grad.row(la_nid);
        }
    } // End of Compute valve to veins distance and direction at left atrium.


    // Compute valve to veins distance and direction at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.TricuspidValveSeptumTag(), atri_tags.TricuspidValveFreeTag(),
                                                atri_tags.RaTopEndoTag(), atri_tags.RaTopEpiTag()};
        std::vector<double> dirichlet_values{1., 1., 0., 0.};

        // Solve valve to veins component.
        Eigen::VectorXd psi_r;
        Eigen::MatrixXd psi_r_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_r);
        this->ComputeGradient(mesh, psi_r, psi_r_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->valve_veins_distance_.coeffRef(ra_nid) = psi_r.coeff(ra_nid);
            this->valve_veins_direction_.row(ra_nid) = psi_r_grad.row(ra_nid);
        }

    } // End of Compute appendage to veins distance and direction at right ventricle.

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriValveVeins(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags)
{
    // Check if at least one of the atria has been extracted.
    if (this->ra_node_ids_.empty() && this->la_node_ids_.empty()) {
        std::string error_msg = "Could not compute atria valve to veins distance and direction. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize valve to veins distance and direction.
    this->valve_veins_distance_ = Eigen::VectorXd::Zero(voro.NodesNum());
    this->valve_veins_direction_ = Eigen::MatrixXd::Zero(voro.NodesNum(),DIM);

    // Compute valve to veins distance and direction at left atrium.
    if (!this->la_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.MitralValveTag(), atri_tags.RightPulmonVeinTag(),
                                                atri_tags.LeftPulmonVeinTag(), atri_tags.LaAppendageTag()};
        std::vector<double> dirichlet_values{1., 0., 0., 0.};

        // Solve valve to veins component.
        Eigen::VectorXd psi_r;
        Eigen::MatrixXd psi_r_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_r);
        this->ComputeGradient(voro, fpm, psi_r, psi_r_grad);

        for (const auto &la_nid : this->la_node_ids_) {
            this->valve_veins_distance_.coeffRef(la_nid) = psi_r.coeff(la_nid);
            this->valve_veins_direction_.row(la_nid) = psi_r_grad.row(la_nid);
        }
    } // End of Compute valve to veins distance and direction at left atrium.


    // Compute valve to veins distance and direction at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.TricuspidValveSeptumTag(), atri_tags.TricuspidValveFreeTag(),
                                                atri_tags.RaTopEndoTag(), atri_tags.RaTopEpiTag()};
        std::vector<double> dirichlet_values{1., 1., 0., 0.};

        // Solve valve to veins component.
        Eigen::VectorXd psi_r;
        Eigen::MatrixXd psi_r_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_r);
        this->ComputeGradient(voro, fpm, psi_r, psi_r_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->valve_veins_distance_.coeffRef(ra_nid) = psi_r.coeff(ra_nid);
            this->valve_veins_direction_.row(ra_nid) = psi_r_grad.row(ra_nid);
        }

    } // End of Compute appendage to veins distance and direction at right ventricle.

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriTricuspid(const IMP::Mesh<DIM,CELL_NODES> &mesh, const AtriTags &atri_tags)
{
    // Check if the right atrium has been extracted.
    if (this->ra_node_ids_.empty()) {
        std::string error_msg = "Could not compute distance and direction between the free and septum parts of the atrial tricuspid valve. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize tricuspid septum to free part distance and direction.
    this->atri_tricuspid_distance_ = Eigen::VectorXd::Zero(mesh.NodesNum());
    this->atri_tricuspid_direction_ = Eigen::MatrixXd::Zero(mesh.NodesNum(),DIM);

    // Compute tricuspid septum to free part distance and direction at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.TricuspidValveSeptumTag(), atri_tags.TricuspidValveFreeTag()};
        std::vector<double> dirichlet_values{1., -1.};

        // Solve tricuspid septum to free part component.
        Eigen::VectorXd psi_w;
        Eigen::MatrixXd psi_w_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_w);
        this->ComputeGradient(mesh, psi_w, psi_w_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->atri_tricuspid_distance_.coeffRef(ra_nid) = psi_w.coeff(ra_nid);
            this->atri_tricuspid_direction_.row(ra_nid) = psi_w_grad.row(ra_nid);
        }
    } // End of tricuspid septum to free part distance and direction at right atrium.
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriTricuspid(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const AtriTags &atri_tags)
{
    // Check if the right atrium has been extracted.
    if (this->ra_node_ids_.empty()) {
        std::string error_msg = "Could not compute distance and direction between the free and septum parts of the atrial tricuspid valve. Extract atrial node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize tricuspid septum to free part distance and direction.
    this->atri_tricuspid_distance_ = Eigen::VectorXd::Zero(voro.NodesNum());
    this->atri_tricuspid_direction_ = Eigen::MatrixXd::Zero(voro.NodesNum(),DIM);

    // Compute tricuspid septum to free part distance and direction at right atrium.
    if (!this->ra_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{atri_tags.TricuspidValveSeptumTag(), atri_tags.TricuspidValveFreeTag()};
        std::vector<double> dirichlet_values{1., -1.};

        // Solve tricuspid septum to free part component.
        Eigen::VectorXd psi_w;
        Eigen::MatrixXd psi_w_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_w);
        this->ComputeGradient(voro, fpm, psi_w, psi_w_grad);

        for (const auto &ra_nid : this->ra_node_ids_) {
            this->atri_tricuspid_distance_.coeffRef(ra_nid) = psi_w.coeff(ra_nid);
            this->atri_tricuspid_direction_.row(ra_nid) = psi_w_grad.row(ra_nid);
        }
    } // End of tricuspid septum to free part distance and direction at right atrium.
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriTransmural(const IMP::Mesh<DIM,CELL_NODES> &mesh, const VentriTags &ventri_tags)
{
    // Check that epiwall and at least one endowall tag are defined.
    if ((ventri_tags.LvEndoWallTag().empty() && ventri_tags.RvEndoWallTag().empty()) || ventri_tags.EpiWallTag().empty()) {
        std::string error_msg = "Could not compute fibers transmural distance. Epicardium and endocardium wall tags should be defined.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Set epicardium boundary tags and values.
    std::vector<std::string> dirichlet_tags{ventri_tags.EpiWallTag()};
    std::vector<double> dirichlet_values{0.};

    // Add left ventricle endocardium boundary tags and values.
    if (!ventri_tags.LvEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(ventri_tags.LvEndoWallTag());
        dirichlet_values.emplace_back(-2.);
    }

    // Add right ventricle endocardium boundary tags and values.
    if (!ventri_tags.RvEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(ventri_tags.RvEndoWallTag());
        dirichlet_values.emplace_back(1.);
    }

    // Compute transmural distance and gradient.
    this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, this->transmural_distance_);
    this->ComputeGradient(mesh, this->transmural_distance_, this->transmural_direction_);
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriTransmural(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const VentriTags &ventri_tags)
{
    // Check that epiwall and at least one endowall tag are defined.
    if ((ventri_tags.LvEndoWallTag().empty() && ventri_tags.RvEndoWallTag().empty()) || ventri_tags.EpiWallTag().empty()) {
        std::string error_msg = "Could not compute fibers transmural distance. Epicardium and endocardium wall tags should be defined.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Set epicardium boundary tags and values.
    std::vector<std::string> dirichlet_tags{ventri_tags.EpiWallTag()};
    std::vector<double> dirichlet_values{0.};

    // Add left ventricle endocardium boundary tags and values.
    if (!ventri_tags.LvEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(ventri_tags.LvEndoWallTag());
        dirichlet_values.emplace_back(-2.);
    }

    // Add right ventricle endocardium boundary tags and values.
    if (!ventri_tags.RvEndoWallTag().empty()) {
        dirichlet_tags.emplace_back(ventri_tags.RvEndoWallTag());
        dirichlet_values.emplace_back(1.);
    }

    // Compute transmural distance and gradient.
    this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, this->transmural_distance_);
    this->ComputeGradient(voro, fpm, this->transmural_distance_, this->transmural_direction_);
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriSeptal(const IMP::Mesh<DIM, CELL_NODES> &mesh, const VentriTags &ventri_tags, double threshold)
{
    // Check that endowall tags are defined.
    if ((ventri_tags.LvEndoWallTag().empty() || ventri_tags.RvEndoWallTag().empty())) {
        std::string error_msg = "Could not compute fibers septal distance. Endocardium wall tags should be defined.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Check that ventricle interface nodes are defined.
    if (this->ventricle_interface_node_ids_.empty()) {
        std::string error_msg = "Could not compute fibers septal distance. Ventricle interface nodes are not available. You should extract ventricle nodal partitions first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Edit mesh node sets
    auto nsets = mesh.NodeSets();
    IMP::NodeSet ventri_interface;
    ventri_interface.Set("interface", this->ventricle_interface_node_ids_);
    nsets["interface"] = ventri_interface;

    // Set boundary tags and values.
    std::vector<std::string> dirichlet_tags{"interface", ventri_tags.LvEndoWallTag(), ventri_tags.RvEndoWallTag()};
    std::vector<double> dirichlet_values{1., 0., 0.};

    // Compute septal distance and gradient.
    this->SolveLaplacian(nsets, dirichlet_tags, dirichlet_values, this->septal_distance_);

    // Set septum node indices.
    this->septum_node_ids_.clear(); this->septum_node_ids_.reserve(mesh.NodesNum());
    for (Eigen::Index i = 0; i != this->septal_distance_.rows(); ++i) {
        if (this->septal_distance_.coeff(i) > threshold) { this->septum_node_ids_.emplace_back(i); }
    }
    this->septum_node_ids_.shrink_to_fit();
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriSeptal(const IMP::Voronoi<DIM> &voro, const VentriTags &ventri_tags, double threshold)
{
    // Check that endowall tags are defined.
    if ((ventri_tags.LvEndoWallTag().empty() || ventri_tags.RvEndoWallTag().empty())) {
        std::string error_msg = "Could not compute fibers septal distance. Endocardium wall tags should be defined.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Check that ventricle interface nodes are defined.
    if (this->ventricle_interface_node_ids_.empty()) {
        std::string error_msg = "Could not compute fibers septal distance. Ventricle interface nodes are not available. You should extract ventricle nodal partitions first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Edit mesh node sets
    auto nsets = voro.NodeSets();
    IMP::NodeSet ventri_interface;
    ventri_interface.Set("interface", this->ventricle_interface_node_ids_);
    nsets["interface"] = ventri_interface;

    // Set boundary tags and values.
    std::vector<std::string> dirichlet_tags{"interface", ventri_tags.LvEndoWallTag(), ventri_tags.RvEndoWallTag()};
    std::vector<double> dirichlet_values{1., 0., 0.};

    // Compute septal distance and gradient.
    this->SolveLaplacian(nsets, dirichlet_tags, dirichlet_values, this->septal_distance_);

    // Set septum node indices.
    this->septum_node_ids_.clear(); this->septum_node_ids_.reserve(voro.NodesNum());
    for (Eigen::Index i = 0; i != this->septal_distance_.rows(); ++i) {
        if (this->septal_distance_.coeff(i) > threshold) { this->septum_node_ids_.emplace_back(i); }
    }
    this->septum_node_ids_.shrink_to_fit();
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriApicoBasal(const IMP::Mesh<DIM,CELL_NODES> &mesh, const VentriTags &ventri_tags)
{
    // Check if at least one of the ventricles has been extracted.
    if (this->rv_node_ids_.empty() && this->lv_node_ids_.empty()) {
        std::string error_msg = "Could not compute apicobasal distance and direction. Extract ventricular node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Check if intraventricular interpolation function has been computed.
    if (this->intraventricular_function_.rows() == 0) {
        std::string error_msg = "Could not compute apicobasal distance and direction. Compute intraventricular interpolation function first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize apicobasal distance.
    this->apicobasal_distance_ = Eigen::VectorXd::Zero(mesh.NodesNum());
    this->apicobasal_direction_ = Eigen::MatrixXd::Zero(mesh.NodesNum(),DIM);

    // Compute apicobasal distance and direction at left ventricle.
    if (!this->lv_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.LvApexTag(), ventri_tags.MitralValveTag()};
        std::vector<double> dirichlet_values{1., 0.};

        // Solve apico-mitral component.
        Eigen::VectorXd psi_am;
        Eigen::MatrixXd psi_am_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_am);
        this->ComputeGradient(mesh, psi_am, psi_am_grad);

        // Solve apico-aortic component if available.
        if (!ventri_tags.AorticValveTag().empty()) {
            Eigen::VectorXd psi_aa;
            Eigen::MatrixXd psi_aa_grad;
            dirichlet_tags[1] = ventri_tags.AorticValveTag();
            this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_aa);
            this->ComputeGradient(mesh, psi_aa, psi_aa_grad);

            for (const auto &lv_nid : this->lv_node_ids_) {
                this->apicobasal_distance_.coeffRef(lv_nid) = this->intraventricular_function_.coeff(lv_nid)*psi_am.coeff(lv_nid) + (1.0-this->intraventricular_function_.coeff(lv_nid))*psi_aa.coeff(lv_nid);
                this->apicobasal_direction_.row(lv_nid) = this->intraventricular_function_.coeff(lv_nid)*psi_am_grad.row(lv_nid) + (1.0-this->intraventricular_function_.coeff(lv_nid))*psi_aa_grad.row(lv_nid);
            }
        } else {
            for (const auto &lv_nid : this->lv_node_ids_) {
                this->apicobasal_distance_.coeffRef(lv_nid) = psi_am.coeff(lv_nid);
                this->apicobasal_direction_.row(lv_nid) = psi_am_grad.row(lv_nid);
            }
        }
    } // End of Compute apicobasal distance and direction at left ventricle.


    // Compute apicobasal distance and direction at right ventricle.
    if (!this->rv_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.RvApexTag(), ventri_tags.TricuspidValveTag()};
        std::vector<double> dirichlet_values{1., 0.};

        // Solve apical-tricuspid component.
        Eigen::VectorXd psi_at;
        Eigen::MatrixXd psi_at_grad;
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_at);
        this->ComputeGradient(mesh, psi_at, psi_at_grad);

        // Solve apical-pulmonary component if available.
        if (!ventri_tags.PulmonaryValveTag().empty()) {
            Eigen::VectorXd psi_ap;
            Eigen::MatrixXd psi_ap_grad;
            dirichlet_tags[1] = ventri_tags.PulmonaryValveTag();
            this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, psi_ap);
            this->ComputeGradient(mesh, psi_ap, psi_ap_grad);

            for (const auto &rv_nid : this->rv_node_ids_) {
                this->apicobasal_distance_.coeffRef(rv_nid) = this->intraventricular_function_.coeff(rv_nid)*psi_at.coeff(rv_nid) + (1.0-this->intraventricular_function_.coeff(rv_nid))*psi_ap.coeff(rv_nid);
                this->apicobasal_direction_.row(rv_nid) = this->intraventricular_function_.coeff(rv_nid)*psi_at_grad.row(rv_nid) + (1.0-this->intraventricular_function_.coeff(rv_nid))*psi_ap_grad.row(rv_nid);
            }
        } else {
            for (const auto &rv_nid : this->rv_node_ids_) {
                this->apicobasal_distance_.coeffRef(rv_nid) = psi_at.coeff(rv_nid);
                this->apicobasal_direction_.row(rv_nid) = psi_at_grad.row(rv_nid);
            }
        }

    } // End of Compute apicobasal distance and direction at right ventricle.

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriApicoBasal(const IMP::Voronoi<DIM> &voro, const CLOUDEA::Fpm<DIM> &fpm, const VentriTags &ventri_tags)
{
    // Check if at least one of the ventricles has been extracted.
    if (this->rv_node_ids_.empty() && this->lv_node_ids_.empty()) {
        std::string error_msg = "Could not compute apicobasal distance and direction. Extract ventricular node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Check if intraventricular interpolation function has been computed.
    if (this->intraventricular_function_.rows() == 0) {
        std::string error_msg = "Could not compute apicobasal distance and direction. Compute intraventricular interpolation function first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize apicobasal distance.
    this->apicobasal_distance_ = Eigen::VectorXd::Zero(voro.NodesNum());
    this->apicobasal_direction_ = Eigen::MatrixXd::Zero(voro.NodesNum(),DIM);

    // Compute apicobasal distance and direction at left ventricle.
    if (!this->lv_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.LvApexTag(), ventri_tags.MitralValveTag()};
        std::vector<double> dirichlet_values{1., 0.};

        // Solve apico-mitral component.
        Eigen::VectorXd psi_am;
        Eigen::MatrixXd psi_am_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_am);
        this->ComputeGradient(voro, fpm, psi_am, psi_am_grad);

        // Solve apico-aortic component if available.
        if (!ventri_tags.AorticValveTag().empty()) {
            Eigen::VectorXd psi_aa;
            Eigen::MatrixXd psi_aa_grad;
            dirichlet_tags[1] = ventri_tags.AorticValveTag();
            this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_aa);
            this->ComputeGradient(voro, fpm, psi_aa, psi_aa_grad);

            for (const auto &lv_nid : this->lv_node_ids_) {
                this->apicobasal_distance_.coeffRef(lv_nid) = this->intraventricular_function_.coeff(lv_nid)*psi_am.coeff(lv_nid) + (1.0-this->intraventricular_function_.coeff(lv_nid))*psi_aa.coeff(lv_nid);
                this->apicobasal_direction_.row(lv_nid) = this->intraventricular_function_.coeff(lv_nid)*psi_am_grad.row(lv_nid) + (1.0-this->intraventricular_function_.coeff(lv_nid))*psi_aa_grad.row(lv_nid);
            }
        } else {
            for (const auto &lv_nid : this->lv_node_ids_) {
                this->apicobasal_distance_.coeffRef(lv_nid) = psi_am.coeff(lv_nid);
                this->apicobasal_direction_.row(lv_nid) = psi_am_grad.row(lv_nid);
            }
        }
    } // End of Compute apicobasal distance and direction at left ventricle.


    // Compute apicobasal distance and direction at right ventricle.
    if (!this->rv_node_ids_.empty()) {

        // Set Dirichlet boundary conditions tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.RvApexTag(), ventri_tags.TricuspidValveTag()};
        std::vector<double> dirichlet_values{1., 0.};

        // Solve apical-tricuspid component.
        Eigen::VectorXd psi_at;
        Eigen::MatrixXd psi_at_grad;
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_at);
        this->ComputeGradient(voro, fpm, psi_at, psi_at_grad);

        // Solve apical-pulmonary component if available.
        if (!ventri_tags.PulmonaryValveTag().empty()) {
            Eigen::VectorXd psi_ap;
            Eigen::MatrixXd psi_ap_grad;
            dirichlet_tags[1] = ventri_tags.PulmonaryValveTag();
            this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, psi_ap);
            this->ComputeGradient(voro, fpm, psi_ap, psi_ap_grad);

            for (const auto &rv_nid : this->rv_node_ids_) {
                this->apicobasal_distance_.coeffRef(rv_nid) = this->intraventricular_function_.coeff(rv_nid)*psi_at.coeff(rv_nid) + (1.0-this->intraventricular_function_.coeff(rv_nid))*psi_ap.coeff(rv_nid);
                this->apicobasal_direction_.row(rv_nid) = this->intraventricular_function_.coeff(rv_nid)*psi_at_grad.row(rv_nid) + (1.0-this->intraventricular_function_.coeff(rv_nid))*psi_ap_grad.row(rv_nid);
            }
        } else {
            for (const auto &rv_nid : this->rv_node_ids_) {
                this->apicobasal_distance_.coeffRef(rv_nid) = psi_at.coeff(rv_nid);
                this->apicobasal_direction_.row(rv_nid) = psi_at_grad.row(rv_nid);
            }
        }

    } // End of Compute apicobasal distance and direction at right ventricle.

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriIntraventricular(const IMP::Mesh<DIM,CELL_NODES> &mesh, const VentriTags &ventri_tags)
{
    // Check if at least one of the ventricles has been extracted.
    if (this->rv_node_ids_.empty() && this->lv_node_ids_.empty()) {
        std::string error_msg = "Could not compute intraventicular interpolation function. Extract ventricular node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize intraventricular interpolation function.
    this->intraventricular_function_ = Eigen::VectorXd::Zero(mesh.NodesNum());

    // Solve left ventricle part.
    if (!this->lv_node_ids_.empty()) {
        // Check that at least the lv apex and mitral valve tags are available.
        if (ventri_tags.LvApexTag().empty() || ventri_tags.MitralValveTag().empty()) {
            std::string error_msg = "Could not compute intraventicular interpolation function at left ventricle."
                                    " At least the left ventricle apex and mitral valve ring tags are required.";
            throw std::runtime_error(Logger::Error(error_msg));
        }

        // Set dirichlet boundary tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.LvApexTag(), ventri_tags.MitralValveTag()};
        std::vector<double> dirichlet_values;
        if (ventri_tags.AorticValveTag().empty()) {
            dirichlet_values.emplace_back(1);  // lv apex value.
            dirichlet_values.emplace_back(0);  // mitral valve value.
        } else {
            dirichlet_tags.emplace_back(ventri_tags.AorticValveTag());
            dirichlet_values.emplace_back(1);  // lv apex value.
            dirichlet_values.emplace_back(1);  // mitral valve value.
            dirichlet_values.emplace_back(0);  // aortic valve value.
        }

        // Solution at left ventricle.
        Eigen::VectorXd lv_intraventricular = Eigen::VectorXd::Zero(mesh.NodesNum());
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, lv_intraventricular);
        for (const auto &lv_nid : this->lv_node_ids_) {
            this->intraventricular_function_.coeffRef(lv_nid) = lv_intraventricular.coeff(lv_nid);
        }

    } // End of Solve left ventricle part.

    // Solve right ventricle part.
    if (!this->rv_node_ids_.empty()) {
        // Check that at least the rv apex and tricuspid valve tags are available.
        if (ventri_tags.RvApexTag().empty() || ventri_tags.TricuspidValveTag().empty()) {
            std::string error_msg = "Could not compute intraventicular interpolation function at right ventricle." 
                                    " At least the right ventricle apex and tricuspid valve ring tags are required.";
            throw std::runtime_error(Logger::Error(error_msg));
        }

        // Set dirichlet boundary tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.RvApexTag(), ventri_tags.TricuspidValveTag()};
        std::vector<double> dirichlet_values;
        if (ventri_tags.PulmonaryValveTag().empty()) {
            dirichlet_values.emplace_back(1);  // rv apex value.
            dirichlet_values.emplace_back(0);  // tricuspid valve value.
        } else {
            dirichlet_tags.emplace_back(ventri_tags.PulmonaryValveTag());
            dirichlet_values.emplace_back(1);  // rv apex value.
            dirichlet_values.emplace_back(1);  // tricuspid valve value.
            dirichlet_values.emplace_back(0);  // pulmonary valve value.
        }

        // Solution at left ventricle.
        Eigen::VectorXd rv_intraventricular = Eigen::VectorXd::Zero(mesh.NodesNum());
        this->SolveLaplacian(mesh.NodeSets(), dirichlet_tags, dirichlet_values, rv_intraventricular);
        for (const auto &rv_nid : this->rv_node_ids_) {
            this->intraventricular_function_.coeffRef(rv_nid) = rv_intraventricular.coeff(rv_nid);
        }

    } // End of Solve right ventricle part.
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriIntraventricular(const IMP::Voronoi<DIM> &voro, const VentriTags &ventri_tags)
{
    // Check if at least one of the ventricles has been extracted.
    if (this->rv_node_ids_.empty() && this->lv_node_ids_.empty()) {
        std::string error_msg = "Could not compute intraventicular interpolation function. Extract ventricular node sets first.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Initialize intraventricular interpolation function.
    this->intraventricular_function_ = Eigen::VectorXd::Zero(voro.NodesNum());

    // Solve left ventricle part.
    if (!this->lv_node_ids_.empty()) {
        // Check that at least the lv apex and mitral valve tags are available.
        if (ventri_tags.LvApexTag().empty() || ventri_tags.MitralValveTag().empty()) {
            std::string error_msg = "Could not compute intraventicular interpolation function at left ventricle."
                                    " At least the left ventricle apex and mitral valve ring tags are required.";
            throw std::runtime_error(Logger::Error(error_msg));
        }

        // Set dirichlet boundary tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.LvApexTag(), ventri_tags.MitralValveTag()};
        std::vector<double> dirichlet_values;
        if (ventri_tags.AorticValveTag().empty()) {
            dirichlet_values.emplace_back(1);  // lv apex value.
            dirichlet_values.emplace_back(0);  // mitral valve value.
        } else {
            dirichlet_tags.emplace_back(ventri_tags.AorticValveTag());
            dirichlet_values.emplace_back(1);  // lv apex value.
            dirichlet_values.emplace_back(1);  // mitral valve value.
            dirichlet_values.emplace_back(0);  // aortic valve value.
        }

        // Solution at left ventricle.
        Eigen::VectorXd lv_intraventricular = Eigen::VectorXd::Zero(voro.NodesNum());
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, lv_intraventricular);
        for (const auto &lv_nid : this->lv_node_ids_) {
            this->intraventricular_function_.coeffRef(lv_nid) = lv_intraventricular.coeff(lv_nid);
        }

    } // End of Solve left ventricle part.

    // Solve right ventricle part.
    if (!this->rv_node_ids_.empty()) {
        // Check that at least the rv apex and tricuspid valve tags are available.
        if (ventri_tags.RvApexTag().empty() || ventri_tags.TricuspidValveTag().empty()) {
            std::string error_msg = "Could not compute intraventicular interpolation function at right ventricle."
                                    " At least the right ventricle apex and tricuspid valve ring tags are required.";
            throw std::runtime_error(Logger::Error(error_msg));
        }

        // Set dirichlet boundary tags and values.
        std::vector<std::string> dirichlet_tags{ventri_tags.RvApexTag(), ventri_tags.TricuspidValveTag()};
        std::vector<double> dirichlet_values;
        if (ventri_tags.PulmonaryValveTag().empty()) {
            dirichlet_values.emplace_back(1);  // rv apex value.
            dirichlet_values.emplace_back(0);  // tricuspid valve value.
        } else {
            dirichlet_tags.emplace_back(ventri_tags.PulmonaryValveTag());
            dirichlet_values.emplace_back(1);  // rv apex value.
            dirichlet_values.emplace_back(1);  // tricuspid valve value.
            dirichlet_values.emplace_back(0);  // pulmonary valve value.
        }

        // Solution at left ventricle.
        Eigen::VectorXd rv_intraventricular = Eigen::VectorXd::Zero(voro.NodesNum());
        this->SolveLaplacian(voro.NodeSets(), dirichlet_tags, dirichlet_values, rv_intraventricular);
        for (const auto &rv_nid : this->rv_node_ids_) {
            this->intraventricular_function_.coeffRef(rv_nid) = rv_intraventricular.coeff(rv_nid);
        }

    } // End of Solve right ventricle part.
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ExtractAtriumNodalPartitions(const IMP::Mesh<DIM,CELL_NODES> &mesh)
{
    if (this->transmural_distance_.rows() != mesh.NodesNum()) {
        std::string error_msg = "Could not extract atrium nodal partitions. Transmural distance is not available for all the nodes of the given mesh.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Set flag for left atrium nodes from transmural distance information.
    std::vector<bool> la_nodes_flag(mesh.NodesNum(), false);
    for (Eigen::Index i = 0; i != this->transmural_distance_.rows(); ++i) {
        if (this->transmural_distance_.coeff(i) < 0.) { la_nodes_flag[i] = true; }
    }

    // Find indices of cells with at least one node in left atrium.
    std::vector<int> la_cell_ids;
    la_cell_ids.reserve(mesh.CellsNum());
    for (const auto &cell : mesh.Cells()) {
        int cid = &cell - &mesh.Cells()[0];
        for (const auto &nid : cell.Connectivity()) {
            if (la_nodes_flag[nid]) { la_cell_ids.emplace_back(cid); break; }
        }
    }
    la_cell_ids.shrink_to_fit();

    // Update left atrium nodes flag.
    for (const auto &cid : la_cell_ids) {
        for (const auto &nid : mesh.Cells(cid).Connectivity()) {
            la_nodes_flag[nid] = true;
        }
    }

    // Set left and right atrium node indices.
    this->la_node_ids_.clear();  this->la_node_ids_.reserve(mesh.NodesNum());
    this->ra_node_ids_.clear();  this->ra_node_ids_.reserve(mesh.NodesNum());
    int it = 0;
    for (const auto &flag : la_nodes_flag) {
        if (flag) this->la_node_ids_.emplace_back(it);
        else this->ra_node_ids_.emplace_back(it);
        it++;
    }
    this->la_node_ids_.shrink_to_fit();
    this->ra_node_ids_.shrink_to_fit();
    la_cell_ids.clear();
    la_cell_ids.shrink_to_fit();
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ExtractVentricleNodalPartitions(const IMP::Mesh<DIM,CELL_NODES> &mesh)
{
    if (this->transmural_distance_.rows() != mesh.NodesNum()) {
        std::string error_msg = "Could not extract ventricle nodal partitions. Transmural distance is not available for all the nodes of the given mesh.";
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Set flag for left ventricle nodes from transmural distance information.
    std::vector<bool> lv_nodes_flag(mesh.NodesNum(), false);
    for (Eigen::Index i = 0; i != this->transmural_distance_.rows(); ++i) {
        if (this->transmural_distance_.coeff(i) < 0) { lv_nodes_flag[i] = true; }
    }

    // Find indices of cells with at least one node in left ventricle.
    std::vector<int> lv_cell_ids;
    lv_cell_ids.reserve(mesh.CellsNum());
    for (const auto &cell : mesh.Cells()) {
        int cid = &cell - &mesh.Cells()[0];
        for (const auto &nid : cell.Connectivity()) {
            if (lv_nodes_flag[nid]) { lv_cell_ids.emplace_back(cid); break; }
        }
    }
    lv_cell_ids.shrink_to_fit();

    // Update left ventricle nodes flag.
    for (const auto &cid : lv_cell_ids) {
        for (const auto &nid : mesh.Cells(cid).Connectivity()) {
            lv_nodes_flag[nid] = true;
        }
    }

    // Set left and right ventricle node indices.
    this->lv_node_ids_.clear();  this->lv_node_ids_.reserve(mesh.NodesNum());
    this->rv_node_ids_.clear();  this->rv_node_ids_.reserve(mesh.NodesNum());
    int it = 0;
    for (const auto &flag : lv_nodes_flag) {
        if (flag) this->lv_node_ids_.emplace_back(it);
        else this->rv_node_ids_.emplace_back(it);
        it++;
    }
    this->lv_node_ids_.shrink_to_fit();
    this->rv_node_ids_.shrink_to_fit();
    lv_cell_ids.clear();
    lv_cell_ids.shrink_to_fit();

    // Flag indices of nodes at the interface of the ventricles.
    std::vector<short> ventri_nodes_flag(mesh.NodesNum(), -1);
    for (const auto &id : this->rv_node_ids_) { ventri_nodes_flag[id] = 1; }
    short flag_sum = 0;
    std::vector<bool> interface_nodes_flag(mesh.NodesNum(), false);
    for (const auto &cell : mesh.Cells()) {
        for (const auto &nid : cell.Connectivity()) { flag_sum += ventri_nodes_flag[nid]; }
        if (std::abs(flag_sum) != CELL_NODES) {
            for (const auto &nid : cell.Connectivity()) { interface_nodes_flag[nid] = true; }
        }
        flag_sum = 0;
    }
    ventri_nodes_flag.clear();
    ventri_nodes_flag.shrink_to_fit();

    // Set the ventricle interface nodes container.
    this->ventricle_interface_node_ids_.clear();
    this->ventricle_interface_node_ids_.reserve(mesh.NodesNum());
    int i = 0;
    for (const auto &it : interface_nodes_flag) {
        if (it) { this->ventricle_interface_node_ids_.emplace_back(i); }
        i++;
    }
    interface_nodes_flag.clear();
    interface_nodes_flag.shrink_to_fit();
    this->ventricle_interface_node_ids_.shrink_to_fit();
}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeAtriFibers(const AtriFiberRules &rules)
{
    // Reset fiber matrices
    Eigen::Index rows_num = this->transmural_direction_.rows();
    Eigen::Index cols_num = this->transmural_direction_.cols();
    this->fiber_direction_ = Eigen::MatrixXd::Zero(rows_num, cols_num);
    this->sheet_direction_ = Eigen::MatrixXd::Zero(rows_num, cols_num);

    // Select atria bundles
    Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(rows_num, cols_num);
    if (!this->la_node_ids_.empty()) {
        // Get the left atrium bundles.
        for (const auto &la_nid : this->la_node_ids_) {
            if (this->valve_veins_distance_.coeff(la_nid) >= rules.MitralValveRule()) {
                psi.row(la_nid) = this->valve_veins_direction_.row(la_nid);
            } else {
                if (this->inter_veins_distance_.coeff(la_nid) >= rules.LeftPulmonVeinRule() ||
                        this->inter_veins_distance_.coeff(la_nid) <= rules.RightPulmonVeinRule()) {
                    psi.row(la_nid) = this->inter_veins_direction_.row(la_nid);
                } else {
                    psi.row(la_nid) = this->appendage_veins_direction_.row(la_nid);
                }
            }
        }
    }
    if (!this->ra_node_ids_.empty()) {
        // Get the right atrium bundles.
        for (const auto &ra_nid : this->ra_node_ids_) {
            if (this->valve_veins_distance_.coeff(ra_nid) >= rules.TricuspidValveRule() ) {
                psi.row(ra_nid) = this->valve_veins_direction_.row(ra_nid);
            } else {
                if (this->valve_veins_distance_.coeff(ra_nid) < rules.RightLateralWallRule()) {

                    if (this->atri_tricuspid_distance_.coeff(ra_nid) >= rules.CristaTerminalMinusRule() &&
                            this->atri_tricuspid_distance_.coeff(ra_nid) <= rules.CristaTerminalPlusRule()) {
                        psi.row(ra_nid) = this->atri_tricuspid_direction_.row(ra_nid);
                    } else if (this->atri_tricuspid_distance_.coeff(ra_nid) < rules.CristaTerminalMinusRule()) {

                        if (this->inter_veins_distance_.coeff(ra_nid) >= rules.InfCavalVeinRule() ||
                                this->inter_veins_distance_.coeff(ra_nid) <= rules.SupCavalVeinRule()) {
                            psi.row(ra_nid) = this->inter_veins_direction_.row(ra_nid);
                        } else {
                            psi.row(ra_nid) = this->appendage_veins_direction_.row(ra_nid);
                        }

                    } else {
                        if (this->inter_veins_distance_.coeff(ra_nid) >= rules.InfCavalVeinRule() ||
                                this->inter_veins_distance_.coeff(ra_nid) <= rules.SupCavalVeinRule()) {
                            psi.row(ra_nid) = this->inter_veins_direction_.row(ra_nid);
                        } else {
                            if (this->atri_tricuspid_distance_.coeff(ra_nid) <= rules.InterCavalBundleRule()) {
                                psi.row(ra_nid) = this->inter_veins_direction_.row(ra_nid);
                            } else if (this->atri_tricuspid_distance_.coeff(ra_nid) >= rules.RightSeptumWallRule()) {
                                psi.row(ra_nid) = this->valve_veins_direction_.row(ra_nid);
                            } else {
                                psi.row(ra_nid) = this->atri_tricuspid_direction_.row(ra_nid);
                            }
                        }
                    }

                } else {
                    if (this->inter_veins_distance_.coeff(ra_nid) >= rules.InfCavalVeinRule() ||
                            this->inter_veins_distance_.coeff(ra_nid) <= rules.SupCavalVeinRule()) {
                        psi.row(ra_nid) = this->inter_veins_direction_.row(ra_nid);
                    } else {
                        if (this->atri_tricuspid_distance_.coeff(ra_nid) >= 0.) {
                            psi.row(ra_nid) = this->valve_veins_direction_.row(ra_nid);
                        } else {
                            psi.row(ra_nid) = this->appendage_veins_direction_.row(ra_nid);
                        }
                    }
                }
            }
        } // End of Get the right atrium bundles.
    }

    // Get transmural direction vector field.
    Eigen::MatrixXd phi = this->transmural_direction_;
    if (!this->la_node_ids_.empty()) {
        // Invert transmural direction at left atrium.
        for (const auto &la_nid : this->la_node_ids_) {
            phi.row(la_nid) = -1.*phi.row(la_nid);
        }
    }

    // Compute fibers.
    Eigen::Vector3d el(0.,0.,0.), et(0.,0.,0.), ec(0.,0.,0.), ephi(0.,0.,0.);
    for (Eigen::Index i = 0; i != rows_num; ++i) {

        // Compute local azimuthal directional axis.
        el = psi.row(i);
        if (el.norm() > 2.*std::numeric_limits<double>::epsilon()) el /= el.norm();

        // Compute local radial directional axis.
        ephi = phi.row(i);
        et = ephi - (el.dot(ephi))*el;
        if (et.norm() > 2.*std::numeric_limits<double>::epsilon()) et /= et.norm();

        // Compute circumferential directional axis.
        ec = et.cross(el);

        this->sheet_direction_.row(i) = et;
        this->fiber_direction_.row(i) = ec;
    }

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::ComputeVentriFibers(const VentriFiberRules &rules)
{
    // Reset fiber matrices
    Eigen::Index rows_num = this->apicobasal_direction_.rows();
    Eigen::Index cols_num = this->apicobasal_direction_.cols();
    this->fiber_direction_ = Eigen::MatrixXd::Zero(rows_num, cols_num);
    this->sheet_direction_ = Eigen::MatrixXd::Zero(rows_num, cols_num);

    // Compute normalized transmural distance and set left ventricle flag.
    Eigen::VectorXd ntd = this->transmural_distance_;
    Eigen::SparseVector<short> lv_flag(rows_num);
    if (!this->lv_node_ids_.empty()) {
        // Normalize values such as epi->endo is 0->1 instead of 0->(-2).
        for (const auto &lv_nid : this->lv_node_ids_) {
            ntd.coeffRef(lv_nid) /= -2.;
            lv_flag.coeffRef(lv_nid) = 1;
        }
    }

    // Set septum flag.
    Eigen::SparseVector<short> septum_flag(rows_num);
    if (!this->septum_node_ids_.empty()) {
        for (const auto &sp_nid : this->septum_node_ids_) {
            septum_flag.coeffRef(sp_nid) = 1;
        }
    }

    // Compute fibers.
    const double pi = 3.14159265358979323846;
    Eigen::Vector3d el(0.,0.,0.), et(0.,0.,0.), ec(0.,0.,0.);
    Eigen::Vector3d phi(0.,0.,0.), psi(0.,0.,0.);
    Eigen::Affine3d rot;
    double a = 0., a_endo = 0., a_epi = 0.;
    double b = 0., b_endo = 0., b_epi = 0.;
    for (Eigen::Index i = 0; i != rows_num; ++i) {
        phi = this->transmural_direction_.row(i);
        psi = this->apicobasal_direction_.row(i);

        // Compute local azimuthal directional axis.
        el = psi;
        if (el.norm() > 2.*std::numeric_limits<double>::epsilon()) el /= el.norm();

        // Compute local radial directional axis.
        et = phi - (el.dot(phi))*el;
        if (et.norm() > 2.*std::numeric_limits<double>::epsilon()) et /= et.norm();

        // Compute circumferential directional axis.
        ec = et.cross(el);

        // Compute endocardial and epicardial rotation angles for given fiber.
        if (lv_flag.coeff(i) == 1) {
            // Angles for left ventricle.
            a_endo = this->intraventricular_function_.coeff(i) * (rules.AlphaLvEndo() - rules.AlphaOtLvEndo()) + rules.AlphaOtLvEndo();
            a_epi = this->intraventricular_function_.coeff(i) * (rules.AlphaLvEpi() - rules.AlphaOtLvEpi()) + rules.AlphaOtLvEpi();

            b_endo = this->intraventricular_function_.coeff(i) * (rules.BetaLvEndo() - rules.BetaOtLvEndo()) + rules.BetaOtLvEndo();
            b_epi = this->intraventricular_function_.coeff(i) * (rules.BetaLvEpi() - rules.BetaOtLvEpi()) + rules.BetaOtLvEpi();
        } else {
            // Angles for right ventricle.
            a_endo = this->intraventricular_function_.coeff(i) * (rules.AlphaRvEndo() - rules.AlphaOtRvEndo()) + rules.AlphaOtRvEndo();
            a_epi = this->intraventricular_function_.coeff(i) * (rules.AlphaRvEpi() - rules.AlphaOtRvEpi()) + rules.AlphaOtRvEpi();

            b_endo = this->intraventricular_function_.coeff(i) * (rules.BetaRvEndo() - rules.BetaOtRvEndo()) + rules.BetaOtRvEndo();
            b_epi = this->intraventricular_function_.coeff(i) * (rules.BetaRvEpi() - rules.BetaOtRvEpi()) + rules.BetaOtRvEpi();
        }

        // Compute final rotation angles of the fiber.
        a = a_epi*(1. - ntd.coeff(i)) + a_endo*ntd.coeff(i);
        if (septum_flag.coeff(i) == 1) {
            a = a*(1. - this->septal_distance_.coeff(i)) + rules.AlphaSeptum()*this->septal_distance_.coeff(i);
        }
        b = b_epi*(1. - ntd.coeff(i)) + b_endo*ntd.coeff(i);

        // Rotation around the tranversal axis.
        rot = Eigen::AngleAxisd(-(a*pi)/180., et);
        ec = rot.linear()*ec;

        // Rotation around the longitudinal axis.
        rot = Eigen::AngleAxisd((b*pi)/180., el);
        ec = rot.linear()*ec;
        et = rot.linear()*et;

        this->sheet_direction_.row(i) = et;
        this->fiber_direction_.row(i) = ec;

    }

}


#include <fstream>
template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::TestSaveAtriFibers(const IMP::Mesh<DIM,CELL_NODES> &mesh, const std::string &outname)
{
    std::ofstream out(outname);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // Write nodes
    out << "POINTS " << mesh.NodesNum() << " float\n";
    for (const auto &node : mesh.Nodes()) {
        out << node[0] << " " << node[1] << " " << node[2] << "\n";
    }

    // Write elements
    out << "\nCELLS " << mesh.CellsNum() << " " << mesh.CellsNum()*5 << "\n";
    for (const auto &cell : mesh.Cells()) {
        out << "4 " << cell.N(0) << " " << cell.N(1) << " " << cell.N(2) << " " << cell.N(3) << "\n";
    }

    out << "\nCELL_TYPES " << mesh.CellsNum() << "\n";
    int el_id = 0;
    while (el_id != mesh.CellsNum()) {
        out << "10\n";
        el_id++;
    }
    out << "\n";

    out << "POINT_DATA " << mesh.NodesNum() << "\n";
    out << "FIELD FieldData 14\n";
    out << "TransmuralDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->transmural_distance_.rows(); ++i) {
        out << this->transmural_distance_.coeff(i) << "\n";
    }
    out << "\nTransmuralDirection 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->transmural_direction_.rows(); ++i) {
        out << this->transmural_direction_.row(i) << "\n";
    }
    out << "\nAppendageVeinsDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->appendage_veins_distance_.rows(); ++i) {
        out << this->appendage_veins_distance_.coeff(i) << "\n";
    }
    out << "\nAppendageVeinsDirection 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->appendage_veins_direction_.rows(); ++i) {
        out << this->appendage_veins_direction_.row(i) << "\n";
    }
    out << "\nInterVeinsDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->inter_veins_distance_.rows(); ++i) {
        out << this->inter_veins_distance_.coeff(i) << "\n";
    }
    out << "\nInterVeinsDirection 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->inter_veins_direction_.rows(); ++i) {
        out << this->inter_veins_direction_.row(i) << "\n";
    }
    out << "\nValveVeinsDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->valve_veins_distance_.rows(); ++i) {
        out << this->valve_veins_distance_.coeff(i) << "\n";
    }
    out << "\nValveVeinsDirection 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->valve_veins_direction_.rows(); ++i) {
        out << this->valve_veins_direction_.row(i) << "\n";
    }
    out << "\nAtriTricuspidDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->atri_tricuspid_distance_.rows(); ++i) {
        out << this->atri_tricuspid_distance_.coeff(i) << "\n";
    }
    out << "\nAtriTricuspidDirection 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->atri_tricuspid_direction_.rows(); ++i) {
        out << this->atri_tricuspid_direction_.row(i) << "\n";
    }

    out << "\nFibers 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->fiber_direction_.rows(); ++i) {
        out << this->fiber_direction_.row(i) << "\n";
    }

    out << "\nSheet 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->sheet_direction_.rows(); ++i) {
        out << this->sheet_direction_.row(i) << "\n";
    }

    out << "\nLAnodes 1 " << mesh.NodesNum() << " DOUBLE\n";
    std::vector<int> la_flag(mesh.NodesNum(), 0);
    for (const auto &id : this->la_node_ids_) { la_flag[id] = 1; }
    for (const auto &flag : la_flag) {
        out << flag << "\n";
    }

    out << "\nRAnodes 1 " << mesh.NodesNum() << " DOUBLE\n";
    std::vector<int> ra_flag(mesh.NodesNum(), 0);
    for (const auto &id : this->ra_node_ids_) { ra_flag[id] = 1; }
    for (const auto &flag : ra_flag) {
        out << flag << "\n";
    }


    out.close();

}


template<short DIM, short CELL_NODES>
void Ldrbm<DIM, CELL_NODES>::TestSaveVentriFibers(const IMP::Mesh<DIM,CELL_NODES> &mesh, const std::string &outname)
{
    std::ofstream out(outname);

    out << "# vtk DataFile Version 3.0\n";
    out << "vtk output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // Write nodes
    out << "POINTS " << mesh.NodesNum() << " float\n";
    for (const auto &node : mesh.Nodes()) {
        out << node[0] << " " << node[1] << " " << node[2] << "\n";
    }

    // Write elements
    out << "\nCELLS " << mesh.CellsNum() << " " << mesh.CellsNum()*5 << "\n";
    for (const auto &cell : mesh.Cells()) {
        out << "4 " << cell.N(0) << " " << cell.N(1) << " " << cell.N(2) << " " << cell.N(3) << "\n";
    }

    out << "\nCELL_TYPES " << mesh.CellsNum() << "\n";
    int el_id = 0;
    while (el_id != mesh.CellsNum()) {
        out << "10\n";
        el_id++;
    }
    out << "\n";

    out << "POINT_DATA " << mesh.NodesNum() << "\n";
    out << "FIELD FieldData 12\n";
    out << "TransmuralDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->transmural_distance_.rows(); ++i) {
        out << this->transmural_distance_.coeff(i) << "\n";
    }
    out << "\nTransmuralDirection 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->transmural_direction_.rows(); ++i) {
        out << this->transmural_direction_.row(i) << "\n";
    }
    out << "\nApicoBasalDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->apicobasal_distance_.rows(); ++i) {
        out << this->apicobasal_distance_.coeff(i) << "\n";
    }
    out << "\nApicoBasalDirection 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->apicobasal_direction_.rows(); ++i) {
        out << this->apicobasal_direction_.row(i) << "\n";
    }
    out << "SeptalDistance 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->septal_distance_.rows(); ++i) {
        out << this->septal_distance_.coeff(i) << "\n";
    }

    out << "\nFibers 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->fiber_direction_.rows(); ++i) {
        out << this->fiber_direction_.row(i) << "\n";
    }

    out << "\nSheet 3 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->sheet_direction_.rows(); ++i) {
        out << this->sheet_direction_.row(i) << "\n";
    }

    out << "\nIntraventricularFunction 1 " << mesh.NodesNum() << " DOUBLE\n";
    for (Eigen::Index i = 0; i != this->intraventricular_function_.rows(); ++i) {
        out << this->intraventricular_function_.coeff(i) << "\n";
    }

    out << "\nLVnodes 1 " << mesh.NodesNum() << " DOUBLE\n";
    std::vector<int> lv_flag(mesh.NodesNum(), 0);
    for (const auto &id : this->lv_node_ids_) { lv_flag[id] = 1; }
    for (const auto &flag : lv_flag) {
        out << flag << "\n";
    }

    out << "\nRVnodes 1 " << mesh.NodesNum() << " DOUBLE\n";
    std::vector<int> rv_flag(mesh.NodesNum(), 0);
    for (const auto &id : this->rv_node_ids_) { rv_flag[id] = 1; }
    for (const auto &flag : rv_flag) {
        out << flag << "\n";
    }

    out << "\nInterface 1 " << mesh.NodesNum() << " DOUBLE\n";
    std::vector<int> interface_flag(mesh.NodesNum(), 0);
    for (const auto &id : this->ventricle_interface_node_ids_) { interface_flag[id] = 1; }
    for (const auto &flag : interface_flag) {
        out << flag << "\n";
    }

    out << "\nSeptum 1 " << mesh.NodesNum() << " DOUBLE\n";
    std::vector<int> septum_flag(mesh.NodesNum(), 0);
    for (const auto &id : this->septum_node_ids_) { septum_flag[id] = 1; }
    for (const auto &flag : septum_flag) {
        out << flag << "\n";
    }


    out.close();

}


} // End of namespace ELECTRA

#endif //ELECTRA_FIBERS_LDRBM_TPP_