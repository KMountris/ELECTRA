/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_TPP_

#include "config_output.hpp"

namespace APP_ELECTRA
{

template<short DIM, short CELL_NODES>
ConfigOutput<DIM, CELL_NODES>::ConfigOutput()
{}


template<short DIM, short CELL_NODES>
ConfigOutput<DIM, CELL_NODES>::~ConfigOutput()
{}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputGeneration(const Parser &parser, const std::vector<IMP::Vec<DIM,double>> &nodes,
        const std::vector<IMP::Cell<DIM,CELL_NODES>> &cells, const std::shared_ptr<ReactionDiffusion<DIM, CELL_NODES>> &react_diff,
        const ELECTRA::PostProcess &post_process, std::ostream &stream) const
{
    // Output if paraview format was requested.
    if (parser.HasAttribute("output.paraview")) {
        this->OutputToParaview(parser, nodes, cells, react_diff, stream);
    }

    // Output if ensight format was requested.
    if (parser.HasAttribute("output.ensight")) {
        this->OutputToEnsight(parser, nodes, cells, react_diff, post_process, stream);
    }

    // Output if ascii format was requested.
    if (parser.HasAttribute("output.ascii")) {
        this->OutputToAscii(parser, post_process, stream);
    }


    // Output if binary format was requested.
    if (parser.HasAttribute("output.binary")) {
        this->OutputToBinary(parser, react_diff, stream);
    }
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputFibers(const Parser &parser, IMP::Mesh<DIM,CELL_NODES> mesh, IMP::Voronoi<DIM> voro, const Ldrbm<DIM,CELL_NODES> &ldrbm, std::ostream &stream) const
{
    // Output partitioned mesh
    if (parser.HasAttribute("output.fibers.partitioned geometry")) {
        std::string method = parser.GetValue<std::string>("numerical approximation.method");
        std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);

        std::string filename = parser.GetValue<std::string>("output.fibers.partitioned geometry");

        // Add ventricular node sets in the mesh.
        if (ldrbm.SeptumNodeIds().size() != 0) {
            // Left ventricle nodes.
            IMP::NodeSet lv_nset; lv_nset.Set("left_ventricle", ldrbm.LeftVentriNodeIds());

            // Right ventricle nodes.
            IMP::NodeSet rv_nset; rv_nset.Set("right_ventricle", ldrbm.RightVentriNodeIds());

            // Interface ventricle nodes.
            IMP::NodeSet intv_nset; intv_nset.Set("ventricle_interface", ldrbm.VentriInterfaceNodeIds());

            // Septum ventricle nodes.
            IMP::NodeSet septum_nset; septum_nset.Set("septum", ldrbm.VentriInterfaceNodeIds());

            if (method == "fem") {
                mesh.AddNodeSet(lv_nset);
                mesh.AddNodeSet(rv_nset);
                mesh.AddNodeSet(intv_nset);
                mesh.AddNodeSet(septum_nset);
            } else {
                voro.AddNodeSet(lv_nset);
                voro.AddNodeSet(rv_nset);
                voro.AddNodeSet(intv_nset);
                voro.AddNodeSet(septum_nset);
            }

        } else {
            // Add atrial node sets in the mesh.
            // Left atrium nodes.
            IMP::NodeSet la_nset; la_nset.Set("left_atrium", ldrbm.LeftAtriNodeIds());

            // Right atrium nodes.
            IMP::NodeSet ra_nset; ra_nset.Set("right_atrium", ldrbm.RightAtriNodeIds());

            if (method == "fem") {
                mesh.AddNodeSet(la_nset);
                mesh.AddNodeSet(ra_nset);
            } else {
                voro.AddNodeSet(la_nset);
                voro.AddNodeSet(ra_nset);
            }
        }

        // Save the mesh.
        // if (method == "fem") {
            mesh.SaveTo(filename);
            stream << Logger::Message("Saved partitioned mesh: " + filename + "\n");
        // } else {
        //     voro.SaveTo(filename);
        //     stream << Logger::Message("Saved partitioned voronoi tesselation: " + filename + "\n");
        // }
    }

    // Output longitudinal fibers.
    if (parser.HasAttribute("output.fibers.longitudinal fibers")) {
        std::string fibs_file = parser.GetValue<std::string>("output.fibers.longitudinal fibers");
        this->SaveFibers(ldrbm.LongFiberDirection(), fibs_file);
        stream << Logger::Message("Saved longitudinal fibers: " + fibs_file + "\n");
    }

    // Output sheet fibers.
    if (parser.HasAttribute("output.fibers.sheet fibers")) {
        std::string fibs_file = parser.GetValue<std::string>("output.fibers.sheet fibers");
        this->SaveFibers(ldrbm.SheetFiberDirection(), fibs_file);
        stream << Logger::Message("Saved sheet fibers: " + fibs_file + "\n");
    }

    // Output transversal fibers.
    if (parser.HasAttribute("output.fibers.transversal fibers")) {
        std::string fibs_file = parser.GetValue<std::string>("output.fibers.transversal fibers");

        // Compute transversal fibers as the outer product of longitudinal and sheet fibers.
        Eigen::MatrixXd trans_fibers = Eigen::MatrixXd::Zero(ldrbm.SheetFiberDirection().rows(), ldrbm.SheetFiberDirection().cols());
        Eigen::Vector3d v_l, v_s;
        for (Eigen::Index i = 0; i != ldrbm.SheetFiberDirection().rows(); ++i) {
            v_l = ldrbm.LongFiberDirection().row(i);
            v_s = ldrbm.SheetFiberDirection().row(i);
            trans_fibers.row(i) = v_l.cross(v_s);
        }

        this->SaveFibers(trans_fibers, fibs_file);
        stream << Logger::Message("Saved transversal fibers: " + fibs_file + "\n");
    }

    // Output transmural distance.
    if (parser.HasAttribute("output.fibers.transmural distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.transmural distance");
        this->SaveDistanceField(ldrbm.TransmuralDistance(), filename);
        stream << Logger::Message("Saved transmural distance: " + filename + "\n");
    }

    // Output appendage-veins distance.
    if (parser.HasAttribute("output.fibers.appendage veins distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.appendage veins distance");
        this->SaveDistanceField(ldrbm.AppendageVeinsDistance(), filename);
        stream << Logger::Message("Saved appendage veins distance: " + filename + "\n");
    }

    // Output interveins distance.
    if (parser.HasAttribute("output.fibers.interveins distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.interveins distance");
        this->SaveDistanceField(ldrbm.InterVeinsDistance(), filename);
        stream << Logger::Message("Saved interveins distance: " + filename + "\n");
    }

    // Output valve-veins distance.
    if (parser.HasAttribute("output.fibers.valve veins distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.valve veins distance");
        this->SaveDistanceField(ldrbm.ValveVeinsDistance(), filename);
        stream << Logger::Message("Saved valve veins distance: " + filename + "\n");
    }

    // Output atrium-tricuspid distance.
    if (parser.HasAttribute("output.fibers.atrium tricuspid distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.atrium tricuspid distance");
        this->SaveDistanceField(ldrbm.AtriTricuspidDistance(), filename);
        stream << Logger::Message("Saved atrium tricuspid distance: " + filename + "\n");
    }

    // Output apicobasal distance.
    if (parser.HasAttribute("output.fibers.apicobasal distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.apicobasal distance");
        this->SaveDistanceField(ldrbm.ApicobasalDistance(), filename);
        stream << Logger::Message("Saved apicobasal distance: " + filename + "\n");
    }

    // Output septal distance.
    if (parser.HasAttribute("output.fibers.septal distance")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.septal distance");
        this->SaveDistanceField(ldrbm.SeptalDistance(), filename);
        stream << Logger::Message("Saved septal distance: " + filename + "\n");
    }

    // Output intraventricular function.
    if (parser.HasAttribute("output.fibers.intraventricular function")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.intraventricular function");
        this->SaveDistanceField(ldrbm.IntraventricularFunction(), filename);
        stream << Logger::Message("Saved intraventricular function: " + filename + "\n");
    }

    // Output transmural direction.
    if (parser.HasAttribute("output.fibers.transmural direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.transmural direction");
        this->SaveDirectionField(ldrbm.TransmuralDirection(), filename);
        stream << Logger::Message("Saved transmural direction: " + filename + "\n");
    }

    // Output appendage-veins direction.
    if (parser.HasAttribute("output.fibers.appendage veins direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.appendage veins direction");
        this->SaveDirectionField(ldrbm.AppendageVeinsDirection(), filename);
        stream << Logger::Message("Saved appendage veins direction: " + filename + "\n");
    }

    // Output interveins direction.
    if (parser.HasAttribute("output.fibers.interveins direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.interveins direction");
        this->SaveDirectionField(ldrbm.InterVeinsDirection(), filename);
        stream << Logger::Message("Saved interveins direction: " + filename + "\n");
    }

    // Output valve-veins direction.
    if (parser.HasAttribute("output.fibers.valve veins direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.valve veins direction");
        this->SaveDirectionField(ldrbm.ValveVeinsDirection(), filename);
        stream << Logger::Message("Saved valve veins direction: " + filename + "\n");
    }

    // Output atrium-tricuspid direction.
    if (parser.HasAttribute("output.fibers.atrium tricuspid direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.atrium tricuspid direction");
        this->SaveDirectionField(ldrbm.AtriTricuspidDirection(), filename);
        stream << Logger::Message("Saved atrium tricuspid direction: " + filename + "\n");
    }

    // Output apicobasal direction.
    if (parser.HasAttribute("output.fibers.apicobasal direction")) {
        std::string filename = parser.GetValue<std::string>("output.fibers.apicobasal direction");
        this->SaveDirectionField(ldrbm.ApicobasalDirection(), filename);
        stream << Logger::Message("Saved apicobasal direction: " + filename + "\n");
    }

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputToParaview(const Parser &parser, const std::vector<IMP::Vec<int(DIM), double>> &nodes,
        const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM, CELL_NODES>> &react_diff, std::ostream &stream) const
{
    std::string states_file = parser.GetValue<std::string>("output.paraview.tissue.states");
    std::filesystem::path states_path(states_file);
    if (states_path.has_extension()) {
        std::filesystem::path e = "";
        states_path.replace_extension(e);
        states_file = states_path.string();
    }

    // Initialize paraview exporte
    ELECTRA::ParaviewExporter<DIM, CELL_NODES> para_export;
    para_export.CreateVtu(nodes, cells);

    // Export react_diff solution states.
    for (std::size_t i = 0; i != react_diff->Vout().size(); ++i) {
        para_export.AddScalarField(react_diff->Vout(i), "Potential");
        para_export.Export(states_file+std::to_string(i)+".vtu");
        para_export.ClearPointDataFields();
    }

    std::string anim_file = parser.GetValue<std::string>("output.paraview.tissue.animation");
    std::string anim_ext = std::filesystem::path(anim_file).extension();
    if (anim_ext != ".pvd")  anim_file += ".pvd";

    // Get states file directory path.
    std::string states_dir = states_path.parent_path().string();

    // Look up for animation directory path.
    std::filesystem::path anim_path(anim_file);
    std::string anim_dir = anim_path.parent_path().string();
    if (anim_dir.empty())  anim_file = states_dir + anim_file;

    para_export.CreatePvdAnimation(states_path, react_diff->Vout().size(), anim_file, react_diff->OutputSteps()*react_diff->Dt());

    stream << ELECTRA::Logger::Warning("Only potential states are stored. If you want to store more variables export in Ensight format\n");

    stream << ELECTRA::Logger::Message("Saved Paraview solutions: "+states_file+".vtu\n");
    stream << ELECTRA::Logger::Message("Saved Paraview animation: "+anim_file+"\n");

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputToEnsight(const Parser &parser, const std::vector<IMP::Vec<int(DIM), double>> &nodes,
        const std::vector<IMP::Cell<DIM, CELL_NODES>> &cells, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM, CELL_NODES>> &react_diff, const ELECTRA::PostProcess &post_process, std::ostream &stream) const
{
    // Initialize ensight exporter.
    ELECTRA::EnsightExporter<DIM, CELL_NODES> ens_export;

    // Save geometry to file.
    std::string geo_file = parser.GetValue<std::string>("output.ensight.tissue.geometry");
    std::string geo_ext = std::filesystem::path(geo_file).extension();
    if (geo_ext != ".geo")  geo_file += ".geo";
    ens_export.SaveGeo(nodes, cells, geo_file);

    // Get states file path and strip down the extension.
    std::string states_file = parser.GetValue<std::string>("output.ensight.tissue.states");
    std::filesystem::path states_path(states_file);
    if (states_path.has_extension()) {
        std::filesystem::path e = "";
        states_path.replace_extension(e);
        states_file = states_path.string();
    }

    // Number of solution states.
    std::size_t states_num = react_diff->Vout().size();

    // Create wildcard string.
    std::size_t wildcard_size = std::to_string(states_num).size();
    std::string wildcard(wildcard_size, '0');

    // Export react_diff solution states.
    std::string state_id = "";
    std::size_t replace_start = 0;
    for (std::size_t i = 0; i != states_num; ++i) {
        state_id = std::to_string(i);
        replace_start = wildcard_size - state_id.size();
        wildcard.replace(replace_start,std::string::npos,state_id);
        ens_export.SaveScalarField(react_diff->Vout(i).head(nodes.size()), states_file+wildcard+".ens");
    }

    // Additional scalar fields and their names.
    std::vector<std::string> scalar_field_files, scalar_field_names;

    // Local activation time scalar field.
    if (parser.HasAttribute("output.ensight.tissue.lat")) {
        std::string lat_file = parser.GetValue<std::string>("output.ensight.tissue.lat");
        std::filesystem::path lat_path(lat_file);
        if (lat_path.has_extension()) {
            std::filesystem::path e = "";
            lat_path.replace_extension(e);
            lat_file = lat_path.string();
        }
        ens_export.SaveScalarField(post_process.LocalActivationTimes().head(nodes.size()), lat_file+".ens");
        scalar_field_files.emplace_back(lat_file+".ens");
        scalar_field_names.emplace_back("LAT");

        stream << ELECTRA::Logger::Message("Saved local activation time: " + lat_file + ".ens\n");
    }

    // Action potential duration scalar field.
    if (parser.HasAttribute("output.ensight.tissue.apd")) {
        std::string apd_file = parser.GetValue<std::string>("output.ensight.tissue.apd");
        std::filesystem::path apd_path(apd_file);
        if (apd_path.has_extension()) {
            std::filesystem::path e = "";
            apd_path.replace_extension(e);
            apd_file = apd_path.string();
        }

        for (const auto &apd : post_process.APDs()) {
            ens_export.SaveScalarField(apd.second.head(nodes.size()), apd_file+std::to_string(apd.first)+".ens");
            scalar_field_files.emplace_back(apd_file+std::to_string(apd.first)+".ens");
            scalar_field_names.emplace_back("APD"+std::to_string(apd.first));
        }

        stream << ELECTRA::Logger::Message("Saved action potential duration: " + apd_file + ".ens\n");
    }

    std::string anim_file = parser.GetValue<std::string>("output.ensight.tissue.animation");
    std::string anim_ext = std::filesystem::path(anim_file).extension();
    if (anim_ext != ".case") { anim_file += ".case"; }

    // Get states file directory path.
    std::string states_dir = states_path.parent_path().string();

    // Look up for animation file path.
    std::filesystem::path anim_path(anim_file);
    if (!anim_path.has_root_path())
        anim_file = states_dir + anim_file;

    ens_export.SaveAnimation(anim_file, geo_file, states_file, scalar_field_files, scalar_field_names, states_num, react_diff->OutputSteps()*react_diff->Dt());

    stream << ELECTRA::Logger::Message("Saved Ensight geometry: " + geo_file + "\n");
    stream << ELECTRA::Logger::Message("Saved Ensight solutions: " + states_file + ".ens\n");
    stream << ELECTRA::Logger::Message("Saved Ensight animation: " + anim_file + "\n");

    if (react_diff->ConductSystem().NodesNum() > 0) {
        // Initialize ensight exporter for conduction systems.
        ELECTRA::EnsightExporter<DIM, 2> cs_ens_export;

        // Save geometry to file.
        if (parser.HasAttribute("output.ensight.conduction system.geometry")) {
            geo_file = parser.GetValue<std::string>("output.ensight.conduction system.geometry");
            std::string geo_ext = std::filesystem::path(geo_file).extension();
            if (geo_ext != ".geo")  geo_file += ".geo";
            cs_ens_export.SaveGeo(react_diff->ConductSystem().Nodes(), react_diff->ConductSystem().Segments(), geo_file);
        }

        if (parser.HasAttribute("output.ensight.conduction system.states")) {
            states_file = parser.GetValue<std::string>("output.ensight.conduction system.states");
            std::filesystem::path states_path(states_file);
            if (states_path.has_extension()) {
                std::filesystem::path e = "";
                states_path.replace_extension(e);
                states_file = states_path.string();
            }

            // Empty wildcard string.
            wildcard.clear(); wildcard.resize(wildcard_size, '0');

            // Export react_diff solution states.
            state_id = "";
            replace_start = 0;
            for (std::size_t i = 0; i != states_num; ++i) {
                state_id = std::to_string(i);
                replace_start = wildcard_size - state_id.size();
                wildcard.replace(replace_start,std::string::npos,state_id);
                cs_ens_export.SaveScalarField(react_diff->Vout(i).tail(react_diff->ConductSystem().NodesNum()), states_file+wildcard+".ens");
            }
        }

        // Reset additional scalar fields and their names.
        scalar_field_files.clear(); scalar_field_names.clear();

        if (parser.HasAttribute("output.ensight.conduction system.lat")) {
            std::string lat_file = parser.GetValue<std::string>("output.ensight.conduction system.lat");
            std::filesystem::path lat_path(lat_file);
            if (lat_path.has_extension()) {
                std::filesystem::path e = "";
                lat_path.replace_extension(e);
                lat_file = lat_path.string();
            }
            cs_ens_export.SaveScalarField(post_process.LocalActivationTimes().tail(react_diff->ConductSystem().NodesNum()), lat_file+".ens");
            scalar_field_files.emplace_back(lat_file+".ens");
            scalar_field_names.emplace_back("LAT");
        }

        if (parser.HasAttribute("output.ensight.conduction system.apd")) {
            std::string apd_file = parser.GetValue<std::string>("output.ensight.conduction system.apd");
            std::filesystem::path apd_path(apd_file);
            if (apd_path.has_extension()) {
                std::filesystem::path e = "";
                apd_path.replace_extension(e);
                apd_file = apd_path.string();
            }

            for (const auto &apd : post_process.APDs()) {
                cs_ens_export.SaveScalarField(apd.second.tail(react_diff->ConductSystem().NodesNum()), apd_file+std::to_string(apd.first)+".ens");
                scalar_field_files.emplace_back(apd_file+std::to_string(apd.first)+".ens");
                scalar_field_names.emplace_back("APD"+std::to_string(apd.first));
            }
        }

        if (parser.HasAttribute("output.ensight.conduction system.animation")) {
            anim_file = parser.GetValue<std::string>("output.ensight.conduction system.animation");
            std::string anim_ext = std::filesystem::path(anim_file).extension();
            if (anim_ext != ".case") { anim_file += ".case"; }

            // Get states file directory path.
            std::filesystem::path cs_states_path(states_file);
            std::string states_dir = cs_states_path.parent_path().string();

            // Look up for animation file path.
            std::filesystem::path cs_anim_path(anim_file);
            if (!cs_anim_path.has_root_path())
                anim_file = states_dir + anim_file;

            cs_ens_export.SaveAnimation(anim_file, geo_file, states_file, scalar_field_files, scalar_field_names, states_num, react_diff->OutputSteps()*react_diff->Dt());
        }
    }
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputToAscii(const Parser &parser, const ELECTRA::PostProcess &post_process, std::ostream &stream) const
{
    ELECTRA::AsciiExporter ascii_export;

    if (parser.HasAttribute("output.ascii.action potential")) {
        std::string ap_file = parser.GetValue<std::string>("output.ascii.action potential");
        std::string ap_ext = std::filesystem::path(ap_file).extension();
        if (ap_ext != ".txt") { ap_file += ".txt"; }

        ascii_export.WriteActionPotentials(post_process.ActionPotentials(), ap_file);
        stream << ELECTRA::Logger::Message("Saved action potential: "+ap_file+"\n");
    }

}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::OutputToBinary(const Parser &parser, const std::shared_ptr<ELECTRA::ReactionDiffusion<DIM, CELL_NODES>> &react_diff, std::ostream &stream) const
{
    ELECTRA::BinaryExporter bin_export;
    if (parser.HasAttribute("output.binary.cells state")) {
        std::string cell_state_file = parser.GetValue<std::string>("output.binary.cells state");
        std::string cell_state_ext = std::filesystem::path(cell_state_file).extension();
        if (cell_state_ext != ".elc") { cell_state_file += ".elc"; }

        bin_export.WriteCellsState(react_diff->Cells(), cell_state_file);
        stream << ELECTRA::Logger::Message("Saved cells state: " + cell_state_file + "\n");
    }
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::SaveFibers(const Eigen::MatrixXd &fibers, const std::string &output_filename) const
{
    // Create the path's directory if it doesn't exist.
    std::filesystem::path path(output_filename);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) { std::filesystem::create_directories(path.parent_path()); }

    // Initialize the output with proper extension.
    std::ofstream output(output_filename, std::ios::out | std::ios::trunc);
    if(!output.good()) {
        std::string error_msg = "ConfigOutput failed to write fibers in file. Check given file: " + output_filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the fibers file.
    output << "##\n";
    output << "# Created by: ElectraPre\n";
    output << "# Format: Fiber X coordinate, Fiber Y coordinate, Fiber Z coordinate\n";
    output << "##\n";
    output << std::setprecision(15);
    for (Eigen::Index i = 0; i != fibers.rows(); ++i)
        output << fibers.coeff(i,0) << ", " << fibers.coeff(i,1) << ", " << fibers.coeff(i,2) << "\n";
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::SaveDistanceField(const Eigen::VectorXd &distance_field, const std::string &output_filename) const
{
    // Create the path's directory if it doesn't exist.
    std::filesystem::path path(output_filename);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) { std::filesystem::create_directories(path.parent_path()); }

    // Initialize the output with proper extension.
    std::ofstream output(output_filename, std::ios::out | std::ios::trunc);
    if(!output.good()) {
        std::string error_msg = "ConfigOutput failed to write distance field in file. Check given file: " + output_filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the distance file.
    output << "##\n";
    output << "# Created by: ElectraPre\n";
    output << "# Distance scalar field\n";
    output << "# Format: A scalar value per node\n";
    output << "##\n";
    output << std::setprecision(15);
    for (Eigen::Index i = 0; i != distance_field.rows(); ++i)
        output << distance_field.coeff(i) << "\n";
}


template<short DIM, short CELL_NODES>
void ConfigOutput<DIM, CELL_NODES>::SaveDirectionField(const Eigen::MatrixXd &direction_field, const std::string &output_filename) const
{
    // Create the path's directory if it doesn't exist.
    std::filesystem::path path(output_filename);
    if (path.has_parent_path() && !std::filesystem::exists(path.parent_path())) { std::filesystem::create_directories(path.parent_path()); }

    // Initialize the output with proper extension.
    std::ofstream output(output_filename, std::ios::out | std::ios::trunc);
    if(!output.good()) {
        std::string error_msg = "ConfigOutput failed to write direction field in file. Check given file: " + output_filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the fibers file.
    output << "##\n";
    output << "# Created by: ElectraPre\n";
    output << "# Direction vector field\n";
    output << "# Format: A " << DIM << "-D vector per node\n";
    output << "##\n";
    output << std::setprecision(15);
    for (Eigen::Index i = 0; i != direction_field.rows(); ++i) {
        for (Eigen::Index j = 0; j != direction_field.cols()-1; ++j)
            output << direction_field.coeff(i,j) << ", ";
        output << direction_field.coeff(i,direction_field.cols()-1) << "\n";
    }
}


} // end of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_OUTPUT_TPP_