/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/exporters/ascii_exporter.hpp"


namespace ELECTRA
{

AsciiExporter::AsciiExporter()
{}


AsciiExporter::~AsciiExporter()
{}


void AsciiExporter::WriteActionPotentials(const std::vector<Eigen::MatrixXd> &action_potentials, const std::string &filename)
{
    std::string output_file = filename.substr(0, filename.size()-4);

    for (const auto &ap : action_potentials) {
        
        auto id = std::distance(&action_potentials[0], &ap);
        
        this->WriteActionPotential(ap, output_file + "_" + std::to_string(id+1) + ".txt");
    }

}


void AsciiExporter::WriteActionPotential(const Eigen::MatrixXd &action_potential, const std::string &filename)
{

    // Position of the last slash in the exporting file's name.
    std::size_t last_slash = filename.find_last_of("/\\");

    // Get the path directory of the exporting file name.
    std::string path = "";
    if (last_slash != std::string::npos) { path = filename.substr(0, last_slash); }

    // Create the path's directory if it doesn't exist.
    boost::filesystem::path dir(path);
    if (!path.empty() && !boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
    }

    // Initialize the output.
    std::ofstream output(filename, std::ios::out);

    // Check if file opened properly.
    if(!output.good()) {
        std::string error_msg = "AsciiExporter failed to write action potential in file. Check given file: " + filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the action potential header file.
    output << "##\n";
    output << "# ELECTRA v0.0.1\n";
    output << "# Action Potential output for a specified node.\n";
    output << "# Node id: XX\n";
    output << "# File format: t (time), AP (action potential)\n";
    output << "# Units: t - ms | AP - mV\n";
    output << "##\n";

    //Set numeric precision to max.
    output.precision(15);

    // Write the time and action potential values in the output file.
    for (int i = 0; i != action_potential.rows(); ++i) {
        output << action_potential.coeff(i,0) << ", " << action_potential.coeff(i,1) << std::endl;
    }

    // Close file.
    output << "##\n";
    output << "# END of file\n";
    output << "##\n";
    output.close();

}


} // End of namespace ELECTRA
