/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/exporters/binary_exporter.hpp"


namespace ELECTRA
{

BinaryExporter::BinaryExporter()
{}


BinaryExporter::~BinaryExporter()
{}


void BinaryExporter::WriteCellsState(const std::vector<std::unique_ptr<EpBasic>> &cells, const std::string &filename)
{
    namespace boost_fs = boost::filesystem;

    // Get the path directory of the filename.
    boost_fs::path p(filename);
    std::string path = p.parent_path().string();

    // Create the path's directory if it doesn't exist.
    boost_fs::path dir(path);
    if (!path.empty() && !boost_fs::exists(dir)) { boost_fs::create_directories(dir); }

    // Search for the extension of the file.
    std::string ext = p.extension().string();

    // Initialize the output with proper extension.
    std::ofstream output;
    if (ext == ".elc") {
        output.open(filename, std::ios::out | std::ios::binary);
    } else {
        output.open(filename+".elc", std::ios::out | std::ios::binary);
    }

    // Check if file opened properly.
    if(!output.good()) {
        std::string error_msg = "BinExporter failed to write cells state in file. Check given file: " + filename;
        throw std::runtime_error(Logger::Error(error_msg));
    }

    // Write the cells state header file.    
    std::string header = "##\n# electra\n# type: cells state\n# format: cells_num cell_type cell_vars cell_pms cell_cur cell_block_cur\n##";
    std::size_t header_size = header.size();

    output.write(reinterpret_cast<char *>(&header_size), sizeof(std::size_t));
	output.write(header.c_str(), header_size);
    
    // Write number of cells.
    std::size_t cells_num = cells.size();
    output.write(reinterpret_cast<char *>(&cells_num), sizeof(std::size_t));

    // Write the state of each of the cells' variables.
    int cell_type_id = 0;
    double val = 0.;
    for (const auto &cell : cells) {
        cell_type_id = static_cast<int>(cell->ModelType());
        output.write(reinterpret_cast<char *>(&cell_type_id), sizeof(int));

        // Write variables.
        for (int i = 0; i != cell->VarNum(); ++i) {
            val = cell->Var(i);
            output.write(reinterpret_cast<char *>(&val), sizeof(double));
        }

        // Write parameters.
        for (int i = 0; i != cell->PrmNum(); ++i) {
            val = cell->Prm(i);
            output.write(reinterpret_cast<char *>(&val), sizeof(double));
        }

        // Write currents.
        for (int i = 0; i != cell->CurrentNum(); ++i) {
            val = cell->Current(i);
            output.write(reinterpret_cast<char *>(&val), sizeof(double));
        }

        // Write block current coefficients.
        for (int i = 0; i != cell->CurrentNum()-1; ++i) {
            val = cell->BlockCoeff(i);
            output.write(reinterpret_cast<char *>(&val), sizeof(double));
        }
    }
    
    // Close the output file.
    output.close();
}


} // End of namespace ELECTRA
