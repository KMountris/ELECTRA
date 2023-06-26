/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_CONFIG_PRE_TPP_
#define ELECTRA_APPS_TOOLS_CONFIG_PRE_TPP_

#include "config_pre.hpp"

namespace APP_ELECTRA {


template<short DIM, short CELL_NODES>
ConfigPre<DIM, CELL_NODES>::ConfigPre()
{}


template<short DIM, short CELL_NODES>
ConfigPre<DIM, CELL_NODES>::~ConfigPre()
{}


template<short DIM, short CELL_NODES>
void ConfigPre<DIM, CELL_NODES>::CheckValid(const Parser &parser)
{
    std::string software = parser.GetValue<std::string>("application");
    std::string author = parser.GetValue<std::string>("author");
    std::string email = parser.GetValue<std::string>("email");
    std::string licence = parser.GetValue<std::string>("license");

    if (software.find("ElectraPre") == std::string::npos ||
        author != "Konstantinos A. Mountris" ||
        email != "konstantinos.mountris@gmail.com" ||
        licence != "all rights reserved") {
        throw std::invalid_argument(Logger::Error("Header info in configuration file is not consistent with ElectraPre."));
    }

    std::string version = software.substr(software.find("v")+1);
    if (version != ELECTRA_VERSION) {
        std::string error_str = "Could not start preprocessing. The configuration file version ["
                                + version + "] does not match the ElectraPre version [" + ELECTRA_VERSION + "]";
        throw std::runtime_error(Logger::Error(error_str));
    }
}


template<short DIM, short CELL_NODES>
void ConfigPre<DIM, CELL_NODES>::Preprocess(const Parser &parser, std::ostream &stream)
{
    // Compute fibers orientation.
    if (parser.HasAttribute("fibers")) { this->Fibers(parser, stream); }
}


template<short DIM, short CELL_NODES>
void ConfigPre<DIM, CELL_NODES>::Fibers(const Parser &parser, std::ostream &stream)
{
    // Check validity of the provided configuration file.
    this->CheckValid(parser);

    // Set up model geometry.
    Timer timer;
    stream << termcolor::green << termcolor::bold << "[*** GEOMETRY SET UP ***]\n" << termcolor::reset;
    ConfigGeo<DIM, CELL_NODES> config_geo;
    IMP::Mesh<DIM, CELL_NODES> mesh;
    IMP::Grid<DIM, CELL_NODES> grid;
    IMP::Voronoi<DIM> voro;
    config_geo.SetGeometryModel(parser, mesh, grid, voro, stream);
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up numerical approximation.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** NUMERICAL APPROXIMATION SET UP ***]\n" << termcolor::reset;
    ConfigApproximation<DIM, CELL_NODES> config_approx;
    CLOUDEA::Fpm<DIM> fpm_approx;
    std::string method = parser.GetValue<std::string>("numerical approximation.method");
    std::transform(std::begin(method), std::end(method), std::begin(method), ::tolower);
    if (method == "fem") {
        config_approx.SetFemApproximation(stream);
    } else if (method == "fpm") {
        config_approx.SetFpmApproximation(parser, voro, fpm_approx, stream);
    } else {
        throw std::invalid_argument(ELECTRA::Logger::Error("Invalid simulation method in configuration file. Expected: FEM | FPM"));
    }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up fiber preprocessing.
    stream << termcolor::green << termcolor::bold << "[*** CARDIAC FIBERS COMPUTATION ***]\n" << termcolor::reset;
    timer.Reset();
    ConfigFibers<DIM, CELL_NODES> config_fibers;
    bool ventri = false;
    if (parser.HasAttribute("fibers.ventricular tags")) ventri = true;

    // Set rules and tags of the geometry.
    AtriFiberRules atri_rules;
    AtriTags atri_tags;
    VentriFiberRules ventri_rules;
    VentriTags ventri_tags;
    if (ventri) {
        if (parser.HasAttribute("fibers.ventricular rules"))
            config_fibers.SetVentriFiberRules(parser, ventri_rules, stream);
        config_fibers.SetVentriTags(parser, ventri_tags, stream);
    } else {
        if (parser.HasAttribute("fibers.atrial rules"))
            config_fibers.SetAtriFiberRules(parser, atri_rules, stream);
        config_fibers.SetAtriTags(parser, atri_tags, stream);
    }

    // Compute fibers direction.
    Ldrbm<DIM, CELL_NODES> ldrbm;
    if (ventri) {
        config_fibers.ComputeVentriFibers(parser, mesh, voro, fpm_approx, ventri_tags, ventri_rules, ldrbm, stream);
    } else {
        config_fibers.ComputeAtriFibers(parser, mesh, voro, fpm_approx, atri_tags, atri_rules, ldrbm, stream);
    }

    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

    // Set up output.
    timer.Reset();
    stream << termcolor::green << termcolor::bold << "[*** OUTPUT ***]\n" << termcolor::reset;
    ConfigOutput<DIM, CELL_NODES> config_output;
    if (method == "fem" || method == "fpm") {
        config_output.OutputFibers(parser, mesh, voro, ldrbm, stream);
    }
    stream << ELECTRA::Logger::Message("Elapsed time: ") << termcolor::cyan << termcolor::bold << timer.PrintElapsedTime() << termcolor::reset << "\n\n";

}


} // End of namespace APP_ELECTRA

#endif //ELECTRA_APPS_TOOLS_CONFIG_PRE_TPP_