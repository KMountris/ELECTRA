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

#include "ELECTRA/ELECTRA"

#include "tools/parser.hpp"
#include "tools/config_sim.hpp"

#include <termcolor/termcolor.hpp>
#include <nlohmann/json.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>
#include <algorithm>

using json = nlohmann::json;
using namespace ELECTRA;
using namespace APP_ELECTRA;

int main(int argc, char *argv[]) {
    try {

        std::cout << termcolor::bold << "\nWelcome to: \n";
        std::cout << termcolor::green << "   _____ _           _            "      << termcolor::magenta << " _____ _ \n";
        std::cout << termcolor::green << "  |  ___| |         | |           "      << termcolor::magenta << "/  ___(_) \n";
        std::cout << termcolor::green << "  | |__ | | ___  ___| |_ _ __ __ _"      << termcolor::magenta << "\\ `--. _ _ __ ___ \n";
        std::cout << termcolor::green << "  |  __|| |/ _ \\/ __| __| '__/ _` |"    << termcolor::magenta << "`--. \\ | '_ ` _ \\ \n";
        std::cout << termcolor::green << "  | |___| |  __/ (__| |_| | | (_| "      << termcolor::magenta << "/\\__/ / | | | | | | \n";
        std::cout << termcolor::green << "  \\____/|_|\\___|\\___|\\__|_|  \\__,_" << termcolor::magenta << "\\____/|_|_| |_| |_| \n";
        std::cout << "\n" << termcolor::reset;
        std::cout << termcolor::bold << "version: " << termcolor::reset << ELECTRA_VERSION << "\n";
        std::cout << termcolor::bold << "licence: " << termcolor::reset << "all rights reserved\n";
        std::cout << termcolor::bold << "author: " << termcolor::reset << "Konstantinos A. Mountris\n";
        std::cout << termcolor::bold << "email: " << termcolor::reset << "konstantinos.mountris@gmail.com\n";
        std::cout << termcolor::bold << "web: " << termcolor::reset << "https://www.mountris.org\n";
        std::cout << "\n";

        // Get electra configuration file.
        std::string electra_filename = "";
        if (argc == 1) {
            std::cout << termcolor::yellow << Logger::Warning("Configuration file was not provided during launching the ElectraSim app.\n"
                                                              "                  Give path to configuration filename.") << termcolor::reset;
            // Read configuration file given by the user.
            std::cout << "Path: ";
            std::cin >> electra_filename;
        }
        else { electra_filename = argv[1]; }

        // Parse input json file.
        Parser sim_parser(electra_filename);

        // Set up simulation for given scale.
        std::string sim_scale = sim_parser.GetValue<std::string>("simulation.scale");
        std::transform(std::begin(sim_scale), std::end(sim_scale), std::begin(sim_scale), ::tolower);
        if (sim_scale == "cell") {

        } else if (sim_scale == "tissue") {
            short dims = sim_parser.GetValue<short>("tissue.geometry.dimensions");
            short cell_nodes = sim_parser.GetValue<short>("tissue.geometry.cell vertices");

            if (dims == 1 && cell_nodes == 2) {
                ConfigSim<1, 2> config_sim; config_sim.Tissue(sim_parser, std::cout);
            } else if (dims == 2 && cell_nodes == 2) {
                ConfigSim<2, 2> config_sim; config_sim.Tissue(sim_parser, std::cout);
            } else if (dims == 2 && cell_nodes == 3) {
                ConfigSim<2, 3> config_sim; config_sim.Tissue(sim_parser, std::cout);
            } else if (dims == 2 && cell_nodes == 4) {
                ConfigSim<2, 4> config_sim; config_sim.Tissue(sim_parser, std::cout);
            } else if (dims == 3 && cell_nodes == 2) {
                ConfigSim<3, 2> config_sim; config_sim.Tissue(sim_parser, std::cout);
            } else if (dims == 3 && cell_nodes == 3) {
                ConfigSim<3, 3> config_sim; config_sim.Tissue(sim_parser, std::cout);
            } else if (dims == 3 && cell_nodes == 4) {
                ConfigSim<3, 4> config_sim; config_sim.Tissue(sim_parser, std::cout);
            } else if (dims == 3 && cell_nodes == 8) {
                ConfigSim<3, 8> config_sim; config_sim.Tissue(sim_parser, std::cout);
            }
        }
        else {
            throw std::invalid_argument(Logger::Error("Invalid simulation scale in configuration file. Expected: cell | tissue"));
        }

        std::cout << termcolor::magenta << termcolor::bold;
        std::cout << Logger::Message("The simulation finished successfully. Thank you for using the ElectraSim app.\n") << termcolor::reset;


    } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (const std::runtime_error &e) {
        std::cerr << "Runtime error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (const std::out_of_range &e) {
        std::cerr << "Out of Range error: " << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (const std::bad_alloc &e) {
        std::cerr << "Bad allocation error:" << termcolor::red << e.what() << termcolor::reset << std::endl;
    } catch (...) {
        std::cerr << termcolor::red << "Unknown exception..." << termcolor::reset << std::endl;
    }

    return EXIT_SUCCESS;
}