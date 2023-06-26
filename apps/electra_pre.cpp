/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/ELECTRA"

#include "tools/parser.hpp"
#include "tools/config_pre.hpp"

#include <termcolor/termcolor.hpp>
#include <nlohmann/json.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <memory>


using json = nlohmann::json;
using namespace ELECTRA;
using namespace APP_ELECTRA;

int main(int argc, char *argv[]) {
    try {

        std::cout << termcolor::bold << "\nWelcome to: \n";
        std::cout << termcolor::green << "   _____ _           _            "      << termcolor::magenta << "______ \n";
        std::cout << termcolor::green << "  |  ___| |         | |           "      << termcolor::magenta << "| ___ \\ \n";
        std::cout << termcolor::green << "  | |__ | | ___  ___| |_ _ __ __ _"      << termcolor::magenta << "| |_/ / __ ___ \n";
        std::cout << termcolor::green << "  |  __|| |/ _ \\/ __| __| '__/ _` "     << termcolor::magenta << "|  __/ '__/ _ \\ \n";
        std::cout << termcolor::green << "  | |___| |  __/ (__| |_| | | (_| "      << termcolor::magenta << "| |  | | |  __/ \n";
        std::cout << termcolor::green << "  \\____/|_|\\___|\\___|\\__|_|  \\__,_" << termcolor::magenta << "\\_|  |_|  \\___| \n";
        std::cout << "\n" << termcolor::reset;
        std::cout << termcolor::bold << "version: " << termcolor::reset << ELECTRA_VERSION << "\n";
        std::cout << termcolor::bold << "licence: " << termcolor::reset << "To be defined\n";
        std::cout << termcolor::bold << "author: " << termcolor::reset << "Konstantinos A. Mountris\n";
        std::cout << termcolor::bold << "email: " << termcolor::reset << "konstantinos.mountris@gmail.com\n";
        std::cout << termcolor::bold << "web: " << termcolor::reset << "https://www.mountris.org\n";
        std::cout << "\n";

        // Get electra configuration file.
        std::string electra_filename = "";
        if (argc == 1) {
            std::cout << termcolor::yellow << Logger::Warning("Configuration file was not provided during launching the ElectraPre app.\n"
                                                              "                  Give path to configuration filename.") << termcolor::reset;
            // Read configuration file given by the user.
            std::cout << "Path: ";
            std::cin >> electra_filename;
        }
        else { electra_filename = argv[1]; }

        // Parse input json file.
        Parser pre_parser(electra_filename);

        // Set up preprocessing.
        short dims = pre_parser.GetValue<short>("tissue.geometry.dimensions");
        if (dims != 3) {
            std::string error_msg = "Could not perform preprocessing. Currently supported only 3D geometries.";
            throw std::runtime_error(Logger::Error(error_msg));
        }

        short cell_nodes = 1;
        if (pre_parser.HasAttribute("tissue.geometry.cell vertices"))
            cell_nodes = pre_parser.GetValue<short>("tissue.geometry.cell vertices");

        // Initialize preprocessing configurator.
        if (dims == 3 && cell_nodes == 4) {
            ConfigPre<3, 4> config_pre;  config_pre.Preprocess(pre_parser, std::cout);
        } else if (dims == 3 && cell_nodes == 8) {
            ConfigPre<3, 8> config_pre;  config_pre.Preprocess(pre_parser, std::cout);
        } else {
            std::string error_str = "Preprocessing is currently available only for 3D meshes composed of tetrahedral elements";
            throw std::invalid_argument(Logger::Error(error_str));
        }

        std::cout << termcolor::magenta << termcolor::bold;
        std::cout << Logger::Message("The preprocessing finished successfully. Thank you for using the ElectraPre app.\n") << termcolor::reset;

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