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

#include <termcolor/termcolor.hpp>
#include <nlohmann/json.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <memory>


using json = nlohmann::json;
using namespace ELECTRA;

int main(int argc, char *argv[]) {
    try {

        std::cout << termcolor::bold << "\nWelcome to: \n"; 
        std::cout << termcolor::green << "   _____ _           _            "      << termcolor::magenta << " _____         _    \n";
        std::cout << termcolor::green << "  |  ___| |         | |           "      << termcolor::magenta << "| ___ \\       | |   \n";
        std::cout << termcolor::green << "  | |__ | | ___  ___| |_ _ __ __ _"      << termcolor::magenta << "| |_/ /__  ___| |_  \n";
        std::cout << termcolor::green << "  |  __|| |/ _ \\/ __| __| '__/ _` "     << termcolor::magenta << "|  __/ _ \\/ __| __| \n";
        std::cout << termcolor::green << "  | |___| |  __/ (__| |_| | | (_| "      << termcolor::magenta << "| | | (_) \\__ \\ |_  \n";
        std::cout << termcolor::green << "  \\____/|_|\\___|\\___|\\__|_|  \\__,_" << termcolor::magenta << "\\_|  \\___/|___/\\__| \n";
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
            std::cout << termcolor::yellow << Logger::Warning("Configuration file was not provided during launching the ELECTRA-console app.\n"
                                                              "                  Give path to configuration filename.") << termcolor::reset;
            // Read configuration file given by the user.
            std::cout << "Path: ";
            std::cin >> electra_filename;
        }
        else { electra_filename = argv[1]; }

        std::ifstream input(electra_filename);
        json electra_json = json::parse(input);
        std::cout << std::setw(4) << electra_json << "\n\n";

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