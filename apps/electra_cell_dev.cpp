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

#include <IMP/IMP>

#include <iostream>
#include <string>
#include <fstream>

using namespace ELECTRA;

int main() {
    try {

        std::string output_file = "/home/mood/Desktop/cell_test/gaur_cxx";

        // Time recorder.
        Timer timer;

        // Set measure units.
        MeasureUnits units;
        units.SetRefTimeSIValue(1.e-3);          // time:        ms
        units.SetRefLengthSIValue(1.e-2);        // length:      cm
        units.SetRefConductanceSIValue(1.e-3);   // conductance: mS
        units.SetRefCapacitanceSIValue(1.e-9);  // capacitance: nF

        // Set cell capacitance and time step
        double dt = 0.01*units["ms"];

        std::cout << Logger::Message("Single cell simulation with timestep: " + std::to_string(dt) + " ms") << std::endl;

        // Create the cell model.
        Gaur2021 cell;
        cell.Initialize(CellType::ventricular);

        // Simulation time and steps.
        int total_time = 3000.*units["ms"];
        int steps = static_cast<int>(std::ceil(total_time/dt));

        std::vector<double> cycle_lengths({630*units["ms"], 630*units["ms"], 620*units["ms"], 620*units["ms"],
                                           610*units["ms"], 610*units["ms"], 600*units["ms"], 600*units["ms"],
                                           590*units["ms"], 590*units["ms"], 580*units["ms"], 580*units["ms"],
                                           570*units["ms"], 570*units["ms"], 560*units["ms"], 560*units["ms"],
                                           550*units["ms"], 550*units["ms"], 540*units["ms"], 540*units["ms"],
                                           500*units["ms"], 500*units["ms"], 500*units["ms"], 500*units["ms"],
                                           500*units["ms"], 500*units["ms"], 500*units["ms"], 500*units["ms"],
                                           500*units["ms"], 500*units["ms"], 540*units["ms"], 540*units["ms"],
                                           550*units["ms"], 550*units["ms"], 560*units["ms"], 560*units["ms"],
                                           570*units["ms"], 570*units["ms"], 580*units["ms"], 580*units["ms"],
                                           590*units["ms"], 590*units["ms"], 600*units["ms"], 600*units["ms"],
                                           610*units["ms"], 610*units["ms"], 620*units["ms"], 620*units["ms"],
                                           630*units["ms"], 630*units["ms"]});

        Stimulus stimulus;
        stimulus.SetStart(0.*units["ms"]);
        stimulus.SetDuration(0.5*units["ms"]);
        // stimulus.SetCycleLengths(cycle_lengths);
        stimulus.SetCycleLengths(1000.);
        stimulus.SetAmplitude(80);
        std::ofstream out(output_file+".txt", std::ofstream::out);
        std::ofstream out_stim(output_file+"_cur.txt", std::ofstream::out);
        
        // Store init state

        timer.Reset();
        double stim_current = 0.;
        for (int i = 1; i <= steps; ++i) {

            if (stimulus.IsActive(i*dt, total_time)) { 
                stim_current = stimulus.Amplitude();
            } else {
                stim_current = 0.;
            }

            cell.Compute(cell.V(), dt, stim_current);
            cell.SetV(ALGORITHM::ForwardEuler(cell.V(), dt, cell.dVdt()));

            // if (i==1) {
            //     std::cout <<"step: " << i << "\n" << cell.dVdt() << "\n";
            //     std::cout << cell.PrintCurrents() << "\n";
            //     break;
            // }

            // Store new state.
            out << i*dt << " " << std::setprecision(15) << cell.V() << std::endl;
            out_stim << i*dt << " " << std::setprecision(15) << cell.Current(1) << std::endl;
        }

        // Close file.
        out.close();
        out_stim.close();
        
        // Finish program.
        std::cout << Logger::Message("Elapsed time for cell compute: " + timer.PrintElapsedTime()) << std::endl;
        std::cout << Logger::Message("Stored in: " + timer.PrintElapsedTime()) << std::endl;


    }
    catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument error: " << e.what() << std::endl;
    }
    catch (const std::runtime_error &e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
    }
    catch (const std::out_of_range &e) {
        std::cerr << "Out of Range error: " << e.what() << std::endl;
    }
    catch (const std::bad_alloc &e) {
        std::cerr << "Bad allocation error:" << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception..." << std::endl;
    }

    return 0;
}
