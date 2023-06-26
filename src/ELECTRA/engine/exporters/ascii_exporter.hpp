/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


/**
   \file ascii_exporter.hpp
   \brief AsciiExporter class header file.
   \author Konstantinos A. Mountris
   \date 06/05/2019
*/

#ifndef ELECTRA_EXPORTERS_ASCII_EXPORTER_HPP_
#define ELECTRA_EXPORTERS_ASCII_EXPORTER_HPP_

#include "ELECTRA/engine/utilities/logger.hpp"

#include <boost/filesystem.hpp>

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>


namespace ELECTRA {

/**
 *  \addtogroup Exporters
 *  \{
 */


/**
 * \class AsciiExporter
 * \brief Class implemmenting the output of models and their solution fields in ascii format.
 */

class AsciiExporter {

public:
    /**
     * \brief AsciiExporter constructor.
     */
    AsciiExporter();
    
    
    /**
     * \brief AsciiExporter destructor.
     */
    virtual ~AsciiExporter();
    
    
    /**
     * \brief
     */
    void WriteActionPotentials(const std::vector<Eigen::MatrixXd> &action_potentials, const std::string &filename);


    void WriteActionPotential(const Eigen::MatrixXd &action_potential, const std::string &filename);
    
};


/** \} End of Doxygen Groups*/

} // End of namespace ELECTRA.

#endif  //ELECTRA_EXPORTERS_ASCII_EXPORTER_HPP_
