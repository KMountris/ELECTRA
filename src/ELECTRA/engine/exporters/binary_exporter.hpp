/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


/**
   \file bin_exporter.hpp
   \brief BinaryExporter class header file.
   \author Konstantinos A. Mountris
   \date 20/06/2020
*/

#ifndef ELECTRA_EXPORTERS_BINARY_EXPORTER_HPP_
#define ELECTRA_EXPORTERS_BINARY_EXPORTER_HPP_

#include "ELECTRA/engine/electrophysiology/ep_factory.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>


namespace ELECTRA {

/** \addtogroup Exporters \{ */

/**
 * \class BinaryExporter
 * \brief Class implemmenting output in binary format.
 */
class BinaryExporter {

public:
    /**
     * \brief BinaryExporter constructor.
     */
    BinaryExporter();


    /**
     * \brief BinaryExporter destructor.
     */
    virtual ~BinaryExporter();


    void WriteCellsState(const std::vector<std::unique_ptr<EpBasic>> &cells, const std::string &filename);

};


/** \} End of Doxygen Groups*/

} // End of namespace ELECTRA.

#endif  //ELECTRA_EXPORTERS_BINARY_EXPORTER_HPP_
