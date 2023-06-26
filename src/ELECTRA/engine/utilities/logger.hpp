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

/*!
   \file logger.hpp
   \brief Logger class header file.
   \author Konstantinos A. Mountris
   \date 07/03/2019
*/

#ifndef ELECTRA_UTILITIES_LOGGER_HPP_
#define ELECTRA_UTILITIES_LOGGER_HPP_

#include <string>

namespace ELECTRA {

/*!
 *  \addtogroup Utilities
 *  @{
 */


/*!
 * \class Logger
 * \brief Class implemmenting output messages for logging of ELECTRA.
 */

class Logger
{
public:

    /*!
     * \brief Logger default constructor.
     */
    Logger() {}


    /*!
     * \brief Logger default destructor.
     */
    virtual ~Logger() {}


    /*!
     * \brief Log an ELECTRA message.
     * \param [in] msg The message to be logged.
     * \return [std::string] The logged message.
     */
    inline static std::string Message(const std::string &msg) {
        return "[ELECTRA] " + msg;
    }


    /*!
     * \brief Log an ELECTRA error.
     * \param [in] err The error to be logged.
     * \return [std::string] The logged error.
     */
    inline static std::string Error(const std::string &err) {
        return "[ELECTRA ERROR] " + err;
    }


    /*!
     * \brief Log an ELECTRA warning.
     * \param [in] wrng The warning to be logged.
     * \return [std::string] The logged warning.
     */
    inline static std::string Warning(const std::string &wrng) {
        return "\n[ELECTRA WARNING] " + wrng + "\n";
    }


        /*!
     * \brief Log an ELECTRA ToDO message.
     * \param [in] todo The ToDo message to be logged.
     * \return [std::string] The logged ToDo message.
     */
    inline static std::string ToDo(const std::string &todo) {
        return "[ELECTRA TODO] " + todo;
    }


};

/*! @} End of Doxygen Groups*/

} //end of namespace ELECTRA

#endif //ELECTRA_UTILITIES_LOGGER_HPP_
