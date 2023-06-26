/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file maccannell.hpp
   \brief Maccannell class header file.
   \author Konstantinos A. Mountris
   \date 27/11/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_MACCANNELL_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_MACCANNELL_HPP_

#include "ELECTRA/engine/electrophysiology/ep_basic.hpp"
#include "ELECTRA/engine/utilities/algorithm.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <limits>


namespace ELECTRA {

/** \addtogroup Electrophysiology \{ */


/**
 * \namespace ELECTRA::McnVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access MacCannell '07 fibroblast action potential model variables.
 */
namespace McnVar {
enum { v,           /**< Cell membrane potential */
       dvdt,        /**< Cell membrane potential's time derivative */
       rf,          /**< Gate parameter r */
       sf,          /**< Gate parameter f */
       gKv,         /**< Potassium transport [nS/pF] */
       gK1,         /**< Potassium transport [nS/pF] */
       gbna,        /**< Sodium transport [nS/pF] */
       ko,          /**< Extracellular potassium concentration [mMol] */
       ki,          /**< Intracellular potassium concentration [mMol] */
       nao,         /**< Extracellular sodium concentration [mMol] */
       nai          /**< Intracellular sodium concentration [mMol] */
      };
} // End of namespace McnVar


/**
 * \namespace ELECTRA::McnPrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access MacCannell '07 fibroblast action potential model parameters.
 */
namespace McnPrm {
enum { R,          /**< Gas constant [J/(K kmol)] */
       T,          /**< Temperature [K] */
       Fdy,        /**< Faraday constant [C/mol] */
       RTF,        /**< R*T*Fdy product [1/mV] */
       iRTF,       /**< R*T*Fdy inverse product [mV] */
       INaKmax,    /**< [pA/pF] */
       KmKf,       /**< [mmol/L] */
       KmNaf,      /**< [mmol/L] */
       Bf,         /**< [mV] */
       Vrev        /**< [mV] */
     };
} // End of namespace McnPrm


/**
 * \namespace ELECTRA::McnCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access MacCannell '07 fibroblast action potential model currents.
 */
namespace McnCur {
enum { INaK,       /**< Sodium-potassium pump current [pA/pF] */
       IKv,        /**< Time and voltage dependent potassium current [pA/pF] */
       IbNa,       /**< Background sodium current [pA/pF] */
       IK1,        /**< Inward rectifying potassium current [pA/pF] */
       Iion        /**< Total ionic current [pA/pF] */
      };
} // End of namespace McnCur


/**
 * \class Maccannell
 * \author Konstantinos A. Mountris
 * \brief MacCannell '07 fibroblast action potential model \cite maccannell2007.
 */
class Maccannell : public EpBasic
{
protected:

    virtual void SetDataMapping();


public:

    /**
     * \brief The default constructor of the Maccannel class.
     */
    Maccannell();


    /**
     * \brief The default destructor of the Maccannel class.
     */
    virtual ~Maccannell();


    /**
     * \brief Initialize the variables and parameters of the MacCannel '07 action potential model. 
     * \return [void]
     */
    void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the MacCannel '07 action potential model.
     * \return [void]
     */
    void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane action potential.
     * \param [in] v The cell membrane action potential value.
     * \return [void]
     */
    inline void SetV(double v) { this->var_[McnVar::v] = v; }


    /**
     * \brief Print to std::string the cell's variables and their values.
     * \return [std::string] The cell's variables and their values.
    */
    virtual std::string PrintVariables() const;


    /**
     * \brief Print to std::string the cell's parameters and their values.
     * \return [std::string] The cell's parameters and their values.
    */
    virtual std::string PrintParameters() const;


    /**
     * \brief Print to std::string the cell's currents and their values.
     * \return [std::string] The cell's currents and their values.
    */
    virtual std::string PrintCurrents() const;


    /**
     * \brief Print to std::string the cell's currents' block coefficients and their values.
     * \return [std::string] The cell's currents' block coefficients and their values.
    */
    virtual std::string PrintBlockCoeffs() const;


    /**
     * \brief Get the membrane potential value of the cell model.
     * \return [double] The membrane potential value.
     */
    inline double V() const { return this->var_[McnVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline double dVdt() const { return this->var_[McnVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[McnCur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_MACCANNELL_HPP_