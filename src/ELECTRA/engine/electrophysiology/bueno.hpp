/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file bueno_cell.hpp
   \brief BuenoCell class header file.
   \author Konstantinos A. Mountris
   \date 21/09/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_BUENO_CELL_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_BUENO_CELL_HPP_

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
 * \namespace ELECTRA::BueVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Bueno-Orovio '07 ventricular cell action potential model variables. 
 */
namespace BueVar {
enum { v,           /**< Cell membrane potential */
       dvdt,        /**< Cell membrane potential's time derivative */
       Vm,
       g_u,
       g_v,
       g_w,
       g_s,
       p,
       tau_s,
       m,
       v_inf,
       q,
       tau_v_minus,
       r,
       w_inf,
       tau_w_minus,
       tau_o,
       tau_so
      };
} // End of namespace BueVar


/**
 * \namespace ELECTRA::BuePrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Bueno-Orovio '07 ventricular cell action potential model parameters. 
 */
namespace BuePrm {
enum { V0,
       Vfi,
       u_m,
       u_p,
       tau_v_plus,
       tau_s1,
       k_s,
       u_s,
       tau_s2,
       tau_v2_minus,
       u_q,
       tau_o1,
       u_r,
       tau_fi,
       u_u,
       tau_so1,
       tau_so2,
       k_so,
       u_so,
       tau_si,
       tau_winf,
       wstar_inf,
       tau_w1_minus,
       tau_w2_minus,
       k_w_minus,
       u_w_minus,
       tau_w_plus,
       tau_v1_minus,
       tau_o2
     };
} // End of namespace BuePrm


/**
 * \namespace ELECTRA::BueCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Bueno-Orovio '07 ventricular cell action potential model currents. 
 */
namespace BueCur {
enum { Ifi,        /**< Fast inward current [dimensionless] */
       Iso,        /**< Slow outward current [dimensionless] */
       Isi,        /**< Slow inward current [dimensionless] */
       Iion        /**< Total ionic current [dimensionless] */
      };
} // End of namespace BueCur


/**
 * \class Bueno
 * \author Konstantinos A. Mountris
 * \brief Class implemmenting the Bueno-Orovio '07, human cardiac ventricular, cell model \cite bueno2008.
 */

class Bueno : public EpBasic
{
protected:

    virtual void SetDataMapping();


public:

    /**
     * \brief The default constructor of the Bueno-Orovio '07 cell model class.
     */
    Bueno();


    /**
     * \brief The default destructor of the Bueno-Orovio '07 cell model class.
     */
    virtual ~Bueno();


    /**
     * \brief Initialize the variables and parameters of the Bueno-Orovio '07 cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Bueno-Orovio '07 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[BueVar::v] = v; }


    /*!
     * \brief Print to std::string the cell's variables and their values.
     * \return [std::string] The cell's variables and their values.
    */
    virtual std::string PrintVariables() const;


    /*!
     * \brief Print to std::string the cell's parameters and their values.
     * \return [std::string] The cell's parameters and their values.
    */
    virtual std::string PrintParameters() const;


    /*!
     * \brief Print to std::string the cell's currents and their values.
     * \return [std::string] The cell's currents and their values.
    */
    virtual std::string PrintCurrents() const;


    /*!
     * \brief Print to std::string the cell's currents' block coefficients and their values.
     * \return [std::string] The cell's currents' block coefficients and their values.
    */
    virtual std::string PrintBlockCoeffs() const;


    /**
     * \brief Get the membrane potential value of the cell model.
     * \return [double] The membrane potential value.
     */
    inline virtual double V() const { return this->var_[BueVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline double dVdt() const { return this->var_[BueVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] The total ionic current.
     */
    inline double Iion() const { return this->cur_[BueCur::Iion]; }

};

/*! \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_BUENO_CELL_HPP_