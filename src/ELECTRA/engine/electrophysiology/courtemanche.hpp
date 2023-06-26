/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file grandi_atri.hpp
   \brief GrandiAtri class header file.
   \author Konstantinos A. Mountris
   \date 23/01/2020
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_COURTEMANCHE_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_COURTEMANCHE_HPP_

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
 * \namespace ELECTRA::CourteVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Courtemanche '98 human atrial action potential model variables. 
 */
namespace CourteVar {
enum { v,           /**< Cell membrane action potential */
       dvdt,
       g_u,
       g_v,
       g_w,
       d,
       f_Ca,
       f,
       h,
       j,
       m,
       Ca_i,
       Ca_rel,
       Ca_up,
       K_i,
       Na_i,
       xr,
       xs,
       oa,
       oi,
       ua,
       ui,
       y
      };

} // End of namespace CourteVar


/**
 * \namespace ELECTRA::CourtePrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Courtemanche '98 human atrial action potential model parameters. 
 */
namespace CourtePrm {
enum { rkr,
       TYPE1,
       TYPE3,
       ISOdlf,
       cAF,
       cAF2,
       SR,
       Ach,
       drugvector,
       Engel,
       Penaranda,
       Cm,              /**< picoF (in membrane)  50pF for Maleckar */
       CMDN_max,        /**<   % millimolar (in Ca_buffers) */
       CSQN_max,        /**<   % millimolar (in Ca_buffers) */
       Km_CMDN,         /**<   % millimolar (in Ca_buffers) */
       Km_CSQN,         /**<   % millimolar (in Ca_buffers) */
       Km_TRPN,         /**<   % millimolar (in Ca_buffers) */
       TRPN_max,        /**<   % millimolar (in Ca_buffers) */
       Ca_up_max,       /**<   % millimolar (in Ca_leak_current_by_the_NSR) */
       K_rel,           /**<   % per_millisecond -1 (in Ca_release_current_from_JSR) */
       I_up_max,        /**<   % millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR) */
       K_up,            /**<   % millimolar (in Ca_uptake_current_by_the_NSR) */
       I_NaCa_max,      /**<   % picoA_per_picoF (in Na_Ca_exchanger_current) */
       K_mCa,           /**<   % millimolar (in Na_Ca_exchanger_current) */
       K_mNa,           /**<   % millimolar (in Na_Ca_exchanger_current) */
       K_sat,           /**<   % dimensionless (in Na_Ca_exchanger_current) */
       gamma,           /**<   % dimensionless (in Na_Ca_exchanger_current) */
       g_B_Ca,          /**< nanoS_per_picoF (in background_currents) */
       g_B_K,           /**< nanoS_per_picoF (in background_currents) */
       g_B_Na,          /**< nanoS_per_picoF (in background_currents) */
       g_Na,            /**< nanoS_per_picoF (in fast_sodium_current) */
       g_Kca_eb,        /**< nanoS_per_picoF */
       g_Kca_pb,        /**< nanoS_per_picoF */
       V_cell,          /**< micrometre_3 (in intracellular_ion_concentrations) */
       V_i,             /**< prm.V_cell*0.68;  %micrometre_3 (in intracellular_ion_concentrations) */
       V_up,            /**< 0.0552*prm.V_cell;  %micrometre_3 (in intracellular_ion_concentrations) */
       V_rel,           /**< 0.0048*prm.V_cell;  %micrometre_3 (in intracellular_ion_concentrations) */
       F,               /**< coulomb_per_millimole (in membrane) */
       R,               /**< joule_per_mole_kelvin (in membrane) */
       T,               /**< kelvin (in membrane) */
       i_CaP_max,       /**< picoA_per_picoF (in sarcolemmal_calcium_pump_current) */
       Km_K_o,          /**< millimolar (in sodium_potassium_pump) */
       Km_Na_i,         /**< millimolar (in sodium_potassium_pump) */
       i_NaK_max,       /**< picoA_per_picoF (in sodium_potassium_pump) */
       Ca_o,            /**< millimolar (in standard_ionic_concentrations) */
       K_o,             /**< millimolar (in standard_ionic_concentrations) */
       Na_o,            /**< millimolar (in standard_ionic_concentrations) */
       tau_tr,          /**< millisecond (in transfer_current_from_NSR_to_JSR) */
       K_Q10,            /**< dimensionless (in transient_outward_K_current) */
       EAD
      };

} // End of namespace CourtePrm


/**
 * \namespace ELECTRA::CourteCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Courtemanche '98 human atrial action potential model currents. 
 */
namespace CourteCur {
enum { INa,
       IK1,
       Ito,
       IKur,
       IKr,
       IKs,
       ICaL,
       INaK,
       IKach,
       IKca_E,
       IKca_P,
       INaCa,
       INab,
       ICab,
       ICap,
       IRel,
       Itr,
       IupLeak,
       Iup,
       Iion
      };
} // End of namespace CourteCur


/**
 * \class Courtemanche
 * \author Konstantinos A. Mountris
 * \brief  Courtemanche '98, cardiac atrial action potential model \cite courtemanche1998.
 */
class Courtemanche : public EpBasic
{

protected:

    virtual void SetDataMapping();


public:

    Courtemanche();


    virtual ~Courtemanche();


    /**
     * \brief Initialize the variables and parameters of the Courtemanche '98 cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Courtemanche '98 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[CourteVar::v] = v; }


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
    inline virtual double V() const { return this->var_[CourteVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[CourteVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[CourteCur::Iion]; }

};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_COURTEMANCHE_HPP_