/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file paci_ventri.hpp
   \brief PaciVentri class header file.
   \author Konstantinos A. Mountris
   \date 31/05/2020
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_PACI_VENTRI_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_PACI_VENTRI_HPP_

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
 * \namespace ELECTRA::PciVtrVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Paci '13 stem cell derived  human ventricular action potential model variables.
 */
namespace PciVtrVar {
enum { v,           /**< Cell membrane action potential */
       dvdt,
       Nai,
       Cai,
       m,
       h,
       j,
       d,
       f1,
       f2,
       fCa,
       Xr1,
       Xr2,
       Xs,
       Xf,
       q,
       r,
       Ca_SR,
       g
      };
} // End of namespace PciVtrVar


/**
 * \namespace ELECTRA::PciVtrPrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Paci '13 stem cell derived  human ventricular action potential model parameters. 
 */
namespace PciVtrPrm {
enum { Cm,
       R,
       T,
       F,
       Nao,
       Cao,
       Ki,
       Ko,
       PkNa,
       g_Na,
       g_CaL,
       tau_fCa,
       g_Kr,
       L0,
       Q,
       g_Ks,
       g_K1,
       g_f,
       E_f,
       g_b_Na,
       g_b_Ca,
       Km_K,
       Km_Na,
       PNaK,
       kNaCa,
       alpha,
       gamma,
       Ksat,
       KmCa,
       KmNai,
       g_PCa,
       KPCa,
       g_to,
       Vc,
       V_SR,
       a_rel,
       b_rel,
       c_rel,
       tau_g,
       Kup,
       Buf_C,
       Buf_SR,
       Kbuf_C,
       Kbuf_SR,
       VmaxUp,
       V_leak,
       E_K,
       constf2,
       V_half,
       m_inf,
       h_inf,
       j_inf,
       d_inf,
       f1_inf,
       f2_inf,
       alpha_fCa,
       Xr1_inf,
       Xr2_inf,
       Xs_inf,
       Xf_inf,
       q_inf,
       r_inf,
       g_inf,
       alpha_m,
       alpha_h,
       alpha_j,
       alpha_d,
       constf1,
       tau_f2,
       beta_fCa,
       alpha_Xr1,
       alpha_Xr2,
       alpha_Xs,
       E_Na,
       tau_Xf,
       tau_q,
       tau_r,
       const2,
       beta_m,
       beta_h,
       beta_j,
       beta_d,
       tau_f1,
       gamma_fCa,
       beta_Xr1,
       beta_Xr2,
       beta_Xs,
       E_Ks,
       tau_m,
       tau_h,
       tau_j,
       gamma_d,
       fCa_inf,
       tau_Xr1,
       tau_Xr2,
       tau_Xs,
       E_Ca,
       tau_d,
       constfCa,
       alpha_K1,
       beta_K1,
       XK1_inf,
       Cai_bufc,
       Ca_SR_bufSR,
       Xf_infinity
      };
} // End of namespace PciVtrPrm


/**
 * \namespace ELECTRA::PciVtrCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Paci '13 stem cell derived human ventricular action potential model currents. 
 */
namespace PciVtrCur {
enum { INa,
       INaK,
       INaCa,
       INab,
       ICaL,
       IK1,
       If,
       IKr,
       IKs,
       Ito,
       IPCa,
       ICab,
       Irel,
       Iup,
       Ileak,
       Iion
      };
} // End of namespace PciVtrCur


/**
 * \class PaciVentri
 * \author Konstantinos A. Mountris
 * \brief  Paci '13, stem cell-derived human cardiac ventricular action potential model \cite paci2013.
 */

class PaciVentri : public EpBasic
{
protected:

    virtual void SetDataMapping();

public:

    /**
     * \brief The default constructor of the PaciVentri class.
     */
    PaciVentri();


    /**
     * \brief The default destructor of the PaciVentri class.
     */
    virtual ~PaciVentri();


    /**
     * \brief Initialize the variables and parameters of the Paci '13 cell model. 
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Paci '13 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[PciVtrVar::v] = v; }


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
    inline virtual double V() const { return this->var_[PciVtrVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[PciVtrVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[PciVtrCur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_PACI_VENTRI_HPP_