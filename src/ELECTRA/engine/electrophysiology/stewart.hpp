/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file stewart.hpp
   \brief Stewart class header file.
   \author Konstantinos A. Mountris
   \date 13/03/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_STEWART_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_STEWART_HPP_

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
 * \namespace ELECTRA::StrtVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Stewart '09 human purkinje action potential model variables. 
 */
namespace StrtVar {
enum { v,           /**< Cell membrane action potential */
       dvdt,
       K_i,
       Na_i,
       Ca_i,
       y,
       Xr1,
       Xr2,
       Xs,
       m,
       h,
       j,
       Ca_ss,
       d,
       f,
       f2,
       fCass,
       s,
       r,
       Ca_SR,
       R_prime,
       V_c,
       P_kna,
       K_o,
       Na_o,
       Ca_o,
       g_f_Na,
       g_f_K,
       g_f_K1,
       g_f_Kr,
       g_f_Ks,
       g_Na,
       g_bna,
       g_CaL,
       g_bca,
       g_to,
       g_sus,
       p_NaK,
       K_mk,
       K_mNa,
       K_NaCa,
       K_sat,
       alpha,
       gamma,
       Km_Ca,
       Km_Nai,
       g_pCa,
       K_pCa,
       g_pK,
       k1_prime,
       k2_prime,
       k3,
       k4,
       EC,
       max_sr,
       min_sr,
       V_rel,
       V_xfer,
       K_up,
       V_leak,
       Vmax_up,
       Buf_c,
       K_buf_c,
       Buf_sr,
       K_buf_sr,
       Buf_ss,
       K_buf_ss,
       V_sr,
       V_ss,
      };

} // End of namespace StrtVar


/**
 * \namespace ELECTRA::StrtPrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Stewart '09 human purkinje action potential model parameters. 
 */
namespace StrtPrm {
enum { R,
       T,
       Frdy,
       Cm,
       y_inf,
       xr1_inf,
       xr2_inf,
       xs_inf,
       m_inf,
       h_inf,
       j_inf,
       d_inf,
       f_inf,
       f2_inf,
       fCass_inf,
       s_inf,
       r_inf,
       E_Na,
       alpha_y,
       alpha_xr1,
       alpha_xr2,
       alpha_xs,
       alpha_m,
       alpha_h,
       alpha_j,
       alpha_d,
       tau_f,
       tau_f2,
       tau_fCass,
       tau_s,
       tau_r,
       E_K,
       beta_y,
       beta_xr1,
       beta_xr2,
       beta_xs,
       beta_m,
       beta_h,
       beta_j,
       beta_d,
       E_Ks,
       tau_y,
       tau_xr1,
       tau_xr2,
       tau_xs,
       tau_m,
       tau_h,
       tau_j,
       gamma_d,
       E_Ca,
       tau_d,
       xK1_inf,
       a,
       kcasr,
       Ca_i_bufc,
       k1,
       k2,
       O,
       Ca_sr_bufsr,
       Ca_ss_bufss
      };
} // End of namespace StrtPrm


/**
 * \namespace ELECTRA::StrtCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Stewart '09 human purkinje action potential model currents. 
 */
namespace StrtCur {
enum { IfNa,
       IfK,
       If,
       Irel,
       Isus,
       INaK,
       INaCa,
       IpCa,
       IpK,
       Iup,
       Ileak,
       Ixfer,
       IK1,
       IKr,
       IKs,
       INa,
       IbNa,
       ICaL,
       IbCa,
       Ito,
       Iion
      };
} // End of namespace StrtCur


/**
 * \class Stewart
 * \author Konstantinos A. Mountris
 * \brief  Stewart '09, human purkinje action potential model \cite stewart2009.
 */

class Stewart : public EpBasic
{
protected:

    virtual void SetDataMapping();

public:

    /**
     * \brief The default constructor of the Stewart class.
     */
    Stewart();


    /**
     * \brief The default destructor of the Stewart class.
     */
    virtual ~Stewart();


    /**
     * \brief Initialize the variables and parameters of the cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[StrtVar::v] = v; }


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
    inline virtual double V() const { return this->var_[StrtVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[StrtVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[StrtCur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_STEWART_HPP_