/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file maleckar2009.hpp
   \brief Maleckar2009 class header file.
   \author Konstantinos A. Mountris
   \date 02/02/2021
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_MALECKAR2009_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_MALECKAR2009_HPP_

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
 * \namespace ELECTRA::MlcrVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Maleckar 2009 human atrial electrophysiology model variables. 
 */
namespace MlcrVar {
enum { v,           /**< Cell membrane action potential */
       dvdt,
       Na_c,
       Na_i,
       m,
       h1,
       h2,
       Ca_d,
       d_L,
       f_L1,
       f_L2,
       K_c,
       K_i,
       r,
       s,
       a_ur,
       i_ur,
       n,
       pa,
       Ca_c,
       Ca_i,
       O_C,
       O_TC,
       O_TMgC,
       O_TMgMg,
       O,
       Ca_rel,
       Ca_up,
       O_Calse,
       F1,
       F2
      };

} // End of namespace MlcrVar


/**
 * \namespace ELECTRA::MlcrPrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Maleckar 2009 human atrial electrophysiology model parameters. 
 */
namespace MlcrPrm {
enum { R,
       T,
       F,
       Cm,
       P_Na,
       g_Ca_L,
       E_Ca_app,
       k_Ca,
       g_t,
       g_kur,
       g_K1,
       g_Ks,
       g_Kr,
       g_B_Na,
       g_B_Ca,
       K_NaK_K,
       i_NaK_max,
       pow_K_NaK_Na_15,
       i_CaP_max,
       k_CaP,
       K_NaCa,
       d_NaCa,
       gamma_Na,
       ACh,
       phi_Na_en,
       Vol_i,
       Vol_d,
       tau_di,
       Mg_i,
       Vol_c,
       tau_Na,
       tau_K,
       tau_Ca,
       Na_b,
       Ca_b,
       K_b,
       I_up_max,
       k_cyca,
       k_srca,
       k_xcs,
       alpha_rel,
       Vol_up,
       Vol_rel,
       r_recov,
       tau_tr,
       k_rel_i,
       k_rel_d
      };

} // End of namespace MlcrPrm


/**
 * \namespace ELECTRA::MlcrCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Maleckar 2009 human atrial electrophysiology model currents. 
 */
namespace MlcrCur {
enum { INa,
       ICaL,
       Ito,
       IKur,
       IK1,
       IKr,
       IKs,
       INab,
       ICab,
       INaK,
       ICaP,
       INaCa,
       IKACh,
       Iion
      };
} // End of namespace MlcrCur


/**
 * \class Maleckar2009
 * \author Konstantinos A. Mountris
 * \brief  Maleckar 2009, cardiac atrial electrophysiology model \cite courtemanche1998.
 */
class Maleckar2009 : public EpBasic
{

protected:

    virtual void SetDataMapping();


public:

    Maleckar2009();


    virtual ~Maleckar2009();


    /**
     * \brief Initialize the variables and parameters of the electrophysiology model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the electrophysiology model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[MlcrVar::v] = v; }


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
    inline virtual double V() const { return this->var_[MlcrVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[MlcrVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[MlcrCur::Iion]; }

};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_MALECKAR2009_HPP_
