/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file gaur.hpp
   \brief Gaur2021 class header file.
   \author Konstantinos A. Mountris
   \date 05/12/2022
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_GAUR_2021_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_GAUR_2021_HPP_

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
 * \namespace ELECTRA::Gaur21Var
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Gaur '21 porcine ventricular action potential model variables.
 */
namespace Gaur21Var {
enum { v,           /**< Cell membrane action potential */
       dvdt,
       CICR__A,
       CICR__Jrel1,
       CICR__Jrel2,
       ionic_concentrations__cacsr,
       ionic_concentrations__cai,
       ionic_concentrations__cai2,
       ionic_concentrations__cajsr,
       ionic_concentrations__cass,
       CICR__tjsrol,
       CaMK__CaMKt,
       ionic_concentrations__ki,
       ionic_concentrations__nai,
       ICaL__d,
       ICaL__fca,
       ICaL__ff,
       ICaL__fs,
       ionic_concentrations__kss,
       ionic_concentrations__nass,
       IKr__xr,
       IKs__xs1,
       IKs__xs2,
       INaL__hl,
       INaL__ml,
       ITo__aa,
       I_Na__h,
       I_Na__j,
       I_Na__m,
       ionic_concentrations__cansr
      };
} // End of namespace Gaur21Var


/**
 * \namespace ELECTRA::Gaur21Prm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Gaur '21 porcine ventricular action potential model parameters. 
 */
namespace Gaur21Prm {
enum { cell__F,  
       CaMK__KmCaMK,
       CICR__SOICR,   
       CICR__grelbarjsrol,
       CICR__tau_gap,
       CICR__tauoff,
       CICR__tauon,
       CaMK__CaMKo,
       CaMK__KmCaM,
       CaMK__PKNa,
       cell__R,
       cell__T,
       CaMK__aCaMK,
       CaMK__bCaMK,
       cell__cli,
       cell__clo,
       cell__ko,
       cell__nao,
       ICaL__PCa,
       ICaL__vhalf_d,
       ICaL__zca,
       IKb__GKb,
       IKs__GKs,
       INaCa_i__KmCaAct,
       INaCa_i__kasymm,
       INaCa_i__kcaoff,
       INaCa_i__kcaon,
       INaCa_i__kna1,
       INaCa_i__kna2,
       INaCa_i__kna3,
       INaCa_i__qca,
       INaCa_i__qna,
       INaCa_i__wca,
       INaCa_i__wna,
       INaCa_i__wnaca,
       INaCa_i__zna,
       INaK__H,
       INaK__Khp,
       INaK__Kki,
       INaK__Kko,
       INaK__Kmgatp,
       INaK__Knai0,
       INaK__Knao0,
       INaK__Knap,
       INaK__Kxkur,
       INaK__MgADP,
       INaK__MgATP,
       INaK__eP,
       INaK__k1m,
       INaK__k1p,
       INaK__k2m,
       INaK__k2p,
       INaK__k3m,
       INaK__k3p2,
       INaK__k4m,
       INaK__k4p2,
       INaK__zk,
       INaL__GNaL,
       INaL__tau_hl,
       INab__PNab,
       I_Na__GNa,
       IpCa__GpCa,
       SR_uptake__BSLmax,
       SR_uptake__BSRmax,
       SR_uptake__KmBSL,
       SR_uptake__KmBSR,
       SR_uptake__cmdnmax,
       SR_uptake__csqnmax,
       SR_uptake__kmcmdn,
       SR_uptake__kmcsqn,
       SR_uptake__kmtrpn,
       SR_uptake__trpnmax,
       cell__L,
       cell__pi,
       cell__rad,
       cell__vmyo1frac,
       stimulus__duration,
       stimulus__offset,
       stimulus__period,
       CaMK__ECl,
       ICaL__PCaK,
       ICaL__PCaNa,
       ICaL__vhalff,
       ICab__PCab,
       IK1__GK1,
       IKr__GKr,
       INaK__delta,
       INaCa_i__Gncx,
       ITo__Gto,
       stimulus__stimulus_amplitude,
       INaCa_i__h10,
       INaCa_i__h11,
       INaCa_i__h12,
       INaCa_i__k2,
       INaCa_i__k5,
       INaCa_ss__h101,
       INaCa_ss__h1111,
       INaCa_ss__h121,
       INaCa_ss__k21,
       INaCa_ss__k51,
       INaK__Pnak,
       INaK__a2,
       INaK__a4,
       INaK__b1,
       cell__Ageo,
       cell__Acap,
       cell__cao,
       INaCa_i__k1,
       INaCa_ss__k11,
       cell__vcell,
       cell__vjsr,
       cell__vcsr,
       cell__vmyo,
       cell__vmyo1,
       cell__vmyo2,
       cell__vnsr,
       cell__vnsr1,
       cell__vnsr2,
       cell__vss
      };
} // End of namespace Gaur21Prm


/**
 * \namespace ELECTRA::Gaur21Cur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Gaur '21 porcine ventricular action potential model currents.
 */
namespace Gaur21Cur {
enum { INa,
       INaL,
       ICaL,
       ICaNa,
       ICaK,
       IKr,
       IKs,
       IK1,
       ITo,
       INaCa_i,
       INaCa_ss,
       INaK,
       INab,
       IKb,
       IpCa,
       ICab,
       Iion
      };
} // End of namespace Gaur21Cur


/**
 * \class Gaur2021
 * \author Konstantinos A. Mountris
 * \brief  Gaur '21 porcine cardiac ventricular action potential model \cite gaur2021computational.
 */

class Gaur2021 : public EpBasic
{
protected:

    virtual void SetDataMapping();

public:

    /**
     * \brief The default constructor of the Gaur2021 class.
     */
    Gaur2021();


    /**
     * \brief The default destructor of the Gaur2021 class.
     */
    virtual ~Gaur2021();


    /**
     * \brief Initialize the variables and parameters of the Gaur2021 Rudy '11 cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Gaur2021 Rudy '11 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[Gaur21Var::v] = v; }


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
    inline virtual double V() const { return this->var_[Gaur21Var::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[Gaur21Var::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[Gaur21Cur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_GAUR_2021_HPP_