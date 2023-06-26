/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file ohara.hpp
   \brief Ohara class header file.
   \author Konstantinos A. Mountris
   \date 13/03/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_OHARA_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_OHARA_HPP_

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
 * \namespace ELECTRA::OhrVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access O'Hara Rudy '11 human ventricular action potential model variables.
 */
namespace OhrVar {
enum { v,           /**< Cell membrane action potential */
       dvdt,
       nai,
       nass,
       ki,
       kss,
       cai,
       cass,
       cansr,
       cajsr,
       m,
       hf,
       hs,
       j,
       hsp,
       jp,
       mL,
       hL,
       hLp,
       a,
       iF,
       iS,
       ap,
       iFp,
       iSp,
       d,
       ff,
       fs,
       fcaf,
       fcas,
       jca,
       nca,
       ffp,
       fcafp,
       xrf,
       xrs,
       xs1,
       xs2,
       xk1,
       Jrelnp,
       Jrelp,
       CaMKt
      };
} // End of namespace OhrVar


/**
 * \namespace ELECTRA::OhrPrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access O'Hara Rudy '11 human ventricular action potential model parameters. 
 */
namespace OhrPrm {
enum { gNaL,
       delta_factor,
       gto,
       pCa,
       gKr,
       gKs,
       gK1,
       gncx,
       pNaK,
       Jrel_0,
       Jrelp_0,
       Jupnp_0,
       Jupp_0,
       Bcai_factor,
       nao,
       cao,
       ko,
       R,
       T,
       Fdy,
       l,
       rad,
       vcell,
       ageo,
       acap,
       vmyo,
       vnsr,
       vjsr,
       vss,
       KmCaMK,
       aCaMK,
       bCaMK,
       CaMKo,
       KmCaM,
       BSRmax,
       KmBSR,
       BSLmax,
       KmBSL,
       cmdnmax,
       Kmcmdn,
       trpnmax,
       Kmtrpn,
       Csqnmax,
       kmcsqn,
       gNa
      };
} // End of namespace OhrPrm


/**
 * \namespace ELECTRA::OhrCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access O'Hara Rudy '11 human ventricular action potential model currents.
 */
namespace OhrCur {
enum { INa,
       INaL,
       Ito,
       ICaL,
       ICaNa,
       ICaK,
       IKr,
       IKs,
       IK1,
       INaCa_i,
       INaCa_ss,
       INaCa,
       INaK,
       IKb,
       INab,
       ICab,
       IpCa,
       Iion
      };
} // End of namespace OhrCur


/**
 * \class Ohara
 * \author Konstantinos A. Mountris
 * \brief  O'hara Rudy '11, human cardiac ventricular action potential model \cite ohara2011.
 */

class Ohara : public EpBasic
{
protected:

    virtual void SetDataMapping();

public:

    /**
     * \brief The default constructor of the Ohara class.
     */
    Ohara();


    /**
     * \brief The default destructor of the Ohara class.
     */
    virtual ~Ohara();


    /**
     * \brief Initialize the variables and parameters of the Ohara Rudy '11 cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Ohara Rudy '11 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[OhrVar::v] = v; }


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
    inline virtual double V() const { return this->var_[OhrVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[OhrVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[OhrCur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_OHARA_HPP_