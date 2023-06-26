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
   \date 17/11/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_GRANDI_ATRI_CELL_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_GRANDI_ATRI_CELL_HPP_

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
 * \namespace ELECTRA::GrdAtrVar
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Grandi '11 human atrial action potential model variables. 
 */
namespace GrdAtrVar {
enum { v,           /**< Cell membrane action potential */
       dvdt,
       mo,
       ho,
       jo,
       d_o,
       fo,
       fcaBjo,
       fcaBslo,
       xtoso,
       ytoso,
       xtofo,
       ytofo,
       xkro,
       xkso,
       RyRro,
       RyRoo,
       RyRio,
       NaBjo,
       NaBslo,
       TnCLo,
       TnCHco,
       TnCHmo,
       CaMo,
       Myoco,
       Myomo,
       SRBo,
       SLLjo,
       SLLslo,
       SLHjo,
       SLHslo,
       Csqnbo
      };

} // End of namespace GrdAtrVar


/**
 * \namespace ELECTRA::GrdAtrPrm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Grandi '11 human atrial action potential model parameters. 
 */
namespace GrdAtrPrm {
enum { Ca_sro,
       Najo,
       Naslo,
       Naio,
       Kio,
       Cajo,
       Caslo,
       Caio,
       rtoso,
       n41,
       C1o,
       C2o,
       C3o,
       C4o,
       C5o,
       C6o,
       C7o,
       C8o,
       C9o,
       C10o,
       C11o,
       C12o,
       C13o,
       C14o,
       C15o,
       O1o,
       rkuro,
       skuro,
       n60,
       n61,
       n62,
       ikcaoo,
       Ach,
       ISO,
       ISOdlf1muM,
       ISOdlf1nM,
       drugvector,
       AF,
       RA,
       R,
       Fdy,
       T,
       l,
       rad,
       Cmem
      };
} // End of namespace GrdAtrPrm


/**
 * \namespace ELECTRA::GrdAtrCur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Grandi '11 human atrial action potential model currents. 
 */
namespace GrdAtrCur {
enum { INa,
       INaL,
       INab,
       INaK,
       IKr,
       IKs,
       IKp,
       IKach,
       IKCa,
       Ito,
       IKur,
       IKi,
       IClb,
       ICa,
       ICaK,
       ICaNa,
       ICatot,
       Incx,
       Ipca,
       ICab,
       Iion
      };
} // End of namespace GrdAtrCur


/**
 * \class GrandiAtri
 * \author Konstantinos A. Mountris
 * \brief  Grandi '11, cardiac atrial action potential model \cite grandi2011.
 */

class GrandiAtri : public EpBasic
{
protected:

    virtual void SetDataMapping();


public:

    GrandiAtri();


    virtual ~GrandiAtri();


    /**
     * \brief Initialize the variables and parameters of the Grandi '11 cell model. 
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Grandi '11 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[GrdAtrVar::v] = v; }


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
    inline virtual double V() const { return this->var_[GrdAtrVar::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[GrdAtrVar::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[GrdAtrCur::Iion]; }

};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_GRANDI_ATRI_CELL_HPP_