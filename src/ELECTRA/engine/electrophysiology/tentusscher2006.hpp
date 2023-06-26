/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file tentusscher2006.hpp
   \brief TenTusscher2006 class header file.
   \author Konstantinos A. Mountris
   \date 17/01/2021
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_TENTUSSCHER2006_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_TENTUSSCHER2006_HPP_

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
 * \namespace ELECTRA::Tnt06Var
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Tentusscher et al 2006 human ventricular electophysiology model variables. 
 */
namespace Tnt06Var {
enum { v,          /**< membrane potential */
       dvdt,       /**< membrane potential derivative */
       ki,
       nai,
       cai,
       Xr1,
       Xr2,
       Xs,
       m,
       h,
       j,
       ca_ss,
       d,
       f,
       f2,
       fCass,
       s,
       r,
       ca_SR,
       Rprime
      };
} // End of namespace Tnt06Var


/**
 * \namespace ELECTRA::Tnt06Prm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Tentusscher et al 2006 human ventricular electophysiology model parameters. 
 */
namespace Tnt06Prm {
enum { R,
       T,
       F,
       Cm,
       Vc,
       Pkna,
       ko,
       nao,
       cao,
       gK1,
       gKr,
       gKs,
       gNa,
       gNab,
       gCaL,
       gCab,
       gto,
       PNaK,
       Kmk,
       Kmna,
       Knaca,
       Ksat,
       alpha,
       gamma,
       Kmca,
       Kmnai,
       gCap,
       Kcap,
       gKp,
       k1_prime,
       k2_prime,
       k3,
       k4,
       EC,
       max_sr,
       min_sr,
       Vrel,
       Vxfer,
       Kup,
       Vleak,
       Vmax_up,
       Buf_c,
       K_buf_c,
       Buf_sr,
       K_buf_sr,
       Buf_ss,
       K_buf_ss,
       Vsr,
       Vss,
       delta
      };
} // End of namespace Tnt06Prm


/**
 * \namespace ELECTRA::Tnt06Cur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Tentusscher 2006 human ventricular electophysiology model currents. 
 */
namespace Tnt06Cur {
enum { INaK,
       INa,
       INab,
       INaCa,
       IK1,
       Ito,
       IKr,
       IKs,
       ICaL,
       ICab,
       IKp,
       ICap,
       Iion
      };
} // End of namespace Tnt06Cur


/**
 * \class TenTusscher2006
 * \author Konstantinos A. Mountris
 * \brief  Tentusscher et al 2006, human cardiac ventricular electrophysiology model \cite ten2006.
 */

class TenTusscher2006 : public EpBasic
{
protected:

    virtual void SetDataMapping();

public:

    /**
     * \brief The default constructor of the TenTusscher2006 class.
     */
    TenTusscher2006();


    /**
     * \brief The default destructor of the TenTusscher2006 class.
     */
    virtual ~TenTusscher2006();


    /**
     * \brief Initialize the variables and parameters of the Tentusscher et al 2006 electophysiology model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Tentusscher et al 2006 electrophysiology model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[Tnt06Var::v] = v; }


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
    inline virtual double V() const { return this->var_[Tnt06Var::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[Tnt06Var::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[Tnt06Cur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_TENTUSSCHER2006_HPP_