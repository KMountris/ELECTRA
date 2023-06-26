/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
   \file ventri_fiber_rules.hpp
   \brief VentriFiberRules class header file.
   \author Konstantinos A. Mountris
   \date 18/03/2021
*/

#ifndef ELECTRA_FIBERS_VENTRI_FIBER_RULES_HPP_
#define ELECTRA_FIBERS_VENTRI_FIBER_RULES_HPP_

namespace ELECTRA
{

/** \addtogroup Fibers \{ */

/**
 * \class VentriFiberRules
 * \brief Class collecting angle rules for cardiac fiber determination with Ldrbm (Laplace-Dirichlet rule based method).
 * Angles are defined in degrees.
 */
class VentriFiberRules
{
private:

    double alpha_lv_endo_;          // The alpha angle value at the left ventricle endocardium

    double alpha_lv_epi_;           // The alpha angle value at the left ventricle epicardium

    double alpha_rv_endo_;          // The alpha angle value at the right ventricle endocardium

    double alpha_rv_epi_;           // The alpha angle value at the right ventricle epicardium

    double alpha_ot_lv_endo_;       // The alpha angle value at the outflow tract left ventricle endocardium

    double alpha_ot_lv_epi_;        // The alpha angle value at the outflow tract left ventricle epicardium

    double alpha_ot_rv_endo_;       // The alpha angle value at the outflow tract right ventricle endocardium

    double alpha_ot_rv_epi_;        // The alpha angle value at the outflow tract right ventricle epicardium

    double alpha_septum_;           // The alpha angle value at the septum

    double beta_lv_endo_;           // The beta angle value at the left ventricle endocardium

    double beta_lv_epi_;            // The beta angle value at the left ventricle epicardium

    double beta_rv_endo_;           // The beta angle value at the right ventricle endocardium

    double beta_rv_epi_;            // The beta angle value at the right ventricle epicardium

    double beta_ot_lv_endo_;        // The beta angle value at the outflow tract left ventricle endocardium

    double beta_ot_lv_epi_;         // The beta angle value at the outflow tract left ventricle epicardium

    double beta_ot_rv_endo_;        // The beta angle value at the outflow tract right ventricle endocardium

    double beta_ot_rv_epi_;         // The beta angle value at the outflow tract right ventricle epicardium

public:

    /**
     * \brief Default constructor of ventricular angle rules.
     */
    inline VentriFiberRules() :
        alpha_lv_endo_(-60.), alpha_lv_epi_(60.), alpha_rv_endo_(90.), alpha_rv_epi_(-25.), alpha_ot_lv_endo_(-90.),
        alpha_ot_lv_epi_(0.), alpha_ot_rv_endo_(90.), alpha_ot_rv_epi_(0.), alpha_septum_(0.), beta_lv_endo_(-20.), beta_lv_epi_(20.),
        beta_rv_endo_(0.), beta_rv_epi_(20.), beta_ot_lv_endo_(0.), beta_ot_lv_epi_(0.), beta_ot_rv_endo_(0.), beta_ot_rv_epi_(0.) {}


    /**
     * \brief Default destructor of ventricular angle rules.
     */
    inline virtual ~VentriFiberRules() {}


    /**
     * \brief Set the alpha angle value at the left ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaLvEndo(double degrees) { this->alpha_lv_endo_ = degrees; }


    /**
     * \brief Set the alpha angle value at the left ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaLvEpi(double degrees) { this->alpha_lv_epi_ = degrees; }


    /**
     * \brief Set the alpha angle value at the right ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaRvEndo(double degrees) { this->alpha_rv_endo_ = degrees; }


    /**
     * \brief Set the alpha angle value at the right ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaRvEpi(double degrees) { this->alpha_rv_epi_ = degrees; }


    /**
     * \brief Set the alpha angle value at the outflow tract left ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaOtLvEndo(double degrees) { this->alpha_ot_lv_endo_ = degrees; }

    /**
     * \brief Set the alpha angle value at the outflow tract left ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaOtLvEpi(double degrees) { this->alpha_ot_lv_epi_ = degrees; }


    /**
     * \brief Set the alpha angle value at the outflow tract right ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaOtRvEndo(double degrees) { this->alpha_ot_rv_endo_ = degrees; }

    /**
     * \brief Set the alpha angle value at the outflow tract right ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaOtRvEpi(double degrees) { this->alpha_ot_rv_epi_ = degrees; }


    /**
     * \brief Set the alpha angle value at the septum (interventricular septal surface).
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetAlphaSeptum(double degrees) { this->alpha_septum_ = degrees; }


    /**
     * \brief Set the beta angle value at the left ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaLvEndo(double degrees) { this->beta_lv_endo_ = degrees; }


    /**
     * \brief Set the beta angle value at the left ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaLvEpi(double degrees) { this->beta_lv_epi_ = degrees; }


    /**
     * \brief Set the beta angle value at the right ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaRvEndo(double degrees) { this->beta_rv_endo_ = degrees; }


    /**
     * \brief Set the beta angle value at the right ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaRvEpi(double degrees) { this->beta_rv_epi_ = degrees; }


    /**
     * \brief Set the beta angle value at the outflow tract left ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaOtLvEndo(double degrees) { this->beta_ot_lv_endo_ = degrees; }


    /**
     * \brief Set the beta angle value at the outflow tract left ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaOtLvEpi(double degrees) { this->beta_ot_lv_epi_ = degrees; }


    /**
     * \brief Set the beta angle value at the outflow tract right ventricle endocardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaOtRvEndo(double degrees) { this->beta_ot_rv_endo_ = degrees; }


    /**
     * \brief Set the beta angle value at the outflow tract right ventricle epicardium.
     * \param [in] degrees The angle value expressed in degrees.
     * \return void
     */
    inline void SetBetaOtRvEpi(double degrees) { this->beta_ot_rv_epi_ = degrees; }


    /**
     * \brief Get the alpha angle value at the left ventricle endocardium.
     * \return [double] The alpha angle value at the left ventricle endocardium.
     */
    double AlphaLvEndo() const { return this->alpha_lv_endo_; }


    /**
     * \brief Get the alpha angle value at the left ventricle epicardium.
     * \return [double] The alpha angle value at the left ventricle epicardium.
     */
    double AlphaLvEpi() const { return this->alpha_lv_epi_; }


    /**
     * \brief Get the alpha angle value at the right ventricle endocardium.
     * \return [double] The alpha angle value at the right ventricle endocardium.
     */
    double AlphaRvEndo() const { return this->alpha_rv_endo_; }


    /**
     * \brief Get the alpha angle value at the right ventricle epicardium.
     * \return [double] The alpha angle value at the right ventricle epicardium.
     */
    double AlphaRvEpi() const { return this->alpha_rv_epi_; }


    /**
     * \brief Get the alpha angle value at the outflow tract left ventricle endocardium.
     * \return [double] The alpha angle value at the outflow tract left ventricle endocardium.
     */
    double AlphaOtLvEndo() const { return this->alpha_ot_lv_endo_; }


    /**
     * \brief Get the alpha angle value at the outflow tract left ventricle epicardium.
     * \return [double] The alpha angle value at the outflow tract left ventricle epicardium.
     */
    double AlphaOtLvEpi() const { return this->alpha_ot_lv_epi_; }


    /**
     * \brief Get the alpha angle value at the outflow tract right ventricle endocardium.
     * \return [double] The alpha angle value at the outflow tract right ventricle endocardium.
     */
    double AlphaOtRvEndo() const { return this->alpha_ot_rv_endo_; }


    /**
     * \brief Get the alpha angle value at the outflow tract right ventricle epicardium.
     * \return [double] The alpha angle value at the outflow tract right ventricle epicardium.
     */
    double AlphaOtRvEpi() const { return this->alpha_ot_rv_epi_; }


    /**
     * \brief Get the alpha angle value at the septum (interventricular septal surface).
     * \return [double] The alpha angle value at the septum (interventricular septal surface).
     */
    double AlphaSeptum() const { return this->alpha_septum_; }


    /**
     * \brief Get the beta angle value at the left endocardium.
     * \return [double] The beta angle value at the left endocardium.
     */
    double BetaLvEndo() const { return this->beta_lv_endo_; }


    /**
     * \brief Get the beta angle value at the left epicardium.
     * \return [double] The beta angle value at the left epicardium.
     */
    double BetaLvEpi() const { return this->beta_lv_epi_; }


    /**
     * \brief Get the beta angle value at the right endocardium.
     * \return [double] The beta angle value at the right endocardium.
     */
    double BetaRvEndo() const { return this->beta_rv_endo_; }


    /**
     * \brief Get the beta angle value at the right epicardium.
     * \return [double] The beta angle value at the right epicardium.
     */
    double BetaRvEpi() const { return this->beta_rv_epi_; }


    /**
     * \brief Get the beta angle value at the outflow tract left ventricle endocardium.
     * \return [double] The beta angle value at the outflow tract left ventricle endocardium.
     */
    double BetaOtLvEndo() const { return this->beta_ot_lv_endo_; }


    /**
     * \brief Get the beta angle value at the outflow tract left ventricle epicardium.
     * \return [double] The beta angle value at the outflow tract left ventricle epicardium.
     */
    double BetaOtLvEpi() const { return this->beta_ot_lv_epi_; }


    /**
     * \brief Get the beta angle value at the outflow tract right ventricle endocardium.
     * \return [double] The beta angle value at the outflow tract right ventricle endocardium.
     */
    double BetaOtRvEndo() const { return this->beta_ot_rv_endo_; }


    /**
     * \brief Get the beta angle value at the outflow tract right ventricle epicardium.
     * \return [double] The beta angle value at the outflow tract right ventricle epicardium.
     */
    double BetaOtRvEpi() const { return this->beta_ot_rv_epi_; }

};

} // End of namespace ELECTRA

#endif //ELECTRA_FIBERS_VENTRI_FIBER_RULES_HPP_
