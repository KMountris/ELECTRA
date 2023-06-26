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
   \file atri_fiber_rules.hpp
   \brief AtriFiberRules class header file.
   \author Konstantinos A. Mountris
   \date 03/05/2021
*/

#ifndef ELECTRA_FIBERS_ATRI_FIBER_RULES_HPP_
#define ELECTRA_FIBERS_ATRI_FIBER_RULES_HPP_

namespace ELECTRA
{

/** \addtogroup Fibers \{ */

/**
 * \class AtriFiberRules
 * \brief Class collecting laplacian rules to distinguish atrial bundles for cardiac fiber determination with Ldrbm (Laplace-Dirichlet rule based method).
 */
class AtriFiberRules
{
private:

    double tricuspid_valve_rule_;          // The tricuspid valve rule at the right atrium

    double inf_caval_vein_rule_;            // The inferior caval vein rule at the right atrium

    double sup_caval_vein_rule_;            // The superior caval vein rule at the right atrium

    double crista_terminal_plus_rule_;      // The positive part of crista terminalis rule at the right atrium

    double crista_terminal_minus_rule_;     // The negative part of crista terminalis rule at the right atrium

    double inter_caval_bundle_rule_;        // The inter caval bundle rule at the right atrium

    double right_septum_wall_rule_;         // The septum wall rule at the right atrium

    double right_lateral_wall_rule_;        // The lateral wall rule at the right atrium

    double mitral_valve_rule_;              // The mitral valve rule at the left atrium

    double left_pulmon_vein_rule_;          // The left pulmonary vein rule at the left atrium

    double right_pulmon_vein_rule_;         // The right pulmonary vein rule at the left atrium


public:

    /**
     * \brief Default constructor of atria fiber rules.
     */
    inline AtriFiberRules() :
        tricuspid_valve_rule_(0.9), inf_caval_vein_rule_(0.85), sup_caval_vein_rule_(0.3),
        crista_terminal_plus_rule_(-0.55), crista_terminal_minus_rule_(-0.6),
        inter_caval_bundle_rule_(-0.25), right_septum_wall_rule_(-0.1), right_lateral_wall_rule_(0.6),
        mitral_valve_rule_(0.85), left_pulmon_vein_rule_(0.85), right_pulmon_vein_rule_(0.2) {}

    /**
     * \brief Default destructor of atria fiber rules.
     */
    inline virtual ~AtriFiberRules() {}


    /**
     * \brief Set the value for the tricuspid valve rule at the right atrium.
     * \param [in] tricuspid_valve_rule The value for the tricuspid valve rule at the right atrium.
     * \return void
     */
    inline void SetTricuspidValveRule(double tricuspid_valve_rule) { this->tricuspid_valve_rule_ = tricuspid_valve_rule; }


    /**
     * \brief Set the value for the inferior caval vein rule at the right atrium.
     * \param [in] inf_caval_vein_rule The value for the inferior caval vein rule at the right atrium.
     * \return void
     */
    inline void SetInfCavalVeinRule(double inf_caval_vein_rule) { this->inf_caval_vein_rule_ = inf_caval_vein_rule; }


    /**
     * \brief Set the value for the superior caval vein rule at the right atrium.
     * \param [in] sup_caval_vein_rule The value for the superior caval vein rule at the right atrium.
     * \return void
     */
    inline void SetSupCavalVeinRule(double sup_caval_vein_rule) { this->sup_caval_vein_rule_ = sup_caval_vein_rule; }


    /**
     * \brief Set the value for the positive part of the crista terminalis rule at the right atrium.
     * \param [in] crista_terminal_plus_rule The value for the positive part of the crista terminalis rule at the right atrium.
     * \return void
     */
    inline void SetCristaTerminalPlusRule(double crista_terminal_plus_rule) { this->crista_terminal_plus_rule_ = crista_terminal_plus_rule; }


    /**
     * \brief Set the value for the negative part of the crista terminalis rule at the right atrium.
     * \param [in] crista_terminal_minus_rule The value for the negative part of the crista terminalis rule at the right atrium.
     * \return void
     */
    inline void SetCristaTerminalMinusRule(double crista_terminal_minus_rule) { this->crista_terminal_minus_rule_ = crista_terminal_minus_rule; }


    /**
     * \brief Set the value for the inter caval bundle rule at the right atrium.
     * \param [in] inter_caval_bundle_rule The value for the inter caval bundle rule at the right atrium.
     * \return void
     */
    inline void SetInterCavalBundleRule(double inter_caval_bundle_rule) { this->inter_caval_bundle_rule_ = inter_caval_bundle_rule; }


    /**
     * \brief Set the value for the septum wall rule at the right atrium.
     * \param [in] right_septum_wall_rule The value for the septum wall rule at the right atrium.
     * \return void
     */
    inline void SetRightSeptumWallRule(double right_septum_wall_rule) { this->right_septum_wall_rule_ = right_septum_wall_rule; }


    /**
     * \brief Set the value for the lateral wall rule at the right atrium.
     * \param [in] right_lateral_wall_rule The value for the lateral wall rule at the right atrium.
     * \return void
     */
    inline void SetRightLateralWallRule(double right_lateral_wall_rule) { this->right_lateral_wall_rule_ = right_lateral_wall_rule; }


    /**
     * \brief Set the value for the mitral valve rule at the left atrium.
     * \param [in] mitral_valve_rule The value for the mitral valve rule at the left atrium.
     * \return void
     */
    inline void SetMitralValveRule(double mitral_valve_rule) { this->mitral_valve_rule_ = mitral_valve_rule; }


    /**
     * \brief Set the value for the left pulmonary vein rule at the left atrium.
     * \param [in] left_pulmon_vein_rule The value for the left pulmonary vein rule at the left atrium.
     * \return void
     */
    inline void SetLeftPulmonVeinRule(double left_pulmon_vein_rule) { this->left_pulmon_vein_rule_ = left_pulmon_vein_rule; }


    /**
     * \brief Set the value for the right pulmonary vein rule at the left atrium.
     * \param [in] right_pulmon_vein_rule The value for the right pulmonary vein rule at the left atrium.
     * \return void
     */
    inline void SetRightPulmonVeinRule(double right_pulmon_vein_rule) { this->right_pulmon_vein_rule_ = right_pulmon_vein_rule; }


    /**
     * \brief Get the value for the tricuspid valve rule at the right atrium.
     * \return [double] The value for the tricuspid valve rule at the right atrium.
     */
    double TricuspidValveRule() const { return this->tricuspid_valve_rule_; }


    /**
     * \brief Get the value for the inferior caval vein rule at the right atrium.
     * \return [double] The value for the inferior caval vein rule at the right atrium.
     */
    double InfCavalVeinRule() const { return this->inf_caval_vein_rule_; }


    /**
     * \brief Get the value for the superior caval vein rule at the right atrium.
     * \return [double] The value for the superior caval vein rule at the right atrium.
     */
    double SupCavalVeinRule() const { return this->sup_caval_vein_rule_; }


    /**
     * \brief Get the value for the positive part of the crista terminalis rule at the right atrium.
     * \return [double] The value for the positive part of the crista terminalis rule at the right atrium.
     */
    double CristaTerminalPlusRule() const { return this->crista_terminal_plus_rule_; }


    /**
     * \brief Get the value for the negative part of the crista terminalis rule at the right atrium.
     * \return [double] The value for the negative part of the crista terminalis rule at the right atrium.
     */
    double CristaTerminalMinusRule() const { return this->crista_terminal_minus_rule_; }


    /**
     * \brief Get the value for the inter caval bundle rule at the right atrium.
     * \return [double] The value for the inter caval bundle rule at the right atrium.
     */
    double InterCavalBundleRule() const { return this->inter_caval_bundle_rule_; }


    /**
     * \brief Get the value for the septum wall rule at the right atrium.
     * \return [double] The value for the septum wall rule at the right atrium.
     */
    double RightSeptumWallRule() const { return this->right_septum_wall_rule_; }


    /**
     * \brief Get the value for the lateral wall rule at the right atrium.
     * \return [double] The value for the lateral wall rule at the right atrium.
     */
    double RightLateralWallRule() const { return this->right_lateral_wall_rule_; }


    /**
     * \brief Get the value for the mitral valve rule at the left atrium.
     * \return [double] The value for the mitral valve rule at the left atrium.
     */
    double MitralValveRule() const { return this->mitral_valve_rule_; }


    /**
     * \brief Get the value for the left pulmonary vein rule at the left atrium.
     * \return [double] The value for the left pulmonary vein rule at the left atrium.
     */
    double LeftPulmonVeinRule() const { return this->left_pulmon_vein_rule_; }


    /**
     * \brief Get the value for the right pulmonary vein rule at the left atrium.
     * \return [double] The value for the right pulmonary vein rule at the left atrium.
     */
    double RightPulmonVeinRule() const { return this->right_pulmon_vein_rule_; }

};

} // End of namespace ELECTRA

#endif //ELECTRA_FIBERS_ATRI_FIBER_RULES_HPP_
