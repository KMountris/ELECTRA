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
   \file ventri_tags.hpp
   \brief VentriTags class header file.
   \author Konstantinos A. Mountris
   \date 30/04/2021
*/

#ifndef ELECTRA_FIBERS_ATRI_TAGS_HPP_
#define ELECTRA_FIBERS_ATRI_TAGS_HPP_

#include <string>

namespace ELECTRA
{

/** \addtogroup Fibers \{ */

/**
 * \class AtriTags
 * \brief Class collecting name tags for atrial partitions that are required during cardiac fiber determination with Ldrbm (Laplace-Dirichlet rule based method).
 */
class AtriTags
{
private:

    std::string la_endo_wall_tag_;              // The name tag of the left atrium endocardium wall node set

    std::string la_appendage_tag_;              // The name tag of the left atrium appendage node set

    std::string ra_endo_wall_tag_;              // The name tag of the right atrium endocardium wall node set

    std::string ra_appendage_tag_;              // The name tag of the right atrium appendage node set

    std::string ra_top_endo_tag_;               // The name tag of the endocardium part of top upper region of the inferior and superior caval vein node set

    std::string ra_top_epi_tag_;                // The name tag of the epicardium part of top upper region of the inferior and superior caval vein node set

    std::string epi_wall_tag_;                  // The name tag of the atrium epicardium wall node set

    std::string mitral_valve_tag_;              // The name tag of the left atrium mitral valve ring node set

    std::string left_pulmon_vein_tag_;          // The name tag of the left atrium left pulmonary vein ring node set

    std::string right_pulmon_vein_tag_;         // The name tag of the left atrium right pulmonary vein ring node set

    std::string superior_caval_vein_tag_;       // The name tag of the right atrium superior caval vein ring node set

    std::string inferior_caval_vein_tag_;       // The name tag of the right atrium inferior caval vein ring node set

    std::string tricuspid_valve_septum_tag_;    // The name tag of the septal part of the right atrium tricuspid valve ring node set

    std::string tricuspid_valve_free_tag_;      // The name tag of the free part of the right atrium tricuspid valve ring node set


public:

    /**
     * \brief Default constructor of atrial name tags.
     */
    inline AtriTags() :
        la_endo_wall_tag_(""), la_appendage_tag_(""), ra_endo_wall_tag_(""), ra_appendage_tag_(""), ra_top_endo_tag_(""),
        ra_top_epi_tag_(""), epi_wall_tag_(""), mitral_valve_tag_(""), left_pulmon_vein_tag_(""), right_pulmon_vein_tag_(""),
        superior_caval_vein_tag_(""), inferior_caval_vein_tag_(""), tricuspid_valve_septum_tag_(""), tricuspid_valve_free_tag_("") {}


    /**
     * \brief Default destructor of atrial name tags.
     */
    inline virtual ~AtriTags() {}


    /**
     * \brief Set the name tag of the left atrium endocardium wall node set.
     * \param [in] la_endo_wall_tag The name tag of the left atrium endocardium wall node set.
     * \return [void]
     */
    inline void SetLaEndoWallTag(const std::string &la_endo_wall_tag) { this->la_endo_wall_tag_ = la_endo_wall_tag; }


    /**
     * \brief Set the name tag of the left atrium appendage node set.
     * \param [in] la_appendage_tag The name tag of the left atrium appendage node set.
     * \return [void]
     */
    inline void SetLaAppendageTag(const std::string &la_appendage_tag) { this->la_appendage_tag_ = la_appendage_tag; }


    /**
     * \brief Set the name tag of the right atrium endocardium wall node set.
     * \param [in] ra_endo_wall_tag The name tag of the right atrium endocardium wall node set.
     * \return [void]
     */
    inline void SetRaEndoWallTag(const std::string &ra_endo_wall_tag) { this->ra_endo_wall_tag_ = ra_endo_wall_tag; }


    /**
     * \brief Set the name tag of the right atrium appendage node set.
     * \param [in] ra_appendage_tag The name tag of the right atrium appendage node set.
     * \return [void]
     */
    inline void SetRaAppendageTag(const std::string &ra_appendage_tag) { this->ra_appendage_tag_ = ra_appendage_tag; }


    /**
     * \brief Set the name tag of the endocardium part of top upper region of the inferior and superior caval vein node set.
     * \param [in] ra_top_endo_tag The name tag of the endocardium part of top upper region of the inferior and superior caval vein node set.
     * \return [void]
     */
    inline void SetRaTopEndoTag(const std::string &ra_top_endo_tag) { this->ra_top_endo_tag_ = ra_top_endo_tag; }


    /**
     * \brief Set the name tag of the epicardium part of top upper region of the inferior and superior caval vein node set.
     * \param [in] ra_top_epi_tag The name tag of the epicardium part of top upper region of the inferior and superior caval vein node set.
     * \return [void]
     */
    inline void SetRaTopEpiTag(const std::string &ra_top_epi_tag) { this->ra_top_epi_tag_ = ra_top_epi_tag; }


    /**
     * \brief Set the name tag of the atrium epicardium wall node set.
     * \param [in] epi_wall_tag The name tag of the atrium epicardium wall node set.
     * \return [void]
     */
    inline void SetEpiWallTag(const std::string &epi_wall_tag) { this->epi_wall_tag_ = epi_wall_tag; }


    /**
     * \brief Set the name tag of the left atrium mitral valve ring node set.
     * \param [in] mitral_valve_tag The name tag of the left atrium mitral valve ring node set.
     * \return [void]
     */
    inline void SetMitralValveTag(const std::string &mitral_valve_tag) { this->mitral_valve_tag_ = mitral_valve_tag; }


    /**
     * \brief Set the name tag of the left atrium left pulmonary vein ring node set.
     * \param [in] left_pulmon_vein_tag The name tag of the left atrium left pulmonary vein ring node set.
     * \return [void]
     */
    inline void SetLeftPulmonVeinTag(const std::string &left_pulmon_vein_tag) { this->left_pulmon_vein_tag_ = left_pulmon_vein_tag; }


    /**
     * \brief Set the name tag of the left atrium right pulmonary vein ring node set.
     * \param [in] right_pulmon_vein_tag The name tag of the left atrium right pulmonary vein ring node set.
     * \return [void]
     */
    inline void SetRightPulmonVeinTag(const std::string &right_pulmon_vein_tag) { this->right_pulmon_vein_tag_ = right_pulmon_vein_tag; }


    /**
     * \brief Set the name tag of the right atrium superior caval vein ring node set.
     * \param [in] superior_caval_vein_tag The name tag of the right atrium superior caval vein ring node set.
     * \return [void]
     */
    inline void SetSuperiorCavalVeinTag(const std::string &superior_caval_vein_tag) { this->superior_caval_vein_tag_ = superior_caval_vein_tag; }


    /**
     * \brief Set the name tag of the right atrium inferior caval vein ring node set.
     * \param [in] inferior_caval_vein_tag The name tag of the right atrium inferior caval vein ring node set.
     * \return [void]
     */
    inline void SetInferiorCavalVeinTag(const std::string &inferior_caval_vein_tag) { this->inferior_caval_vein_tag_ = inferior_caval_vein_tag; }


    /**
     * \brief Set the name tag of the septum part of the right tricuspid valve ring node set.
     * \param [in] tricuspid_valve_septum_tag The name tag of the septum part of the right tricuspid valve ring node set.
     * \return [void]
     */
    inline void SetTricuspidValveSeptumTag(const std::string &tricuspid_valve_septum_tag) { this->tricuspid_valve_septum_tag_ = tricuspid_valve_septum_tag; }


    /**
     * \brief Set the name tag of the free part of the right atrium tricuspid valve ring node set.
     * \param [in] tricuspid_valve_free_tag The name tag of the free part of the right atrium tricuspid valve ring node set.
     * \return [void]
     */
    inline void SetTricuspidValveFreeTag(const std::string &tricuspid_valve_free_tag) { this->tricuspid_valve_free_tag_ = tricuspid_valve_free_tag; }


    /**
     * \brief Get the name tag of the left atrium endocardium wall node set.
     * \return [const std::string&] The name tag of the left atrium endocardium wall node set.
     */
    const std::string & LaEndoWallTag() const { return this->la_endo_wall_tag_; }


    /**
     * \brief Get the name tag of the left atrium appendage node set.
     * \return [const std::string&] The name tag of the left atrium appendage node set.
     */
    const std::string & LaAppendageTag() const { return this->la_appendage_tag_; }


    /**
     * \brief Get the name tag of the right atrium endocardium wall node set.
     * \return [const std::string&] The name tag of the right atrium endocardium wall node set.
     */
    const std::string & RaEndoWallTag() const { return this->ra_endo_wall_tag_; }


    /**
     * \brief Get the name tag of the right atrium appendage node set.
     * \return [const std::string&] The name tag of the right atrium appendage node set.
     */
    const std::string & RaAppendageTag() const { return this->ra_appendage_tag_; }


    /**
     * \brief Get the name tag of the endocardium part of top upper region of the inferior and superior caval vein node set.
     * \return [const std::string&] The name tag of the endocardium part of top upper region of the inferior and superior caval vein node set.
     */
    const std::string & RaTopEndoTag() const { return this->ra_top_endo_tag_; }


    /**
     * \brief Get the name tag of the epicardium part of top upper region of the inferior and superior caval vein node set.
     * \return [const std::string&] The name tag of the epicardium part of top upper region of the inferior and superior caval vein node set.
     */
    const std::string & RaTopEpiTag() const { return this->ra_top_epi_tag_; }


    /**
     * \brief Get the name tag of the atrium epicardium wall node set.
     * \return [const std::string&] The name tag of the atrium epicardium wall node set.
     */
    const std::string & EpiWallTag() const { return this->epi_wall_tag_; }


    /**
     * \brief Get the name tag of the left atrium mitral valve ring node set.
     * \return [const std::string&] The name tag of the left atrium mitral valve ring node set.
     */
    const std::string & MitralValveTag() const { return this->mitral_valve_tag_; }


    /**
     * \brief Get the name tag of the left atrium left pulmonary vein ring node set.
     * \return [const std::string&] The name tag of the left atrium left pulmonary vein ring node set.
     */
    const std::string & LeftPulmonVeinTag() const { return this->left_pulmon_vein_tag_; }


    /**
     * \brief Get the name tag of the left atrium right pulmonary vein ring node set.
     * \return [const std::string&] The name tag of the left atrium right pulmonary vein ring node set.
     */
    const std::string & RightPulmonVeinTag() const { return this->right_pulmon_vein_tag_; }


    /**
     * \brief Get the name tag of the right atrium superior caval vein ring node set.
     * \return [const std::string&] The name tag of the right atrium superior caval vein ring node set.
     */
    const std::string & SuperiorCavalVeinTag() const { return this->superior_caval_vein_tag_; }


    /**
     * \brief Get the name tag of the right atrium inferior caval vein ring node set.
     * \return [const std::string&] The name tag of the right atrium inferior caval vein ring node set.
     */
    const std::string & InferiorCavalVeinTag() const { return this->inferior_caval_vein_tag_; }


    /**
     * \brief Get the name tag of the septum part of the right tricuspid valve ring node set.
     * \return [const std::string&] The name tag of the septum part of the right tricuspid valve ring node set.
     */
    const std::string & TricuspidValveSeptumTag() const { return this->tricuspid_valve_septum_tag_; }


    /**
     * \brief Get the name tag of the free part of the right tricuspid valve ring node set.
     * \return [const std::string&] The name tag of the free part of the right tricuspid valve ring node set.
     */
    const std::string & TricuspidValveFreeTag() const { return this->tricuspid_valve_free_tag_; }


};

} // End of namespace ELECTRA

#endif //ELECTRA_FIBERS_ATRI_TAGS_HPP_
