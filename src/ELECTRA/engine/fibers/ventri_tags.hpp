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
   \date 18/03/2021
*/

#ifndef ELECTRA_FIBERS_VENTRI_TAGS_HPP_
#define ELECTRA_FIBERS_VENTRI_TAGS_HPP_

#include <string>

namespace ELECTRA
{

/** \addtogroup Fibers \{ */

/**
 * \class VentriTags
 * \brief Class collecting name tags for ventricular partitions that are required during cardiac fiber determination with Ldrbm (Laplace-Dirichlet rule based method).
 */
class VentriTags
{
private:

    std::string lv_endo_wall_tag_;      // The name tag of the left ventricle's endocardium wall node set

    std::string lv_apex_tag_;           // The name tag of the left ventricle's apex node set

    std::string rv_endo_wall_tag_;      // The name tag of the right ventricle's endocardium wall node set

    std::string rv_apex_tag_;           // The name tag of the right ventricle's apex node set

    std::string epi_wall_tag_;          // The name tag of the ventricles' epicardium wall node set

    std::string mitral_valve_tag_;      // The name tag of the left ventricle's mitral valve ring node set

    std::string aortic_valve_tag_;      // The name tag of the left ventricle's aortic valve ring node set

    std::string tricuspid_valve_tag_;   // The name tag of the right ventricle's tricuspid valve ring node set

    std::string pulmonary_valve_tag_;   // The name tag of the left ventricle's pulmonary valve ring node set


public:

    /**
     * \brief Default constructor of ventricular name tags.
     */
    inline VentriTags() :
        lv_endo_wall_tag_(""), lv_apex_tag_(""), rv_endo_wall_tag_(""), rv_apex_tag_(""), epi_wall_tag_(""),
        mitral_valve_tag_(""), aortic_valve_tag_(""), tricuspid_valve_tag_(""), pulmonary_valve_tag_("") {}


    /**
     * \brief Default destructor of ventricular name tags.
     */
    inline virtual ~VentriTags() {}


    /**
     * \brief Set the name tag of the left ventricle's endocardium wall node set.
     * \param [in] lv_endo_wall_tag The name tag of the left ventricle's endocardium wall node set.
     * \return [void]
     */
    inline void SetLvEndoWallTag(const std::string &lv_endo_wall_tag) { this->lv_endo_wall_tag_ = lv_endo_wall_tag; }


    /**
     * \brief Set the name tag of the left ventricles' apex node set.
     * \param [in] lv_apex_tag The name tag of the left ventricles' apex node set.
     * \return [void]
     */
    inline void SetLvApexTag(const std::string &lv_apex_tag) { this->lv_apex_tag_ = lv_apex_tag; }


    /**
     * \brief Set the name tag of the right ventricle's endocardium wall node set.
     * \param [in] rv_endo_wall_tag The name tag of the right ventricle's endocardium wall node set.
     * \return [void]
     */
    inline void SetRvEndoWallTag(const std::string &rv_endo_wall_tag) { this->rv_endo_wall_tag_ = rv_endo_wall_tag; }


    /**
     * \brief Set the name tag of the right ventricles' apex node set.
     * \param [in] rv_apex_tag The name tag of the right ventricles' apex node set.
     * \return [void]
     */
    inline void SetRvApexTag(const std::string &rv_apex_tag) { this->rv_apex_tag_ = rv_apex_tag; }



    /**
     * \brief Set the name tag of the ventricles' epicardium wall node set.
     * \param [in] epi_wall_tag The name tag of the ventricles' epicardium wall node set.
     * \return [void]
     */
    inline void SetEpiWallTag(const std::string &epi_wall_tag) { this->epi_wall_tag_ = epi_wall_tag; }


    /**
     * \brief Set the name tag of the left ventricles' mitral valve ring node set.
     * \param [in] mitral_valve_tag The name tag of the left ventricles' mitral valve ring node set.
     * \return [void]
     */
    inline void SetMitralValveTag(const std::string &mitral_valve_tag) { this->mitral_valve_tag_ = mitral_valve_tag; }


    /**
     * \brief Set the name tag of the left ventricles' aortic valve ring node set.
     * \param [in] aortic_valve_tag The name tag of the left ventricles' aortic valve ring node set.
     * \return [void]
     */
    inline void SetAorticValveTag(const std::string &aortic_valve_tag) { this->aortic_valve_tag_ = aortic_valve_tag; }


    /**
     * \brief Set the name tag of the right ventricles' tricuspid valve ring node set.
     * \param [in] tricuspid_valve_tag The name tag of the right ventricles' tricuspid valve ring node set.
     * \return [void]
     */
    inline void SetTricuspidValveTag(const std::string &tricuspid_valve_tag) { this->tricuspid_valve_tag_ = tricuspid_valve_tag; }


    /**
     * \brief Set the name tag of the right ventricles' pulmonary valve ring node set.
     * \param [in] pulmonary_valve_tag_ The name tag of the right ventricles' pulmonary valve ring node set.
     * \return [void]
     */
    inline void SetPulmonaryValveTag(const std::string &pulmonary_valve_tag) { this->pulmonary_valve_tag_ = pulmonary_valve_tag; }


    /**
     * \brief Get the name tag of the left ventricle's endocardium wall node set.
     * \return [const std::string&] The name tag of the left ventricle's endocardium wall node set.
     */
    const std::string & LvEndoWallTag() const { return this->lv_endo_wall_tag_; }


    /**
     * \brief Get the name tag of the left ventricles' apex node set.
     * \return [const std::string&] The name tag of the left ventricles' apex node set.
     */
    const std::string & LvApexTag() const { return this->lv_apex_tag_; }


    /**
     * \brief Get the name tag of the right ventricle's endocardium wall node set.
     * \return [const std::string&] The name tag of the right ventricle's endocardium wall node set.
     */
    const std::string & RvEndoWallTag() const { return this->rv_endo_wall_tag_; }


    /**
     * \brief Get the name tag of the right ventricles' apex node set.
     * \return [const std::string&] The name tag of the right ventricles' apex node set.
     */
    const std::string & RvApexTag() const { return this->rv_apex_tag_; }


    /**
     * \brief Get the name tag of the ventricles' epicardium wall node set.
     * \return [const std::string&] The name tag of the ventricles' epicardium wall node set.
     */
    const std::string & EpiWallTag() const { return this->epi_wall_tag_; }


    /**
     * \brief Get the name tag of the left ventricles' mitral valve ring node set.
     * \return [const std::string&] The name tag of the left ventricles' mitral valve ring node set.
     */
    const std::string & MitralValveTag() const { return this->mitral_valve_tag_; }


    /**
     * \brief Get the name tag of the left ventricles' aortic valve ring node set.
     * \return [const std::string&] The name tag of the left ventricles' aortic valve ring node set.
     */
    const std::string & AorticValveTag() const { return this->aortic_valve_tag_; }


    /**
     * \brief Get the name tag of the right ventricles' tricuspid valve ring node set.
     * \return [double] The name tag of the right ventricles' tricuspid valve ring node set.
     */
    const std::string & TricuspidValveTag() const { return this->tricuspid_valve_tag_; }


    /**
     * \brief Get the name tag of the right ventricles' pulmonary valve ring node set.
     * \return [const std::string&] The name tag of the right ventricles' pulmonary valve ring node set.
     */
    const std::string & PulmonaryValveTag() const { return this->pulmonary_valve_tag_; }

};

} // End of namespace ELECTRA

#endif //ELECTRA_FIBERS_VENTRI_TAGS_HPP_