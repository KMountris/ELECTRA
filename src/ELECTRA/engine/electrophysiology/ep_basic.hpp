/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file ep_basic.hpp
   \brief EpBasic class header file.
   \author Konstantinos A. Mountris
   \date 13/03/2019
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_EP_BASIC_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_EP_BASIC_HPP_

#include "ELECTRA/engine/electrophysiology/ep_varying_params.hpp"
#include "ELECTRA/engine/utilities/measure_units.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"
#include "ELECTRA/engine/utilities/algorithm.hpp"

#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <stdexcept>
#include <exception>
#include <algorithm>
#include <sstream>


namespace ELECTRA {

/** \addtogroup Electrophysiology \{ */


/**
 * \enum CellType
 * \author Konstantinos A. Mountris
 * \brief Enumeration to declare the type of a cell.
 */
enum class CellType {
    ventricular,  /**< General ventricular cell type */
    endo,         /**< Endocardium ventricular cell type */
    mid,          /**< Mid-myocardium ventricular cell type */
    epi,          /**< Epicardium ventricular cell type */
    fibro,        /**< Fibroblast cell type */
    atrial,       /**< Generic atrial cell type */
    left_atrial,  /**< Left atrial cell type */
    right_atrial, /**< Right atrial cell type */
    purkinje      /**< Purkinje cell type */
};


/**
 * \enum EpModelType
 * \author Konstantinos A. Mountris
 * \brief Enumeration to declare the type of the electrophysiology model.
 */
enum class EpModelType {
    unknown,
    MacCannell,        /**< MacCannell '07 fibroblast cell electrophysiology model */
    Bueno,             /**< Bueno-Orovio '07 human ventricular cell electrophysiology model */
    OHara,             /**< O'Hara Rudy '11 human ventricular cell electrophysiology model */
    Gaur2021,          /**< Gaur '21 porcine ventricular cell electrophysiology model */
    Gong2020,          /**< Gong et al 2020 human ventricular cell electrophysiology model */
    Gong2020m,         /**< Gong et al 2020 human ventricular cell electrophysiology model with modified INa gates from TenTusscher2006 */
    PaciVentri,        /**< Paci '13 stem cell derived ventricular cell electrophysiology model */
    TenTusscher2006,   /**< TenTusscher et al 2006 human ventricular cell electrophysiology model */
    Courtemanche,      /**< Courtemanche '98 atrial cell electrophysiology model */
    GrandiAtri,        /**< Grandi '11 atrial cell electrophysiology model */
    Maleckar2009,      /**< Maleckar et al 2009 atrial cell electrophysiology model */
    Stewart            /**< Stewart '09 purkinje cell electrophysiology model */
};


/**
 * \class EpBasic
 * \author Konstantinos A. Mountris
 * \brief Generic abstract electrophysiology model.
 */
class EpBasic
{

protected:

    EpModelType model_type_ = EpModelType::unknown;  /**< The type of the electrophysiology model. Initialized to unknown */

    double dt_stable_ = 0.;            /**< The maximum stable time step of the electrophysiology model */

    double integ_time_ = 0.;           /**< The elapsed integration time during electrophysiology computation */

    std::vector<double> var_;          /**< The variables of the electrophysiology model */

    std::vector<double> prm_;          /**< The parameters of the electrophysiology model */

    std::vector<double> cur_;          /**< The ionic currents of the electrophysiology model */

    std::vector<double> block_coeff_;  /**< The ionic current block coefficients of the electrophysiology model */

    std::unordered_map<std::string, std::size_t> mapped_data_;  /**< Mapping for accessing ap data efficiently */

    /**
     * \brief Set the mapping for data access.
     */
    virtual void SetDataMapping() = 0;


public:

    /**
     * \brief The default constructor of the EpBasic class.
     */
    explicit EpBasic() {}


    /**
     * \brief The destructor of the EpBasic class.
     */
    virtual ~EpBasic() {}


    /**
     * \brief To be used by the inheriting electrophysiology model to initialize the variables and parameters.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type) = 0;


    /**
     * \brief To be used by the inheriting electrophysiology model to compute the temporal update.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current) = 0;


    virtual void UpdateTimeVaryingPrms(const EpVaryingParams &var_params, double dt);


    inline virtual void ResetIntegrationTime() { this->integ_time_ = 0.; }


    /**
     * \brief To be used by the inheriting electrophysiology model to set the membrane electrophysiology value.
     * \return [void]
     */
    virtual void SetV(double v) = 0;


    /**
     * \brief Set a specific variable of the electrophysiology model.
     * \param id
     * \param value
     * \return [void]
     */
    inline virtual void SetVar(std::size_t id, double value) { this->var_.at(id) = value; }


    /**
     * \brief Set a specific parameter of the electrophysiology model.
     * \param id
     * \param value
     * \return [void]
     */
    inline virtual void SetPrm(std::size_t id, double value) { this->prm_.at(id) = value; }


    /**
     * \brief Set a specific current of the electrophysiology model.
     * \param id
     * \param value
     * \return [void]
     */
    inline virtual void SetCurrent(std::size_t id, double value) { this->cur_.at(id) = value; }


    /**
     * \brief Set the inhibition coefficient for a specified current of the electrophysiology model. 
     * \param [in] id The index of the desired current. Use ap model enumeration for safe access. 
     * \param [in] value The inhibition coefficient percentage value expressed in decimal form.
     * \return [void]
     */
    inline virtual void SetBlockCoeff(std::size_t id, double value) { this->block_coeff_.at(id) = value; }


    /*!
     * \brief To be used by the inheriting electrophysiology model to print to std::string its variables and their values.
     * \return [std::string] The electrophysiology model's variables and their values.
    */
    virtual std::string PrintVariables() const = 0;


    /*!
     * \brief To be used by the inheriting electrophysiology model to print to std::string its parameters and their values.
     * \return [std::string] The electrophysiology model's parameters and their values.
    */
    virtual std::string PrintParameters() const = 0;


    /*!
     * \brief To be used by the inheriting electrophysiology model to print to std::string its currents and their values.
     * \return [std::string] The electrophysiology model's currents and their values.
    */
    virtual std::string PrintCurrents() const = 0;


    /*!
     * \brief To be used by the inheriting electrophysiology model to print to std::string its currents' inhibition coefficients and their values.
     * \return [std::string] The electrophysiology model's currents' inhibition coefficients and their values.
    */
    virtual std::string PrintBlockCoeffs() const = 0;


    /**
     * \brief To be used by the inheriting electrophysiology model to get the membrane potential value.
     * \return [double] The membrane electrophysiology value.
     */
    virtual double V() const = 0;


    /**
     * \brief To be used by the inheriting electrophysiology model to get the membrane potential's time derivative.
     * \return [double] The membrane potential's time derivative.
     */
    virtual double dVdt() const = 0;


    /**
     * \brief To be used by the inheriting electrophysiology model to get the total ionic current value.
     * \return [double] The total ionic current value.
     */
    virtual double Iion() const = 0;


    /**
     * \brief
     * \return int
     */
    inline int VarNum() const { return static_cast<int>(this->var_.size()); }


    /**
     * \brief
     * \return int
     */
    inline int PrmNum() const { return static_cast<int>(this->prm_.size()); }


    /**
     * \brief
     * \return int
     */
    inline int CurrentNum() const { return static_cast<int>(this->cur_.size()); }


    /**
     * \brief Get the electrophysiology model's type.
     * \return [ELECTRA::EpModelType] The electrophysiology model's type.
     */
    inline EpModelType ModelType() const { return this->model_type_; }


    /**
     * \brief Get the value of a variable of the electrophysiology model.
     *
     * Fast access, no range check.
     *
     * \param [in] id The index of the desired variable. Use ap model enumeration for safe access. 
     * \return [double] The variable value.
     */
    inline double Var(std::size_t id) { return this->var_[id]; }


    /**
     * \brief Get the value of a variable of the electrophysiology model.
     *
     * Slower access, with range check.
     *
     * \param [in] id The index of the desired variable. Use ap model enumeration for safe access. 
     * \return [double] The variable value.
     */
    inline double VarAt(std::size_t id) { return this->var_.at(id); }


    /**
     * \brief Get the value of a parameter of the electrophysiology model.
     *
     * Fast access, no range check.
     *
     * \param [in] id The index of the desired parameter. Use ap model enumeration for safe access. 
     * \return [double] The parameter value.
     */
    inline double Prm(std::size_t id) { return this->prm_[id]; }


    /**
     * \brief Get the value of a parameter of the electrophysiology model.
     *
     * Slower access, with range check.
     *
     * \param [in] id The index of the desired parameter. Use ap model enumeration for safe access. 
     * \return [double] The parameter value.
     */
    inline double PrmAt(std::size_t id) { return this->prm_.at(id); }


    /**
     * \brief Get the value of a current of the electrophysiology model.
     *
     * Fast access, no range check.
     *
     * \param [in] id The index of the desired current. Use ap model enumeration for safe access. 
     * \return [double] The current value.
     */
    inline double Current(std::size_t id) { return this->cur_[id]; }


    /**
     * \brief Get the value of a current of the electrophysiology model.
     *
     * Slower access, with range check.
     *
     * \param [in] id The index of the desired current. Use ap model enumeration for safe access. 
     * \return [double] The current value.
     */
    inline double CurrentAt(std::size_t id) { return this->cur_.at(id); }


    /**
     * \brief Get the block coefficient for a specified current of the electrophysiology model.
     *
     * Fast access, no range check.
     *
     * \param [in] id The index of the desired current. Use ap model enumeration for safe access. 
     * \return [double] The block coefficient percentage value expressed in decimal form for the required current with index cur.
     */
    inline double BlockCoeff(std::size_t id) { return this->block_coeff_[id]; }


    /**
     * \brief Get the block coefficient for a specified current of the electrophysiology model.
     *
     * Slower access, with range check.
     *
     * \param [in] id The index of the desired current. Use ap model enumeration for safe access. 
     * \return [double] The block coefficient percentage value expressed in decimal form for the required current with index cur.
     */
    inline double BlockCoeffAt(std::size_t id) const { return this->block_coeff_.at(id); }


    /**
     * \brief Get a mapped data entry of the electrophysiology model.
     * \param key The name of the required data entry.
     * \return [std::size_t] The mapped index to access the data container.
     */
    inline std::size_t MappedDataEntry(const std::string &key) const { return this->mapped_data_.at(key); }


    /**
     * \brief Check if the electrophysiology model has the specified data by the given key.
     * \param [in] key The key of the data to be checked if exist in the electrophysiology model.
     * \return [true] The data with the specified key exist in the model.
     * \return [false] The data with the specified key does not exist in the model.
     */
    inline bool HasDataEntry(const std::string &key) const { if (this->mapped_data_.find(key) != this->mapped_data_.end()) { return true; } return false; }

};


/*! \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_EP_BASIC_HPP_