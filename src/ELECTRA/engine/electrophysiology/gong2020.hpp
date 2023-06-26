/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

/**
   \file gong2020.hpp
   \brief Gong2020 class header file.
   \author Konstantinos A. Mountris
   \date 22/12/2020
*/

#ifndef ELECTRA_ELECTROPHYSIOLOGY_MODELS_GONG2020_HPP_
#define ELECTRA_ELECTROPHYSIOLOGY_MODELS_GONG2020_HPP_

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
 * \namespace ELECTRA::Gng20Var
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access Gong 2020 human ventricular action potential model variables. 
 */
namespace Gng20Var {
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
       CaMKt,
       hPf,   // INaF PKA gates and 2 more for Both-P
       hPs,
       jP,
       hBPf,
       hBPs,
       jBP,
       dP,   // ICaL PKA gates and 2 more for Both-P
       fPf,
       fPs,
       fcaPf,
       fcaPs,
       fBPf,
       fcaBPf,
       xs1P,  // IKs PKA gates
       xs2P,
       JrelP,
       JrelBP,
       beta_cav_Gs_aGTP,
       beta_eca_Gs_aGTP,
       beta_cyt_Gs_aGTP,
       beta_cav_Gs_bg,
       beta_eca_Gs_bg,
       beta_cyt_Gs_bg,
       beta_cav_Gs_aGDP,
       beta_eca_Gs_aGDP,
       beta_cyt_Gs_aGDP,
       cAMP_cav,
       cAMP_eca,
       cAMP_cyt,
       beta_cav_Rb1_pka_tot,
       beta_eca_Rb1_pka_tot,
       beta_cyt_Rb1_pka_tot,
       beta_cav_Rb1_grk_tot,
       beta_eca_Rb1_grk_tot,
       beta_cyt_Rb1_grk_tot,
       pka_cav_ARC,
       pka_cav_A2RC,
       pka_cav_A2R,
       pka_cav_C,
       pka_cav_PKIC,
       pka_eca_ARC,
       pka_eca_A2RC,
       pka_eca_A2R,
       pka_eca_C,
       pka_eca_PKIC,
       pka_cyt_ARC,
       pka_cyt_A2RC,
       pka_cyt_A2R,
       pka_cyt_C,
       pka_cyt_PKIC,
       PDE3_P_cav,
       PDE3_P_cyt,
       PDE4_P_cav,
       PDE4_P_eca,
       PDE4_P_cyt,
       inhib1_p,
       ICaLp,
       IKsp,
       iup_f_plb,
       f_tni,
       ina_f_ina,
       f_inak,
       RyRp,
       f_ikur,
       beta_cav_Rb2_pka_tot,
       beta_cav_Rb2_grk_tot,
       beta_cav_Gi_aGTP,
       beta_cav_Gi_bg,
       beta_cav_Gi_aGDP,
       beta_eca_Rb2_pka_tot,
       beta_eca_Rb2_grk_tot,
       beta_eca_Gi_aGTP,
       beta_eca_Gi_bg,
       beta_eca_Gi_aGDP
      };
} // End of namespace Gng20Var


/**
 * \namespace ELECTRA::Gng20Prm
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access O'Hara Rudy '11 human ventricular action potential model parameters. 
 */
namespace Gng20Prm {
enum {  celltype,    
        signaling,    
        electrophys,    
        cipa,    
        GKS_mul,    
        VP_d,    
        VP_f,    
        GCal_mul,    
        CiPA_scale_IKr,    
        CiPA_scale_IKs,    
        CiPA_scale_IK1,    
        CiPA_scale_ICaL,   
        CiPA_scale_INaL,   
        nao,    
        cao,    
        ko,     
        R,    
        T,    
        F,    
        L,    
        rad,    
        vcell,   
        Ageo,    
        Acap,    
        vmyo,    
        vnsr,    
        vjsr,    
        vss,     
        aCaMK,  
        bCaMK,  
        CaMKo,  
        PKNa,   
        GNa,    
        GNaL,   
        Gto,    
        GKr,    
        GKs,    
        GK1,    
        Gncx,    
        GKb,     
        GpCa,    
        PCa,     
        Pnak,    
        PNab,    
        PCab,    
        KmCaMK,   
        KmCaM,    
        SERCA_total,  
        RyR_total,    
        Leak_total,    
        Trans_total,    
        iso_L,    
        ibmx,    
        phys_T,    
        phys_R,    
        length,    
        cell_pi,   
        radius,    
        geoArea,   
        volume,    
        v_cav,    
        v_eca,    
        v_cyt,    
        vr_cav,    
        vr_eca,    
        vr_cyt,    
        j_cav_eca,    
        j_cav_cyt,    
        j_eca_cyt,    
        PKA_tot,      
        f_cav,    
        f_eca,    
        f_cyt,    
        PKA_eca,    
        PKA_cav,    
        PKA_cyt,    
        PKI_tot,    
        f_pki_cav,    
        f_pki_eca,    
        f_pki_cyt,    
        PKI_cav,    
        PKI_eca,    
        PKI_cyt,    
        f_pki,    
        K_pki,    
        b_pki,    
        pka_cav_f1,    
        pka_cav_f2,    
        pka_cav_f3,    
        pka_cav_K1,    
        pka_cav_K2,    
        pka_cav_K3,    
        pka_cav_b1,    
        pka_cav_b2,    
        pka_cav_b3,    
        pka_eca_f1,    
        pka_eca_f2,    
        pka_eca_f3,    
        pka_eca_K1,    
        pka_eca_K2,    
        pka_eca_K3,    
        pka_eca_b1,    
        pka_eca_b2,    
        pka_eca_b3,    
        pka_cyt_f1,    
        pka_cyt_f2,    
        pka_cyt_f3,    
        pka_cyt_K1,    
        pka_cyt_K2,    
        pka_cyt_K3,    
        pka_cyt_b1,    
        pka_cyt_b2,    
        pka_cyt_b3,    
        f,    
        kp,    
        kdp,   
        Kp,    
        Kdp,    
        PP1_eca,    
        PP1_cav,    
        PP1_cyt,    
        pp1_K,      
        inhib1_tot,    
        PP2A,    
        PDE2_tot,    
        f_pde2_cav,    
        f_pde2_eca,    
        f_pde2_cyt,    
        f_pde2_part,   
        f_pde_part,    
        f_pde4_part,   
        f_pde4_cav,    
        f_pde4_eca,   
        f_pde4_cyt,    
        kPDE2,    
        kPDE3,    
        kPDE4,    
        KmPDE2,    
        KmPDE3,    
        KmPDE4,    
        KmIbmxPde2,    
        h_ibmx_pde2,    
        h_ibmx_pde3,    
        h_ibmx_pde4,    
        KmIbmxPde3,     
        KmIbmxPde4,    
        KPDEp,    
        delta_k_pde34,  
        ff_pde3_cyt,    
        kfPDEp,    
        r_pde34_frac, 
        r_pde3_cyt,
        kbPDEp,    
        pde_PDE3_tot_alpha,   
        pde_PDE3_tot_beta,    
        PDE3_tot,    
        PDE4_tot,    
        ibmx_h2,    
        ibmx_h3,    
        ibmx_h4,    
        ibmx2,    
        ibmx3,    
        ibmx4,    
        f_pde3_cav,    
        f_pde3_cyt,    
        PDE2_cav,    
        PDE2_eca,    
        PDE2_cyt,    
        PDE3_cav,    
        PDE3_cyt,    
        PDE4_cav,    
        PDE4_eca,    
        PDE4_cyt,    
        beta_R_b1_tot,    
        beta_R_b2_tot,    
        Gs_tot,    
        Gi_tot,    
        f_Gs_eca,    
        f_Gs_cav,    
        f_Gs_cyt,    
        f_Gi_cav,    
        f_Gi_eca,    
        f_Rb1_cav,    
        f_Rb1_eca,    
        f_Rb1_cyt,    
        f_Rb2_cav,    
        f_Rb2_eca,    
        k_b1_l,    
        k_b1_c,    
        k_b1_h,    
        k_b2_n,    
        k_b2_h,    
        k_b2_a,    
        k_b2_c,    
        k_b2_l,    
        k_b2_f,    
        k_act1_Gs,  
        k_act2_Gs,  
        k_act1_Gi,  
        k_act2_Gi,  
        k_hydr_Gs,  
        k_hydr_Gi,  
        k_reas_Gs,  
        k_reas_Gi,  
        rate_bds,   
        k_grk_dp,   
        k_grk_p,    
        k_pka_p,    
        k_pka_dp,    
        beta_cav_GRK,    
        beta_cav_R_b1_tot,    
        beta_cav_k_GsAct_b2,  
        beta_cav_R_b2_tot,    
        beta_cav_Rb2_pka_f_a, 
        beta_cav_Gs_f_c22,
        beta_cav_Gs_f_a,  
        beta_cav_Gs_f_c33,
        beta_cav_Gs_f_c11,
        beta_eca_GRK,    
        beta_eca_k_GsAct_b2, 
        beta_eca_R_b2_tot,   
        beta_eca_R_b1_tot,   
        beta_eca_Rb2_pka_f_a, 
        beta_eca_Gs_f_c11, 
        beta_eca_Gs_f_c33, 
        beta_eca_Gs_f_a,   
        beta_eca_Gs_f_c22,   
        beta_cyt_R_b1_tot,   
        beta_cyt_GRK,    
        beta_cyt_Rb1_np_f_a, 
        ATP,    
        KmATP,    
        AC_tot,    
        f_AC47_eca,  
        f_AC56_cav,  
        f_AC56_AC47, 
        hGsAC47,   
        hGsAC56,   
        hGsGiAC56, 
        KmGsAC47,  
        KmGsAC56,  
        KmGiAC56,  
        KmGsGiAC56,
        basalAC47, 
        basalAC56, 
        afAC47,   
        afAC56,   
        vGsGiAC56, 
        AC47_cyt, 
        AC56_cav, 
        fATP,   
        AC47_eca, 
        AC56_cyt,
        ka_inak, 
        kp_inak, 
        Ka_inak, 
        Kp_inak, 
        Ka_ina,  
        Kp_ina,  
        ka_ina,  
        kp_ina,  
        ka_plb,  
        kp_plb,  
        Ka_plb,  
        Kp_plb,  
        ka_ikur, 
        kp_ikur, 
        Ka_ikur, 
        Kp_ikur, 
        ka_iks,  
        kp_iks,  
        Ka_iks,  
        Kp_iks,  
        M,    
        iks_sig_L, 
        iks_sig_K, 
        Yotiao,    
        IKs_tot,    
        iks_sig_PKAf_sum,
        iks_sig_PKAf,    
        iks_sig_PP1f_eca_sum,
        PP1f_eca,    
        iks_sig_IKsf_sum,
        IKsf,    
        Yotiaof,  
        IKs_arn,  
        IKs_arp,  
        ka_tni,   
        kp_tni,   
        Ka_tni,   
        Kp_tni,   
        RyR_tot,  
        RyR_akap, 
        ka_ryr,   
        kp_ryr,   
        Ka_ryr,   
        Kp_ryr,   
        Mr,    
        Lr,    
        Kr,    
        akap_sig_RyRf_sum,
        RyRf,    
        ICaL_tot, 
        ICaL_akap,
        ka_ical, 
        kp_ical, 
        Ka_ical,
        Kp_ical,
        Mi,    
        Li,    
        Ki,    
        akap_sig_ICaLf_sum,
        ICaLf,    
        akap_sig_PP1f_cav_b,    
        akap_sig_PP1f_cav_c,    
        akap_sig_PP1f_cav_d,    
        akap_sig_PP1f_cav_rr,   
        akap_sig_PP1f_cav_yi,   
        akap_sig_PP1f_cav_yr,   
        akap_sig_PP1f_cav_mag,  
        akap_sig_PP1f_cav_arg,  
        akap_sig_PP1f_cav_x,    
        PP1f_cav,    
        akap_sig_PKAf_d,   
        akap_sig_PKAf_b,   
        akap_sig_PKAf_c,   
        akap_sig_PKAf_rr,  
        akap_sig_PKAf_yr,  
        akap_sig_PKAf_yi,  
        akap_sig_PKAf_mag, 
        akap_sig_PKAf_arg, 
        akap_sig_PKAf_x,
        akap_sig_PKAf, 
        RyR_akapf, 
        RyR_arn,   
        RyR_arp,   
        ICaL_akapf,
        ICaL_arn,  
        ICaL_arp,  
        beta_0,    
        irel_fhat_ratio,    
        ical_f_hat_ratio,   
        iks_f_hat_ratio,
        Whole_cell_PP1,
        PP1block,
        fICaLP,
        fIKsP,
        fPLBP,
        fTnIP,
        fINaP,
        fINaKP,
        fRyRP,
        fIKurP
      };
} // End of namespace Gng20Prm


/**
 * \namespace ELECTRA::Gng20Cur
 * \author Konstantinos A. Mountris
 * \brief Encaplulates enumeration to access O'Hara Rudy '11 human ventricular action potential model currents. 
 */
namespace Gng20Cur {
enum {
       INa,
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
       INaK,
       IKb,
       INab,
       ICab,
       IpCa,
       Iion
      };
} // End of namespace Gng20Cur


/**
 * \class Gong2020
 * \author Konstantinos A. Mountris
 * \brief  Gong et al 2020, human cardiac ventricular action potential model with \beta-adrenergic singaling \cite gong2020.
 */

class Gong2020 : public EpBasic
{
protected:

    virtual void SetDataMapping();

    virtual void ComputeBetaAdrenergicSignaling(double dt);

    virtual void ComputeEffectiveFraction();

    virtual void ComputeElectrophysiology(double v_new, double dt, double stim_current);

public:

    /**
     * \brief The default constructor of the Gong2020 class.
     */
    Gong2020();


    /**
     * \brief The default destructor of the Gong2020 class.
     */
    virtual ~Gong2020();


    /**
     * \brief Initialize the variables and parameters of the Gong2020 cell model.
     * \return [void]
     */
    virtual void Initialize(CellType cell_type);


    /**
     * \brief Compute the temporal update of the Gong2020 cell model.
     * \return [void]
     */
    virtual void Compute(double v_new, double dt, double stim_current);


    /**
     * \brief Set the cell membrane potential.
     * \param [in] v The cell membrane potential value.
     * \return [void]
     */
    virtual void SetV(double v) { this->var_[Gng20Var::v] = v; }


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
    inline virtual double V() const { return this->var_[Gng20Var::v]; }


    /**
     * \brief Get the membrane potential time derivative of the cell model.
     * \return [double] The membrane potential time derivative.
     */
    inline virtual double dVdt() const { return this->var_[Gng20Var::dvdt]; }


    /**
     * \brief Get the total ionic current.
     * \return [double] i_ion The total ionic current.
     */
    inline double Iion() const { return this->cur_[Gng20Cur::Iion]; }

};

/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ELECTROPHYSIOLOGY_MODELS_GONG2020_HPP_