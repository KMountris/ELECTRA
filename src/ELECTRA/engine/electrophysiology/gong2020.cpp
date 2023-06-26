/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#include "ELECTRA/engine/electrophysiology/gong2020.hpp"


namespace ELECTRA {

void Gong2020::SetDataMapping()
{
    using namespace Gng20Var;
    using namespace Gng20Prm;
    using namespace Gng20Cur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["nai"] = static_cast<std::size_t>(nai);
    this->mapped_data_["nass"] = static_cast<std::size_t>(nass);
    this->mapped_data_["ki"] = static_cast<std::size_t>(ki);
    this->mapped_data_["kss"] = static_cast<std::size_t>(kss);
    this->mapped_data_["cai"] = static_cast<std::size_t>(cai);
    this->mapped_data_["cass"] = static_cast<std::size_t>(cass);
    this->mapped_data_["cansr"] = static_cast<std::size_t>(cansr);
    this->mapped_data_["cajsr"] = static_cast<std::size_t>(cajsr);
    this->mapped_data_["m"] = static_cast<std::size_t>(m);
    this->mapped_data_["hf"] = static_cast<std::size_t>(hf);
    this->mapped_data_["hs"] = static_cast<std::size_t>(hs);
    this->mapped_data_["j"] = static_cast<std::size_t>(j);
    this->mapped_data_["hsp"] = static_cast<std::size_t>(hsp);
    this->mapped_data_["jp"] = static_cast<std::size_t>(jp);
    this->mapped_data_["mL"] = static_cast<std::size_t>(mL);
    this->mapped_data_["hL"] = static_cast<std::size_t>(hL);
    this->mapped_data_["hLp"] = static_cast<std::size_t>(hLp);
    this->mapped_data_["a"] = static_cast<std::size_t>(a);
    this->mapped_data_["iF"] = static_cast<std::size_t>(iF);
    this->mapped_data_["iS"] = static_cast<std::size_t>(iS);
    this->mapped_data_["ap"] = static_cast<std::size_t>(ap);
    this->mapped_data_["iFp"] = static_cast<std::size_t>(iFp);
    this->mapped_data_["iSp"] = static_cast<std::size_t>(iSp);
    this->mapped_data_["d"] = static_cast<std::size_t>(d);
    this->mapped_data_["ff"] = static_cast<std::size_t>(ff);
    this->mapped_data_["fs"] = static_cast<std::size_t>(fs);
    this->mapped_data_["fcaf"] = static_cast<std::size_t>(fcaf);
    this->mapped_data_["fcas"] = static_cast<std::size_t>(fcas);
    this->mapped_data_["jca"] = static_cast<std::size_t>(jca);
    this->mapped_data_["nca"] = static_cast<std::size_t>(nca);
    this->mapped_data_["ffp"] = static_cast<std::size_t>(ffp);
    this->mapped_data_["fcafp"] = static_cast<std::size_t>(fcafp);
    this->mapped_data_["xrf"] = static_cast<std::size_t>(xrf);
    this->mapped_data_["xrs"] = static_cast<std::size_t>(xrs);
    this->mapped_data_["xs1"] = static_cast<std::size_t>(xs1);
    this->mapped_data_["xs2"] = static_cast<std::size_t>(xs2);
    this->mapped_data_["xk1"] = static_cast<std::size_t>(xk1);
    this->mapped_data_["Jrelnp"] = static_cast<std::size_t>(Jrelnp);
    this->mapped_data_["Jrelp"] = static_cast<std::size_t>(Jrelp);
    this->mapped_data_["CaMKt"] = static_cast<std::size_t>(CaMKt);
    this->mapped_data_["hPf"] = static_cast<std::size_t>(hPf);
    this->mapped_data_["hPs"] = static_cast<std::size_t>(hPs);
    this->mapped_data_["jP"] = static_cast<std::size_t>(jP);
    this->mapped_data_["hBPf"] = static_cast<std::size_t>(hBPf);
    this->mapped_data_["hBPs"] = static_cast<std::size_t>(hBPs);
    this->mapped_data_["jBP"] = static_cast<std::size_t>(jBP);
    this->mapped_data_["dP"] = static_cast<std::size_t>(dP);
    this->mapped_data_["fPf"] = static_cast<std::size_t>(fPf);
    this->mapped_data_["fPs"] = static_cast<std::size_t>(fPs);
    this->mapped_data_["fcaPf"] = static_cast<std::size_t>(fcaPf);
    this->mapped_data_["fcaPs"] = static_cast<std::size_t>(fcaPs);
    this->mapped_data_["fBPf"] = static_cast<std::size_t>(fBPf);
    this->mapped_data_["fcaBPf"] = static_cast<std::size_t>(fcaBPf);
    this->mapped_data_["xs1P"] = static_cast<std::size_t>(xs1P);
    this->mapped_data_["xs2P"] = static_cast<std::size_t>(xs2P);
    this->mapped_data_["JrelP"] = static_cast<std::size_t>(JrelP);
    this->mapped_data_["JrelBP"] = static_cast<std::size_t>(JrelBP);
    this->mapped_data_["beta_cav_Gs_aGTP"] = static_cast<std::size_t>(beta_cav_Gs_aGTP);
    this->mapped_data_["beta_eca_Gs_aGTP"] = static_cast<std::size_t>(beta_eca_Gs_aGTP);
    this->mapped_data_["beta_cyt_Gs_aGTP"] = static_cast<std::size_t>(beta_cyt_Gs_aGTP);
    this->mapped_data_["beta_cav_Gs_bg"] = static_cast<std::size_t>(beta_cav_Gs_bg);
    this->mapped_data_["beta_eca_Gs_bg"] = static_cast<std::size_t>(beta_eca_Gs_bg);
    this->mapped_data_["beta_cyt_Gs_bg"] = static_cast<std::size_t>(beta_cyt_Gs_bg);
    this->mapped_data_["beta_cav_Gs_aGDP"] = static_cast<std::size_t>(beta_cav_Gs_aGDP);
    this->mapped_data_["beta_eca_Gs_aGDP"] = static_cast<std::size_t>(beta_eca_Gs_aGDP);
    this->mapped_data_["beta_cyt_Gs_aGDP"] = static_cast<std::size_t>(beta_cyt_Gs_aGDP);
    this->mapped_data_["cAMP_cav"] = static_cast<std::size_t>(cAMP_cav);
    this->mapped_data_["cAMP_eca"] = static_cast<std::size_t>(cAMP_eca);
    this->mapped_data_["cAMP_cyt"] = static_cast<std::size_t>(cAMP_cyt);
    this->mapped_data_["beta_cav_Rb1_pka_tot"] = static_cast<std::size_t>(beta_cav_Rb1_pka_tot);
    this->mapped_data_["beta_eca_Rb1_pka_tot"] = static_cast<std::size_t>(beta_eca_Rb1_pka_tot);
    this->mapped_data_["beta_cyt_Rb1_pka_tot"] = static_cast<std::size_t>(beta_cyt_Rb1_pka_tot);
    this->mapped_data_["beta_cav_Rb1_grk_tot"] = static_cast<std::size_t>(beta_cav_Rb1_grk_tot);
    this->mapped_data_["beta_eca_Rb1_grk_tot"] = static_cast<std::size_t>(beta_eca_Rb1_grk_tot);
    this->mapped_data_["beta_cyt_Rb1_grk_tot"] = static_cast<std::size_t>(beta_cyt_Rb1_grk_tot);
    this->mapped_data_["pka_cav_ARC"] = static_cast<std::size_t>(pka_cav_ARC);
    this->mapped_data_["pka_cav_A2RC"] = static_cast<std::size_t>(pka_cav_A2RC);
    this->mapped_data_["pka_cav_A2R"] = static_cast<std::size_t>(pka_cav_A2R);
    this->mapped_data_["pka_cav_C"] = static_cast<std::size_t>(pka_cav_C);
    this->mapped_data_["pka_cav_PKIC"] = static_cast<std::size_t>(pka_cav_PKIC);
    this->mapped_data_["pka_eca_ARC"] = static_cast<std::size_t>(pka_eca_ARC);
    this->mapped_data_["pka_eca_A2RC"] = static_cast<std::size_t>(pka_eca_A2RC);
    this->mapped_data_["pka_eca_A2R"] = static_cast<std::size_t>(pka_eca_A2R);
    this->mapped_data_["pka_eca_C"] = static_cast<std::size_t>(pka_eca_C);
    this->mapped_data_["pka_eca_PKIC"] = static_cast<std::size_t>(pka_eca_PKIC);
    this->mapped_data_["pka_cyt_ARC"] = static_cast<std::size_t>(pka_cyt_ARC);
    this->mapped_data_["pka_cyt_A2RC"] = static_cast<std::size_t>(pka_cyt_A2RC);
    this->mapped_data_["pka_cyt_A2R"] = static_cast<std::size_t>(pka_cyt_A2R);
    this->mapped_data_["pka_cyt_C"] = static_cast<std::size_t>(pka_cyt_C);
    this->mapped_data_["pka_cyt_PKIC"] = static_cast<std::size_t>(pka_cyt_PKIC);
    this->mapped_data_["PDE3_P_cav"] = static_cast<std::size_t>(PDE3_P_cav);
    this->mapped_data_["PDE3_P_cyt"] = static_cast<std::size_t>(PDE3_P_cyt);
    this->mapped_data_["PDE4_P_cav"] = static_cast<std::size_t>(PDE4_P_cav);
    this->mapped_data_["PDE4_P_eca"] = static_cast<std::size_t>(PDE4_P_eca);
    this->mapped_data_["PDE4_P_cyt"] = static_cast<std::size_t>(PDE4_P_cyt);
    this->mapped_data_["inhib1_p"] = static_cast<std::size_t>(inhib1_p);
    this->mapped_data_["ICaLp"] = static_cast<std::size_t>(ICaLp);
    this->mapped_data_["IKsp"] = static_cast<std::size_t>(IKsp);
    this->mapped_data_["iup_f_plb"] = static_cast<std::size_t>(iup_f_plb);
    this->mapped_data_["f_tni"] = static_cast<std::size_t>(f_tni);
    this->mapped_data_["ina_f_ina"] = static_cast<std::size_t>(ina_f_ina);
    this->mapped_data_["f_inak"] = static_cast<std::size_t>(f_inak);
    this->mapped_data_["RyRp"] = static_cast<std::size_t>(RyRp);
    this->mapped_data_["f_ikur"] = static_cast<std::size_t>(f_ikur);
    this->mapped_data_["beta_cav_Rb2_pka_tot"] = static_cast<std::size_t>(beta_cav_Rb2_pka_tot);
    this->mapped_data_["beta_cav_Rb2_grk_tot"] = static_cast<std::size_t>(beta_cav_Rb2_grk_tot);
    this->mapped_data_["beta_cav_Gi_aGTP"] = static_cast<std::size_t>(beta_cav_Gi_aGTP);
    this->mapped_data_["beta_cav_Gi_bg"] = static_cast<std::size_t>(beta_cav_Gi_bg);
    this->mapped_data_["beta_cav_Gi_aGDP"] = static_cast<std::size_t>(beta_cav_Gi_aGDP);
    this->mapped_data_["beta_eca_Rb2_pka_tot"] = static_cast<std::size_t>(beta_eca_Rb2_pka_tot);
    this->mapped_data_["beta_eca_Rb2_grk_tot"] = static_cast<std::size_t>(beta_eca_Rb2_grk_tot);
    this->mapped_data_["beta_eca_Gi_aGTP"] = static_cast<std::size_t>(beta_eca_Gi_aGTP);
    this->mapped_data_["beta_eca_Gi_bg"] = static_cast<std::size_t>(beta_eca_Gi_bg);
    this->mapped_data_["beta_eca_Gi_aGDP"] = static_cast<std::size_t>(beta_eca_Gi_aGDP);

    // Set parameters mapping.
    this->mapped_data_["signaling"] = static_cast<std::size_t>(signaling);
    this->mapped_data_["electrophys"] = static_cast<std::size_t>(electrophys);
    this->mapped_data_["cipa"] = static_cast<std::size_t>(cipa);
    this->mapped_data_["GKS_mul"] = static_cast<std::size_t>(GKS_mul);
    this->mapped_data_["VP_d"] = static_cast<std::size_t>(VP_d);
    this->mapped_data_["VP_f"] = static_cast<std::size_t>(VP_f);
    this->mapped_data_["GCal_mul"] = static_cast<std::size_t>(GCal_mul);
    this->mapped_data_["CiPA_scale_IKr"] = static_cast<std::size_t>(CiPA_scale_IKr);
    this->mapped_data_["CiPA_scale_IKs"] = static_cast<std::size_t>(CiPA_scale_IKs);
    this->mapped_data_["CiPA_scale_IK1"] = static_cast<std::size_t>(CiPA_scale_IK1);
    this->mapped_data_["CiPA_scale_ICaL"] = static_cast<std::size_t>(CiPA_scale_ICaL);
    this->mapped_data_["CiPA_scale_INaL"] = static_cast<std::size_t>(CiPA_scale_INaL);
    this->mapped_data_["nao"] = static_cast<std::size_t>(nao);
    this->mapped_data_["cao"] = static_cast<std::size_t>(cao);
    this->mapped_data_["ko"] = static_cast<std::size_t>(ko);
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["F"] = static_cast<std::size_t>(F);
    this->mapped_data_["L"] = static_cast<std::size_t>(L);
    this->mapped_data_["rad"] = static_cast<std::size_t>(rad);
    this->mapped_data_["vcell"] = static_cast<std::size_t>(vcell);
    this->mapped_data_["Ageo"] = static_cast<std::size_t>(Ageo);
    this->mapped_data_["Acap"] = static_cast<std::size_t>(Acap);
    this->mapped_data_["vmyo"] = static_cast<std::size_t>(vmyo);
    this->mapped_data_["vnsr"] = static_cast<std::size_t>(vnsr);
    this->mapped_data_["vjsr"] = static_cast<std::size_t>(vjsr);
    this->mapped_data_["vss"] = static_cast<std::size_t>(vss);
    this->mapped_data_["aCaMK"] = static_cast<std::size_t>(aCaMK);
    this->mapped_data_["bCaMK"] = static_cast<std::size_t>(bCaMK);
    this->mapped_data_["CaMKo"] = static_cast<std::size_t>(CaMKo);
    this->mapped_data_["PKNa"] = static_cast<std::size_t>(PKNa);
    this->mapped_data_["GNa"] = static_cast<std::size_t>(GNa);
    this->mapped_data_["GNaL"] = static_cast<std::size_t>(GNaL);
    this->mapped_data_["Gto"] = static_cast<std::size_t>(Gto);
    this->mapped_data_["GKr"] = static_cast<std::size_t>(GKr);
    this->mapped_data_["GKs"] = static_cast<std::size_t>(GKs);
    this->mapped_data_["GK1"] = static_cast<std::size_t>(GK1);
    this->mapped_data_["Gncx"] = static_cast<std::size_t>(Gncx);
    this->mapped_data_["GKb"] = static_cast<std::size_t>(GKb);
    this->mapped_data_["GpCa"] = static_cast<std::size_t>(GpCa);
    this->mapped_data_["PCa"] = static_cast<std::size_t>(PCa);
    this->mapped_data_["Pnak"] = static_cast<std::size_t>(Pnak);
    this->mapped_data_["PNab"] = static_cast<std::size_t>(PNab);
    this->mapped_data_["PCab"] = static_cast<std::size_t>(PCab);
    this->mapped_data_["KmCaMK"] = static_cast<std::size_t>(KmCaMK);
    this->mapped_data_["KmCaM"] = static_cast<std::size_t>(KmCaM);
    this->mapped_data_["SERCA_total"] = static_cast<std::size_t>(SERCA_total);
    this->mapped_data_["RyR_total"] = static_cast<std::size_t>(RyR_total);
    this->mapped_data_["Leak_total"] = static_cast<std::size_t>(Leak_total);
    this->mapped_data_["Trans_total"] = static_cast<std::size_t>(Trans_total);
    this->mapped_data_["iso_L"] = static_cast<std::size_t>(iso_L);
    this->mapped_data_["ibmx"] = static_cast<std::size_t>(ibmx);
    this->mapped_data_["phys_T"] = static_cast<std::size_t>(phys_T);
    this->mapped_data_["phys_R"] = static_cast<std::size_t>(phys_R);
    this->mapped_data_["length"] = static_cast<std::size_t>(length);
    this->mapped_data_["cell_pi"] = static_cast<std::size_t>(cell_pi);
    this->mapped_data_["radius"] = static_cast<std::size_t>(radius);
    this->mapped_data_["geoArea"] = static_cast<std::size_t>(geoArea);
    this->mapped_data_["volume"] = static_cast<std::size_t>(volume);
    this->mapped_data_["v_cav"] = static_cast<std::size_t>(v_cav);
    this->mapped_data_["v_eca"] = static_cast<std::size_t>(v_eca);
    this->mapped_data_["v_cyt"] = static_cast<std::size_t>(v_cyt);
    this->mapped_data_["vr_cav"] = static_cast<std::size_t>(vr_cav);
    this->mapped_data_["vr_eca"] = static_cast<std::size_t>(vr_eca);
    this->mapped_data_["vr_cyt"] = static_cast<std::size_t>(vr_cyt);
    this->mapped_data_["j_cav_eca"] = static_cast<std::size_t>(j_cav_eca);
    this->mapped_data_["j_cav_cyt"] = static_cast<std::size_t>(j_cav_cyt);
    this->mapped_data_["j_eca_cyt"] = static_cast<std::size_t>(j_eca_cyt);
    this->mapped_data_["PKA_tot"] = static_cast<std::size_t>(PKA_tot);
    this->mapped_data_["f_cav"] = static_cast<std::size_t>(f_cav);
    this->mapped_data_["f_eca"] = static_cast<std::size_t>(f_eca);
    this->mapped_data_["f_cyt"] = static_cast<std::size_t>(f_cyt);
    this->mapped_data_["PKA_eca"] = static_cast<std::size_t>(PKA_eca);
    this->mapped_data_["PKA_cav"] = static_cast<std::size_t>(PKA_cav);
    this->mapped_data_["PKA_cyt"] = static_cast<std::size_t>(PKA_cyt);
    this->mapped_data_["PKI_tot"] = static_cast<std::size_t>(PKI_tot);
    this->mapped_data_["f_pki_cav"] = static_cast<std::size_t>(f_pki_cav);
    this->mapped_data_["f_pki_eca"] = static_cast<std::size_t>(f_pki_eca);
    this->mapped_data_["f_pki_cyt"] = static_cast<std::size_t>(f_pki_cyt);
    this->mapped_data_["PKI_cav"] = static_cast<std::size_t>(PKI_cav);
    this->mapped_data_["PKI_eca"] = static_cast<std::size_t>(PKI_eca);
    this->mapped_data_["PKI_cyt"] = static_cast<std::size_t>(PKI_cyt);
    this->mapped_data_["f_pki"] = static_cast<std::size_t>(f_pki);
    this->mapped_data_["K_pki"] = static_cast<std::size_t>(K_pki);
    this->mapped_data_["b_pki"] = static_cast<std::size_t>(b_pki);
    this->mapped_data_["pka_cav_f1"] = static_cast<std::size_t>(pka_cav_f1);
    this->mapped_data_["pka_cav_f2"] = static_cast<std::size_t>(pka_cav_f2);
    this->mapped_data_["pka_cav_f3"] = static_cast<std::size_t>(pka_cav_f3);
    this->mapped_data_["pka_cav_K1"] = static_cast<std::size_t>(pka_cav_K1);
    this->mapped_data_["pka_cav_K2"] = static_cast<std::size_t>(pka_cav_K2);
    this->mapped_data_["pka_cav_K3"] = static_cast<std::size_t>(pka_cav_K3);
    this->mapped_data_["pka_cav_b1"] = static_cast<std::size_t>(pka_cav_b1);
    this->mapped_data_["pka_cav_b2"] = static_cast<std::size_t>(pka_cav_b2);
    this->mapped_data_["pka_cav_b3"] = static_cast<std::size_t>(pka_cav_b3);
    this->mapped_data_["pka_eca_f1"] = static_cast<std::size_t>(pka_eca_f1);
    this->mapped_data_["pka_eca_f2"] = static_cast<std::size_t>(pka_eca_f2);
    this->mapped_data_["pka_eca_f3"] = static_cast<std::size_t>(pka_eca_f3);
    this->mapped_data_["pka_eca_K1"] = static_cast<std::size_t>(pka_eca_K1);
    this->mapped_data_["pka_eca_K2"] = static_cast<std::size_t>(pka_eca_K2);
    this->mapped_data_["pka_eca_K3"] = static_cast<std::size_t>(pka_eca_K3);
    this->mapped_data_["pka_eca_b1"] = static_cast<std::size_t>(pka_eca_b1);
    this->mapped_data_["pka_eca_b2"] = static_cast<std::size_t>(pka_eca_b2);
    this->mapped_data_["pka_eca_b3"] = static_cast<std::size_t>(pka_eca_b3);
    this->mapped_data_["pka_cyt_f1"] = static_cast<std::size_t>(pka_cyt_f1);
    this->mapped_data_["pka_cyt_f2"] = static_cast<std::size_t>(pka_cyt_f2);
    this->mapped_data_["pka_cyt_f3"] = static_cast<std::size_t>(pka_cyt_f3);
    this->mapped_data_["pka_cyt_K1"] = static_cast<std::size_t>(pka_cyt_K1);
    this->mapped_data_["pka_cyt_K2"] = static_cast<std::size_t>(pka_cyt_K2);
    this->mapped_data_["pka_cyt_K3"] = static_cast<std::size_t>(pka_cyt_K3);
    this->mapped_data_["pka_cyt_b1"] = static_cast<std::size_t>(pka_cyt_b1);
    this->mapped_data_["pka_cyt_b2"] = static_cast<std::size_t>(pka_cyt_b2);
    this->mapped_data_["pka_cyt_b3"] = static_cast<std::size_t>(pka_cyt_b3);
    this->mapped_data_["f"] = static_cast<std::size_t>(f);
    this->mapped_data_["kp"] = static_cast<std::size_t>(kp);
    this->mapped_data_["kdp"] = static_cast<std::size_t>(kdp);
    this->mapped_data_["Kp"] = static_cast<std::size_t>(Kp);
    this->mapped_data_["Kdp"] = static_cast<std::size_t>(Kdp);
    this->mapped_data_["PP1_eca"] = static_cast<std::size_t>(PP1_eca);
    this->mapped_data_["PP1_cav"] = static_cast<std::size_t>(PP1_cav);
    this->mapped_data_["PP1_cyt"] = static_cast<std::size_t>(PP1_cyt);
    this->mapped_data_["pp1_K"] = static_cast<std::size_t>(pp1_K);
    this->mapped_data_["inhib1_tot"] = static_cast<std::size_t>(inhib1_tot);
    this->mapped_data_["PP2A"] = static_cast<std::size_t>(PP2A);
    this->mapped_data_["PDE2_tot"] = static_cast<std::size_t>(PDE2_tot);
    this->mapped_data_["f_pde2_cav"] = static_cast<std::size_t>(f_pde2_cav);
    this->mapped_data_["f_pde2_eca"] = static_cast<std::size_t>(f_pde2_eca);
    this->mapped_data_["f_pde2_cyt"] = static_cast<std::size_t>(f_pde2_cyt);
    this->mapped_data_["f_pde2_part"] = static_cast<std::size_t>(f_pde2_part);
    this->mapped_data_["f_pde_part"] = static_cast<std::size_t>(f_pde_part);
    this->mapped_data_["f_pde4_part"] = static_cast<std::size_t>(f_pde4_part);
    this->mapped_data_["f_pde4_cav"] = static_cast<std::size_t>(f_pde4_cav);
    this->mapped_data_["f_pde4_eca"] = static_cast<std::size_t>(f_pde4_eca);
    this->mapped_data_["f_pde4_cyt"] = static_cast<std::size_t>(f_pde4_cyt);
    this->mapped_data_["kPDE2"] = static_cast<std::size_t>(kPDE2);
    this->mapped_data_["kPDE3"] = static_cast<std::size_t>(kPDE3);
    this->mapped_data_["kPDE4"] = static_cast<std::size_t>(kPDE4);
    this->mapped_data_["KmPDE2"] = static_cast<std::size_t>(KmPDE2);
    this->mapped_data_["KmPDE3"] = static_cast<std::size_t>(KmPDE3);
    this->mapped_data_["KmPDE4"] = static_cast<std::size_t>(KmPDE4);
    this->mapped_data_["KmIbmxPde2"] = static_cast<std::size_t>(KmIbmxPde2);
    this->mapped_data_["h_ibmx_pde2"] = static_cast<std::size_t>(h_ibmx_pde2);
    this->mapped_data_["h_ibmx_pde3"] = static_cast<std::size_t>(h_ibmx_pde3);
    this->mapped_data_["h_ibmx_pde4"] = static_cast<std::size_t>(h_ibmx_pde4);
    this->mapped_data_["KmIbmxPde3"] = static_cast<std::size_t>(KmIbmxPde3);
    this->mapped_data_["KmIbmxPde4"] = static_cast<std::size_t>(KmIbmxPde4);
    this->mapped_data_["KPDEp"] = static_cast<std::size_t>(KPDEp);
    this->mapped_data_["delta_k_pde34"] = static_cast<std::size_t>(delta_k_pde34);
    this->mapped_data_["ff_pde3_cyt"] = static_cast<std::size_t>(ff_pde3_cyt);
    this->mapped_data_["kfPDEp"] = static_cast<std::size_t>(kfPDEp);
    this->mapped_data_["r_pde34_frac"] = static_cast<std::size_t>(r_pde34_frac);
    this->mapped_data_["r_pde3_cyt"] = static_cast<std::size_t>(r_pde3_cyt);
    this->mapped_data_["kbPDEp"] = static_cast<std::size_t>(kbPDEp);
    this->mapped_data_["pde_PDE3_tot_alpha"] = static_cast<std::size_t>(pde_PDE3_tot_alpha);
    this->mapped_data_["pde_PDE3_tot_beta"] = static_cast<std::size_t>(pde_PDE3_tot_beta);
    this->mapped_data_["PDE3_tot"] = static_cast<std::size_t>(PDE3_tot);
    this->mapped_data_["PDE4_tot"] = static_cast<std::size_t>(PDE4_tot);
    this->mapped_data_["ibmx_h2"] = static_cast<std::size_t>(ibmx_h2);
    this->mapped_data_["ibmx_h3"] = static_cast<std::size_t>(ibmx_h3);
    this->mapped_data_["ibmx_h4"] = static_cast<std::size_t>(ibmx_h4);
    this->mapped_data_["ibmx2"] = static_cast<std::size_t>(ibmx2);
    this->mapped_data_["ibmx3"] = static_cast<std::size_t>(ibmx3);
    this->mapped_data_["ibmx4"] = static_cast<std::size_t>(ibmx4);
    this->mapped_data_["f_pde3_cav"] = static_cast<std::size_t>(f_pde3_cav);
    this->mapped_data_["f_pde3_cyt"] = static_cast<std::size_t>(f_pde3_cyt);
    this->mapped_data_["PDE2_cav"] = static_cast<std::size_t>(PDE2_cav);
    this->mapped_data_["PDE2_eca"] = static_cast<std::size_t>(PDE2_eca);
    this->mapped_data_["PDE2_cyt"] = static_cast<std::size_t>(PDE2_cyt);
    this->mapped_data_["PDE3_cav"] = static_cast<std::size_t>(PDE3_cav);
    this->mapped_data_["PDE3_cyt"] = static_cast<std::size_t>(PDE3_cyt);
    this->mapped_data_["PDE4_cav"] = static_cast<std::size_t>(PDE4_cav);
    this->mapped_data_["PDE4_eca"] = static_cast<std::size_t>(PDE4_eca);
    this->mapped_data_["PDE4_cyt"] = static_cast<std::size_t>(PDE4_cyt);
    this->mapped_data_["beta_R_b1_tot"] = static_cast<std::size_t>(beta_R_b1_tot);
    this->mapped_data_["beta_R_b2_tot"] = static_cast<std::size_t>(beta_R_b2_tot);
    this->mapped_data_["Gs_tot"] = static_cast<std::size_t>(Gs_tot);
    this->mapped_data_["Gi_tot"] = static_cast<std::size_t>(Gi_tot);
    this->mapped_data_["f_Gs_eca"] = static_cast<std::size_t>(f_Gs_eca);
    this->mapped_data_["f_Gs_cav"] = static_cast<std::size_t>(f_Gs_cav);
    this->mapped_data_["f_Gs_cyt"] = static_cast<std::size_t>(f_Gs_cyt);
    this->mapped_data_["f_Gi_cav"] = static_cast<std::size_t>(f_Gi_cav);
    this->mapped_data_["f_Gi_eca"] = static_cast<std::size_t>(f_Gi_eca);
    this->mapped_data_["f_Rb1_cav"] = static_cast<std::size_t>(f_Rb1_cav);
    this->mapped_data_["f_Rb1_eca"] = static_cast<std::size_t>(f_Rb1_eca);
    this->mapped_data_["f_Rb1_cyt"] = static_cast<std::size_t>(f_Rb1_cyt);
    this->mapped_data_["f_Rb2_cav"] = static_cast<std::size_t>(f_Rb2_cav);
    this->mapped_data_["f_Rb2_eca"] = static_cast<std::size_t>(f_Rb2_eca);
    this->mapped_data_["k_b1_l"] = static_cast<std::size_t>(k_b1_l);
    this->mapped_data_["k_b1_c"] = static_cast<std::size_t>(k_b1_c);
    this->mapped_data_["k_b1_h"] = static_cast<std::size_t>(k_b1_h);
    this->mapped_data_["k_b2_n"] = static_cast<std::size_t>(k_b2_n);
    this->mapped_data_["k_b2_h"] = static_cast<std::size_t>(k_b2_h);
    this->mapped_data_["k_b2_a"] = static_cast<std::size_t>(k_b2_a);
    this->mapped_data_["k_b2_c"] = static_cast<std::size_t>(k_b2_c);
    this->mapped_data_["k_b2_l"] = static_cast<std::size_t>(k_b2_l);
    this->mapped_data_["k_b2_f"] = static_cast<std::size_t>(k_b2_f);
    this->mapped_data_["k_act1_Gs"] = static_cast<std::size_t>(k_act1_Gs);
    this->mapped_data_["k_act2_Gs"] = static_cast<std::size_t>(k_act2_Gs);
    this->mapped_data_["k_act1_Gi"] = static_cast<std::size_t>(k_act1_Gi);
    this->mapped_data_["k_act2_Gi"] = static_cast<std::size_t>(k_act2_Gi);
    this->mapped_data_["k_hydr_Gs"] = static_cast<std::size_t>(k_hydr_Gs);
    this->mapped_data_["k_hydr_Gi"] = static_cast<std::size_t>(k_hydr_Gi);
    this->mapped_data_["k_reas_Gs"] = static_cast<std::size_t>(k_reas_Gs);
    this->mapped_data_["k_reas_Gi"] = static_cast<std::size_t>(k_reas_Gi);
    this->mapped_data_["rate_bds"] = static_cast<std::size_t>(rate_bds);
    this->mapped_data_["k_grk_dp"] = static_cast<std::size_t>(k_grk_dp);
    this->mapped_data_["k_grk_p"] = static_cast<std::size_t>(k_grk_p);
    this->mapped_data_["k_pka_p"] = static_cast<std::size_t>(k_pka_p);
    this->mapped_data_["k_pka_dp"] = static_cast<std::size_t>(k_pka_dp);
    this->mapped_data_["beta_cav_GRK"] = static_cast<std::size_t>(beta_cav_GRK);
    this->mapped_data_["beta_cav_R_b1_tot"] = static_cast<std::size_t>(beta_cav_R_b1_tot);
    this->mapped_data_["beta_cav_k_GsAct_b2"] = static_cast<std::size_t>(beta_cav_k_GsAct_b2);
    this->mapped_data_["beta_cav_R_b2_tot"] = static_cast<std::size_t>(beta_cav_R_b2_tot);
    this->mapped_data_["beta_cav_Rb2_pka_f_a"] = static_cast<std::size_t>(beta_cav_Rb2_pka_f_a);
    this->mapped_data_["beta_cav_Gs_f_c22"] = static_cast<std::size_t>(beta_cav_Gs_f_c22);
    this->mapped_data_["beta_cav_Gs_f_a"] = static_cast<std::size_t>(beta_cav_Gs_f_a);
    this->mapped_data_["beta_cav_Gs_f_c33"] = static_cast<std::size_t>(beta_cav_Gs_f_c33);
    this->mapped_data_["beta_cav_Gs_f_c11"] = static_cast<std::size_t>(beta_cav_Gs_f_c11);
    this->mapped_data_["beta_eca_GRK"] = static_cast<std::size_t>(beta_eca_GRK);
    this->mapped_data_["beta_eca_k_GsAct_b2"] = static_cast<std::size_t>(beta_eca_k_GsAct_b2);
    this->mapped_data_["beta_eca_R_b2_tot"] = static_cast<std::size_t>(beta_eca_R_b2_tot);
    this->mapped_data_["beta_eca_R_b1_tot"] = static_cast<std::size_t>(beta_eca_R_b1_tot);
    this->mapped_data_["beta_eca_Rb2_pka_f_a"] = static_cast<std::size_t>(beta_eca_Rb2_pka_f_a);
    this->mapped_data_["beta_eca_Gs_f_c11"] = static_cast<std::size_t>(beta_eca_Gs_f_c11);
    this->mapped_data_["beta_eca_Gs_f_c33"] = static_cast<std::size_t>(beta_eca_Gs_f_c33);
    this->mapped_data_["beta_eca_Gs_f_a"] = static_cast<std::size_t>(beta_eca_Gs_f_a);
    this->mapped_data_["beta_eca_Gs_f_c22"] = static_cast<std::size_t>(beta_eca_Gs_f_c22);
    this->mapped_data_["beta_cyt_R_b1_tot"] = static_cast<std::size_t>(beta_cyt_R_b1_tot);
    this->mapped_data_["beta_cyt_GRK"] = static_cast<std::size_t>(beta_cyt_GRK);
    this->mapped_data_["beta_cyt_Rb1_np_f_a"] = static_cast<std::size_t>(beta_cyt_Rb1_np_f_a);
    this->mapped_data_["ATP"] = static_cast<std::size_t>(ATP);
    this->mapped_data_["KmATP"] = static_cast<std::size_t>(KmATP);
    this->mapped_data_["AC_tot"] = static_cast<std::size_t>(AC_tot);
    this->mapped_data_["f_AC47_eca"] = static_cast<std::size_t>(f_AC47_eca);
    this->mapped_data_["f_AC56_cav"] = static_cast<std::size_t>(f_AC56_cav);
    this->mapped_data_["f_AC56_AC47"] = static_cast<std::size_t>(f_AC56_AC47);
    this->mapped_data_["hGsAC47"] = static_cast<std::size_t>(hGsAC47);
    this->mapped_data_["hGsAC56"] = static_cast<std::size_t>(hGsAC56);
    this->mapped_data_["hGsGiAC56"] = static_cast<std::size_t>(hGsGiAC56);
    this->mapped_data_["KmGsAC47"] = static_cast<std::size_t>(KmGsAC47);
    this->mapped_data_["KmGsAC56"] = static_cast<std::size_t>(KmGsAC56);
    this->mapped_data_["KmGiAC56"] = static_cast<std::size_t>(KmGiAC56);
    this->mapped_data_["KmGsGiAC56"] = static_cast<std::size_t>(KmGsGiAC56);
    this->mapped_data_["basalAC47"] = static_cast<std::size_t>(basalAC47);
    this->mapped_data_["basalAC56"] = static_cast<std::size_t>(basalAC56);
    this->mapped_data_["afAC47"] = static_cast<std::size_t>(afAC47);
    this->mapped_data_["afAC56"] = static_cast<std::size_t>(afAC56);
    this->mapped_data_["vGsGiAC56"] = static_cast<std::size_t>(vGsGiAC56);
    this->mapped_data_["AC47_cyt"] = static_cast<std::size_t>(AC47_cyt);
    this->mapped_data_["AC56_cav"] = static_cast<std::size_t>(AC56_cav);
    this->mapped_data_["fATP"] = static_cast<std::size_t>(fATP);
    this->mapped_data_["AC47_eca"] = static_cast<std::size_t>(AC47_eca);
    this->mapped_data_["AC56_cyt"] = static_cast<std::size_t>(AC56_cyt);
    this->mapped_data_["ka_inak"] = static_cast<std::size_t>(ka_inak);
    this->mapped_data_["kp_inak"] = static_cast<std::size_t>(kp_inak);
    this->mapped_data_["Ka_inak"] = static_cast<std::size_t>(Ka_inak);
    this->mapped_data_["Kp_inak"] = static_cast<std::size_t>(Kp_inak);
    this->mapped_data_["Ka_ina"] = static_cast<std::size_t>(Ka_ina);
    this->mapped_data_["Kp_ina"] = static_cast<std::size_t>(Kp_ina);
    this->mapped_data_["ka_ina"] = static_cast<std::size_t>(ka_ina);
    this->mapped_data_["kp_ina"] = static_cast<std::size_t>(kp_ina);
    this->mapped_data_["ka_plb"] = static_cast<std::size_t>(ka_plb);
    this->mapped_data_["kp_plb"] = static_cast<std::size_t>(kp_plb);
    this->mapped_data_["Ka_plb"] = static_cast<std::size_t>(Ka_plb);
    this->mapped_data_["Kp_plb"] = static_cast<std::size_t>(Kp_plb);
    this->mapped_data_["ka_ikur"] = static_cast<std::size_t>(ka_ikur);
    this->mapped_data_["kp_ikur"] = static_cast<std::size_t>(kp_ikur);
    this->mapped_data_["Ka_ikur"] = static_cast<std::size_t>(Ka_ikur);
    this->mapped_data_["Kp_ikur"] = static_cast<std::size_t>(Kp_ikur);
    this->mapped_data_["ka_iks"] = static_cast<std::size_t>(ka_iks);
    this->mapped_data_["kp_iks"] = static_cast<std::size_t>(kp_iks);
    this->mapped_data_["Ka_iks"] = static_cast<std::size_t>(Ka_iks);
    this->mapped_data_["Kp_iks"] = static_cast<std::size_t>(Kp_iks);
    this->mapped_data_["M"] = static_cast<std::size_t>(M);
    this->mapped_data_["iks_sig_L"] = static_cast<std::size_t>(iks_sig_L);
    this->mapped_data_["iks_sig_K"] = static_cast<std::size_t>(iks_sig_K);
    this->mapped_data_["Yotiao"] = static_cast<std::size_t>(Yotiao);
    this->mapped_data_["IKs_tot"] = static_cast<std::size_t>(IKs_tot);
    this->mapped_data_["iks_sig_PKAf_sum"] = static_cast<std::size_t>(iks_sig_PKAf_sum);
    this->mapped_data_["iks_sig_PKAf"] = static_cast<std::size_t>(iks_sig_PKAf);
    this->mapped_data_["iks_sig_PP1f_eca_sum"] = static_cast<std::size_t>(iks_sig_PP1f_eca_sum);
    this->mapped_data_["PP1f_eca"] = static_cast<std::size_t>(PP1f_eca);
    this->mapped_data_["iks_sig_IKsf_sum"] = static_cast<std::size_t>(iks_sig_IKsf_sum);
    this->mapped_data_["IKsf"] = static_cast<std::size_t>(IKsf);
    this->mapped_data_["Yotiaof"] = static_cast<std::size_t>(Yotiaof);
    this->mapped_data_["IKs_arn"] = static_cast<std::size_t>(IKs_arn);
    this->mapped_data_["IKs_arp"] = static_cast<std::size_t>(IKs_arp);
    this->mapped_data_["ka_tni"] = static_cast<std::size_t>(ka_tni);
    this->mapped_data_["kp_tni"] = static_cast<std::size_t>(kp_tni);
    this->mapped_data_["Ka_tni"] = static_cast<std::size_t>(Ka_tni);
    this->mapped_data_["Kp_tni"] = static_cast<std::size_t>(Kp_tni);
    this->mapped_data_["RyR_tot"] = static_cast<std::size_t>(RyR_tot);
    this->mapped_data_["RyR_akap"] = static_cast<std::size_t>(RyR_akap);
    this->mapped_data_["ka_ryr"] = static_cast<std::size_t>(ka_ryr);
    this->mapped_data_["kp_ryr"] = static_cast<std::size_t>(kp_ryr);
    this->mapped_data_["Ka_ryr"] = static_cast<std::size_t>(Ka_ryr);
    this->mapped_data_["Kp_ryr"] = static_cast<std::size_t>(Kp_ryr);
    this->mapped_data_["Mr"] = static_cast<std::size_t>(Mr);
    this->mapped_data_["Lr"] = static_cast<std::size_t>(Lr);
    this->mapped_data_["Kr"] = static_cast<std::size_t>(Kr);
    this->mapped_data_["akap_sig_RyRf_sum"] = static_cast<std::size_t>(akap_sig_RyRf_sum);
    this->mapped_data_["RyRf"] = static_cast<std::size_t>(RyRf);
    this->mapped_data_["ICaL_tot"] = static_cast<std::size_t>(ICaL_tot);
    this->mapped_data_["ICaL_akap"] = static_cast<std::size_t>(ICaL_akap);
    this->mapped_data_["ka_ical"] = static_cast<std::size_t>(ka_ical);
    this->mapped_data_["kp_ical"] = static_cast<std::size_t>(kp_ical);
    this->mapped_data_["Ka_ical"] = static_cast<std::size_t>(Ka_ical);
    this->mapped_data_["Kp_ical"] = static_cast<std::size_t>(Kp_ical);
    this->mapped_data_["Mi"] = static_cast<std::size_t>(Mi);
    this->mapped_data_["Li"] = static_cast<std::size_t>(Li);
    this->mapped_data_["Ki"] = static_cast<std::size_t>(Ki);
    this->mapped_data_["akap_sig_ICaLf_sum"] = static_cast<std::size_t>(akap_sig_ICaLf_sum);
    this->mapped_data_["ICaLf"] = static_cast<std::size_t>(ICaLf);
    this->mapped_data_["akap_sig_PP1f_cav_b"] = static_cast<std::size_t>(akap_sig_PP1f_cav_b);
    this->mapped_data_["akap_sig_PP1f_cav_c"] = static_cast<std::size_t>(akap_sig_PP1f_cav_c);
    this->mapped_data_["akap_sig_PP1f_cav_d"] = static_cast<std::size_t>(akap_sig_PP1f_cav_d);
    this->mapped_data_["akap_sig_PP1f_cav_rr"] = static_cast<std::size_t>(akap_sig_PP1f_cav_rr);
    this->mapped_data_["akap_sig_PP1f_cav_yi"] = static_cast<std::size_t>(akap_sig_PP1f_cav_yi);
    this->mapped_data_["akap_sig_PP1f_cav_yr"] = static_cast<std::size_t>(akap_sig_PP1f_cav_yr);
    this->mapped_data_["akap_sig_PP1f_cav_mag"] = static_cast<std::size_t>(akap_sig_PP1f_cav_mag);
    this->mapped_data_["akap_sig_PP1f_cav_arg"] = static_cast<std::size_t>(akap_sig_PP1f_cav_arg);
    this->mapped_data_["akap_sig_PP1f_cav_x"] = static_cast<std::size_t>(akap_sig_PP1f_cav_x);
    this->mapped_data_["PP1f_cav"] = static_cast<std::size_t>(PP1f_cav);
    this->mapped_data_["akap_sig_PKAf_d"] = static_cast<std::size_t>(akap_sig_PKAf_d);
    this->mapped_data_["akap_sig_PKAf_b"] = static_cast<std::size_t>(akap_sig_PKAf_b);
    this->mapped_data_["akap_sig_PKAf_c"] = static_cast<std::size_t>(akap_sig_PKAf_c);
    this->mapped_data_["akap_sig_PKAf_rr"] = static_cast<std::size_t>(akap_sig_PKAf_rr);
    this->mapped_data_["akap_sig_PKAf_yr"] = static_cast<std::size_t>(akap_sig_PKAf_yr);
    this->mapped_data_["akap_sig_PKAf_yi"] = static_cast<std::size_t>(akap_sig_PKAf_yi);
    this->mapped_data_["akap_sig_PKAf_mag"] = static_cast<std::size_t>(akap_sig_PKAf_mag);
    this->mapped_data_["akap_sig_PKAf_arg"] = static_cast<std::size_t>(akap_sig_PKAf_arg);
    this->mapped_data_["akap_sig_PKAf_x"] = static_cast<std::size_t>(akap_sig_PKAf_x);
    this->mapped_data_["akap_sig_PKAf"] = static_cast<std::size_t>(akap_sig_PKAf);
    this->mapped_data_["RyR_akapf"] = static_cast<std::size_t>(RyR_akapf);
    this->mapped_data_["RyR_arn"] = static_cast<std::size_t>(RyR_arn);
    this->mapped_data_["RyR_arp"] = static_cast<std::size_t>(RyR_arp);
    this->mapped_data_["ICaL_akapf"] = static_cast<std::size_t>(ICaL_akapf);
    this->mapped_data_["ICaL_arn"] = static_cast<std::size_t>(ICaL_arn);
    this->mapped_data_["ICaL_arp"] = static_cast<std::size_t>(ICaL_arp);
    this->mapped_data_["beta_0"] = static_cast<std::size_t>(beta_0);
    this->mapped_data_["irel_fhat_ratio"] = static_cast<std::size_t>(irel_fhat_ratio);
    this->mapped_data_["ical_f_hat_ratio"] = static_cast<std::size_t>(ical_f_hat_ratio);
    this->mapped_data_["iks_f_hat_ratio"] = static_cast<std::size_t>(iks_f_hat_ratio);
    this->mapped_data_["Whole_cell_PP1"] = static_cast<std::size_t>(Whole_cell_PP1);
    this->mapped_data_["PP1block"] = static_cast<std::size_t>(PP1block);
    this->mapped_data_["fICaLP"] = static_cast<std::size_t>(fICaLP);
    this->mapped_data_["fIKsP"] = static_cast<std::size_t>(fIKsP);
    this->mapped_data_["fPLBP"] = static_cast<std::size_t>(fPLBP);
    this->mapped_data_["fTnIP"] = static_cast<std::size_t>(fTnIP);
    this->mapped_data_["fINaP"] = static_cast<std::size_t>(fINaP);
    this->mapped_data_["fINaKP"] = static_cast<std::size_t>(fINaKP);
    this->mapped_data_["fRyRP"] = static_cast<std::size_t>(fRyRP);
    this->mapped_data_["fIKurP"] = static_cast<std::size_t>(fIKurP);

    // Set currents mapping.
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["INaL"] = static_cast<std::size_t>(INaL);
    this->mapped_data_["Ito"] = static_cast<std::size_t>(Ito);
    this->mapped_data_["ICaL"] = static_cast<std::size_t>(ICaL);
    this->mapped_data_["ICaNa"] = static_cast<std::size_t>(ICaNa);
    this->mapped_data_["ICaK"] = static_cast<std::size_t>(ICaK);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);
    this->mapped_data_["INaCa_i"] = static_cast<std::size_t>(INaCa_i);
    this->mapped_data_["INaCa_ss"] = static_cast<std::size_t>(INaCa_ss);
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["IKb"] = static_cast<std::size_t>(IKb);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);
    this->mapped_data_["IpCa"] = static_cast<std::size_t>(IpCa);
    this->mapped_data_["Iion"] = static_cast<std::size_t>(Gng20Cur::Iion);
}


Gong2020::Gong2020()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Gong2020;
    this->dt_stable_ = 0.002;
    this->var_.resize(116, 0.);
    this->prm_.resize(353, 0.);
    this->cur_.resize(17, 0.);
    this->block_coeff_.resize(16, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


Gong2020::~Gong2020()
{}


void Gong2020::Initialize(CellType cell_type)
{
    using namespace Gng20Var;
    using namespace Gng20Prm;

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(116, 0.);
    this->prm_.clear();           this->prm_.resize(353, 0.);
    this->cur_.clear();           this->cur_.resize(17, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(16, 0.);

    // Set initial values for the variables.
    this->var_[v] = -87.;
    this->var_[dvdt] = 0.;
    this->var_[nai] = 7.;
    this->var_[nass] = 7.;
    this->var_[ki] = 145.;
    this->var_[kss] = 145.;
    this->var_[cai] = 1.0e-4;
    this->var_[cass] = 1.0e-4;
    this->var_[cansr] = 1.2;
    this->var_[cajsr] = 1.2;
    this->var_[m] = 0.;
    this->var_[hf] = 1.;
    this->var_[hs] = 1.;
    this->var_[j] = 1.;
    this->var_[hsp] = 1.;
    this->var_[jp] = 1.;
    this->var_[mL] = 0.;
    this->var_[hL] = 1.;
    this->var_[hLp] = 1.;
    this->var_[a] = 0.;
    this->var_[iF] = 1.;
    this->var_[iS] = 1.;
    this->var_[ap] = 0.;
    this->var_[iFp] = 1.;
    this->var_[iSp] = 1.;
    this->var_[d] = 0.;
    this->var_[ff] = 1.;
    this->var_[fs] = 1.;
    this->var_[fcaf] = 1.;
    this->var_[fcas] = 1.;
    this->var_[jca] = 1.;
    this->var_[nca] = 0.;
    this->var_[ffp] = 1.;
    this->var_[fcafp] = 1.;
    this->var_[xrf] = 0.;
    this->var_[xrs] = 0.;
    this->var_[xs1] = 0.;
    this->var_[xs2] = 0.;
    this->var_[xk1] = 1.;
    this->var_[Jrelnp] = 0.;
    this->var_[Jrelp] = 0.;
    this->var_[CaMKt] = 0.;
    this->var_[hPf] = 1.;
    this->var_[hPs] = 1.;
    this->var_[jP] = 1.;
    this->var_[hBPf] = 1.;
    this->var_[hBPs] = 1.;
    this->var_[jBP] = 1.;
    this->var_[dP] = 0.;
    this->var_[fPf] = 1.;
    this->var_[fPs] = 1.;
    this->var_[fcaPf] = 1.;
    this->var_[fcaPs] = 1.;
    this->var_[fBPf] = 1.;
    this->var_[fcaBPf] = 1.;
    this->var_[xs1P] = 0.;
    this->var_[xs2P] = 0.;
    this->var_[JrelP] = 0.;
    this->var_[JrelBP] = 0.;
    this->var_[beta_cav_Gs_aGTP] = 0.00685041638458665;
    this->var_[beta_eca_Gs_aGTP] = 0.0184627603007976;
    this->var_[beta_cyt_Gs_aGTP] = 0.000731420577213056;
    this->var_[beta_cav_Gs_bg] = 0.00745773247314215;
    this->var_[beta_eca_Gs_bg] = 0.0191017987408719;
    this->var_[beta_cyt_Gs_bg] = 0.00115141243826747;
    this->var_[beta_cav_Gs_aGDP] = 0.000607316088556676;
    this->var_[beta_eca_Gs_aGDP] = 0.000639038440072507;
    this->var_[beta_cyt_Gs_aGDP] = 0.000419991861054322;
    this->var_[cAMP_cav] = 0.347102959606005;
    this->var_[cAMP_eca] = 9.62359241535767;
    this->var_[cAMP_cyt] = 0.474081735738211;
    this->var_[beta_cav_Rb1_pka_tot] = 0.0149041813757831;
    this->var_[beta_eca_Rb1_pka_tot] = 0.203016833596288;
    this->var_[beta_cyt_Rb1_pka_tot] = 0.00944463350378086;
    this->var_[beta_cav_Rb1_grk_tot] = 2.49592854373432e-10;
    this->var_[beta_eca_Rb1_grk_tot] = 1.18055788874765e-09;
    this->var_[beta_cyt_Rb1_grk_tot] = 7.07824478944671e-11;
    this->var_[pka_cav_ARC] = 0.0904820284659604;
    this->var_[pka_cav_A2RC] = 0.00276490711096605;
    this->var_[pka_cav_A2R] = 0.225475702283053;
    this->var_[pka_cav_C] = 0.0326565916584703;
    this->var_[pka_cav_PKIC] = 0.192819110624505;
    this->var_[pka_eca_ARC] = 0.205444874210056;
    this->var_[pka_eca_A2RC] = 0.174057375932567;
    this->var_[pka_eca_A2R] = 0.817161796756964;
    this->var_[pka_eca_C] = 0.567249910261073;
    this->var_[pka_eca_PKIC] = 0.249911886495890;
    this->var_[pka_cyt_ARC] = 0.0646928309115710;
    this->var_[pka_cyt_A2RC] = 0.0664997605558791;
    this->var_[pka_cyt_A2R] = 0.489063888619456;
    this->var_[pka_cyt_C] = 0.362113356111496;
    this->var_[pka_cyt_PKIC] = 0.126950532507959;
    this->var_[PDE3_P_cav] = 0.0236821659448037;
    this->var_[PDE3_P_cyt] = 0.0128402905095188;
    this->var_[PDE4_P_cav] = 0.00637363047239019;
    this->var_[PDE4_P_eca] = 4.29171113639322e-05;
    this->var_[PDE4_P_cyt] = 0.00917039986149184;
    this->var_[inhib1_p] = 0.0282662056977524;
    this->var_[ICaLp] = 0.000673713947839317;
    this->var_[IKsp] = 0.000765988420110534;
    this->var_[iup_f_plb] = 0.592167467082831;
    this->var_[f_tni] = 0.673518785672382;
    this->var_[ina_f_ina] = 0.239479458960528;
    this->var_[f_inak] = 0.126345311579566;
    this->var_[RyRp] = 0.00410693810508171;
    this->var_[f_ikur] = 0.0589379755147718;
    this->var_[beta_cav_Rb2_pka_tot] = 0.0275455839709412;
    this->var_[beta_cav_Rb2_grk_tot] = 8.91799633266019e-10;
    this->var_[beta_cav_Gi_aGTP] = 0.00159632196638178;
    this->var_[beta_cav_Gi_bg] = 0.00209911481235842;
    this->var_[beta_cav_Gi_aGDP] = 0.000502792845976641;
    this->var_[beta_eca_Rb2_pka_tot] = 0.0110248953370551;
    this->var_[beta_eca_Rb2_grk_tot] = 1.13428924662652e-10;
    this->var_[beta_eca_Gi_aGTP] = 0.000364315164237569;
    this->var_[beta_eca_Gi_bg] = 0.000705656306851923;
    this->var_[beta_eca_Gi_aGDP] = 0.000341341142614041;

    // Set initial values for the parameters.
    if (cell_type == CellType::endo) this->prm_[celltype] = 0.;
    if (cell_type == CellType::epi) this->prm_[celltype] = 1.;
    if (cell_type == CellType::mid) this->prm_[celltype] = 2.;
    this->prm_[signaling] = 1.;
    this->prm_[electrophys] = 1.;
    this->prm_[cipa] = 1.;
    this->prm_[GKS_mul] = 8.;
    this->prm_[VP_d] = 12.;
    this->prm_[VP_f] = 8.;
    this->prm_[GCal_mul] = 1.9;

    if (static_cast<int>(this->prm_[cipa]) == 1) {
        this->prm_[CiPA_scale_IKr] = 1.119;
        this->prm_[CiPA_scale_IKs] = 1.648;
        this->prm_[CiPA_scale_IK1] = 1.414;
        this->prm_[CiPA_scale_ICaL] = 1.018;
        this->prm_[CiPA_scale_INaL] = 2.274;
    } else {
        this->prm_[CiPA_scale_IKr] = 1.;
        this->prm_[CiPA_scale_IKs] = 1.;
        this->prm_[CiPA_scale_IK1] = 1.;
        this->prm_[CiPA_scale_ICaL] = 1.;
        this->prm_[CiPA_scale_INaL] = 1.;
    }

    this->prm_[nao] = 140.0;
    this->prm_[cao] = 1.8;
    this->prm_[ko] = 5.4;
    this->prm_[R] = 8314.0;
    this->prm_[T] = 310.0;
    this->prm_[F] = 96485.0;
    this->prm_[L] = 0.01;
    this->prm_[rad] = 0.0011;
    this->prm_[vcell] = 1000*3.14*this->prm_[rad]*this->prm_[rad]*this->prm_[L];
    this->prm_[Ageo] = 2.*3.14*this->prm_[rad]*this->prm_[rad]+2.*3.14*this->prm_[rad]*this->prm_[L];
    this->prm_[Acap] = 2.*this->prm_[Ageo];
    this->prm_[vmyo] = 0.68*this->prm_[vcell];
    this->prm_[vnsr] = 0.0552*this->prm_[vcell];
    this->prm_[vjsr] = 0.0048*this->prm_[vcell];
    this->prm_[vss] = 0.02*this->prm_[vcell];
    this->prm_[aCaMK] = 0.05;
    this->prm_[bCaMK] = 0.00068;
    this->prm_[CaMKo] = 0.05;
    this->prm_[PKNa] = 0.01833;
    this->prm_[GNa] = 75.;
    this->prm_[GNaL] = 0.0075 * this->prm_[CiPA_scale_INaL];
    this->prm_[Gto] = 0.02;
    this->prm_[GKr] = 0.046 * this->prm_[CiPA_scale_IKr];
    this->prm_[GKs] = 0.0034 * this->prm_[CiPA_scale_IKs];
    this->prm_[GK1] = 0.1908 * this->prm_[CiPA_scale_IK1];
    this->prm_[Gncx] = 0.0008;
    this->prm_[GKb] = 0.003;
    this->prm_[GpCa] = 0.0005;
    this->prm_[PCa] = 0.0001 * this->prm_[CiPA_scale_ICaL];
    this->prm_[Pnak] = 30.;

    if (static_cast<int>(this->prm_[celltype]) == 1) {
        this->prm_[GNaL] *= 0.6;
        this->prm_[Gto] *= 4.0;
        this->prm_[GKr] *= 1.3;
        this->prm_[GKs] *= 1.4;
        this->prm_[GK1] *= 1.2;
        this->prm_[Gncx] *= 1.1;
        this->prm_[GKb] *= 0.6;
        this->prm_[PCa] *= 1.2;
        this->prm_[Pnak] *= 0.9;
    } else if (static_cast<int>(this->prm_[celltype]) == 2) {
        this->prm_[Gto] *= 4.0;
        this->prm_[GKr] *= 0.8;
        this->prm_[GK1] *= 1.3;
        this->prm_[Gncx] *= 1.4;
        this->prm_[PCa] *= 2.5;
        this->prm_[Pnak] *= 0.7;
    }

    this->prm_[PNab] = 3.75e-10;
    this->prm_[PCab] = 2.5e-8;
    this->prm_[KmCaMK] = 0.15;
    this->prm_[KmCaM] = 0.0015;
    this->prm_[SERCA_total]  = 1.;
    this->prm_[RyR_total]  = 1.;
    this->prm_[Leak_total]  = 1.;
    this->prm_[Trans_total]  = 1.;
    this->prm_[iso_L] = 0.;
    this->prm_[ibmx] = 0.;
    this->prm_[phys_T] = 310.;
    this->prm_[F] = 96487.;
    this->prm_[phys_R] = 8314.;
    this->prm_[length] = 0.01;
    this->prm_[cell_pi] =  3.14159265358979312e+00;
    this->prm_[radius] = 0.0011;
    this->prm_[geoArea] = 2.0 * this->prm_[cell_pi] * this->prm_[radius] * (this->prm_[radius] + this->prm_[length]);
    this->prm_[volume] = 1000.0 * this->prm_[cell_pi] * this->prm_[radius] * this->prm_[radius] * this->prm_[length];
    this->prm_[v_cav] = 0.02 * this->prm_[volume];
    this->prm_[v_eca] = 0.04 * this->prm_[volume];
    this->prm_[v_cyt] = this->prm_[volume] * 0.678;
    this->prm_[vr_cav] = this->prm_[volume] / this->prm_[v_cav];
    this->prm_[vr_eca] = this->prm_[volume] / this->prm_[v_eca];
    this->prm_[vr_cyt] = this->prm_[volume] / this->prm_[v_cyt];
    this->prm_[j_cav_eca] = 5e-15 * 1000000.0;
    this->prm_[j_cav_cyt] = 7.5e-14 * 1000000.0;
    this->prm_[j_eca_cyt] = 9e-15 * 1000000.0;
    this->prm_[PKA_tot] = 0.5;
    this->prm_[f_cav] = 0.0388;
    this->prm_[f_eca] = 0.1;
    this->prm_[f_cyt] = 1.0 - this->prm_[f_cav] - this->prm_[f_eca];
    this->prm_[PKA_eca] = this->prm_[f_eca] * this->prm_[PKA_tot] * this->prm_[vr_eca];
    this->prm_[PKA_cav] = this->prm_[f_cav] * this->prm_[PKA_tot] * this->prm_[vr_cav];
    this->prm_[PKA_cyt] = this->prm_[f_cyt] * this->prm_[PKA_tot] * this->prm_[vr_cyt];
    this->prm_[PKI_tot] = 0.2 * this->prm_[PKA_tot];
    this->prm_[f_pki_cav] = this->prm_[f_cav];
    this->prm_[f_pki_eca] = this->prm_[f_eca];
    this->prm_[f_pki_cyt] = 1.0 - this->prm_[f_pki_cav] - this->prm_[f_pki_eca];
    this->prm_[PKI_cav] = this->prm_[f_pki_cav] * this->prm_[PKI_tot] * this->prm_[vr_cav];
    this->prm_[PKI_eca] = this->prm_[f_pki_eca] * this->prm_[PKI_tot] * this->prm_[vr_eca];
    this->prm_[PKI_cyt] = this->prm_[f_pki_cyt] * this->prm_[PKI_tot] * this->prm_[vr_cyt];
    this->prm_[f_pki] = 50.0;
    this->prm_[K_pki] = 0.01 / 50.0;
    this->prm_[b_pki] = this->prm_[f_pki] * this->prm_[K_pki];
    this->prm_[pka_cav_f1] = 100.0;
    this->prm_[pka_cav_f2] = 100.0;
    this->prm_[pka_cav_f3] = 100.0;
    this->prm_[pka_cav_K1] = 2.4984;
    this->prm_[pka_cav_K2] = 11.359;
    this->prm_[pka_cav_K3] = 0.3755;
    this->prm_[pka_cav_b1] = this->prm_[pka_cav_f1] * this->prm_[pka_cav_K1];
    this->prm_[pka_cav_b2] = this->prm_[pka_cav_f2] * this->prm_[pka_cav_K2];
    this->prm_[pka_cav_b3] = this->prm_[pka_cav_f3] * this->prm_[pka_cav_K3];
    this->prm_[pka_eca_f1] = this->prm_[pka_cav_f1];
    this->prm_[pka_eca_f2] = this->prm_[pka_cav_f2];
    this->prm_[pka_eca_f3] = this->prm_[pka_cav_f3];
    this->prm_[pka_eca_K1] = this->prm_[pka_cav_K1];
    this->prm_[pka_eca_K2] = this->prm_[pka_cav_K2];
    this->prm_[pka_eca_K3] = this->prm_[pka_cav_K3];
    this->prm_[pka_eca_b1] = this->prm_[pka_eca_f1] * this->prm_[pka_eca_K1];
    this->prm_[pka_eca_b2] = this->prm_[pka_eca_f2] * this->prm_[pka_eca_K2];
    this->prm_[pka_eca_b3] = this->prm_[pka_eca_f3] * this->prm_[pka_eca_K3];
    this->prm_[pka_cyt_f1] = this->prm_[pka_cav_f1];
    this->prm_[pka_cyt_f2] = this->prm_[pka_cav_f2];
    this->prm_[pka_cyt_f3] = this->prm_[pka_cav_f3];
    this->prm_[pka_cyt_K1] = 0.1088;
    this->prm_[pka_cyt_K2] = 0.4612;
    this->prm_[pka_cyt_K3] = 0.3755;
    this->prm_[pka_cyt_b1] = this->prm_[pka_cyt_f1] * this->prm_[pka_cyt_K1];
    this->prm_[pka_cyt_b2] = this->prm_[pka_cyt_f2] * this->prm_[pka_cyt_K2];
    this->prm_[pka_cyt_b3] = this->prm_[pka_cyt_f3] * this->prm_[pka_cyt_K3];
    this->prm_[f] = 0.3;
    this->prm_[kp] = 0.010145;
    this->prm_[kdp] = 0.0035731;
    this->prm_[Kp] = 0.001469;
    this->prm_[Kdp] =  1.95259999999999991e-05;
    this->prm_[PP1_eca] = 0.1;
    this->prm_[PP1_cav] = 0.25;
    this->prm_[PP1_cyt] = 0.2;
    this->prm_[pp1_K] = 0.001;
    this->prm_[inhib1_tot] = this->prm_[f] / (1.0 - this->prm_[f]) * this->prm_[pp1_K] + this->prm_[f] * this->prm_[PP1_cyt];
    this->prm_[PP2A] = 1.0;
    this->prm_[PDE2_tot] = 0.029268;
    this->prm_[f_pde2_cav] = 0.16957;
    this->prm_[f_pde2_eca] =  2.12570000000000006e-04;
    this->prm_[f_pde2_cyt] = 1.0 - this->prm_[f_pde2_cav] - this->prm_[f_pde2_eca];
    this->prm_[f_pde2_part] = this->prm_[f_pde2_cav] + this->prm_[f_pde2_eca];
    this->prm_[f_pde_part] = 0.2;
    this->prm_[f_pde4_part] = 0.125;
    this->prm_[f_pde4_cav] = 0.12481;
    this->prm_[f_pde4_eca] = this->prm_[f_pde4_part] - this->prm_[f_pde4_cav];
    this->prm_[f_pde4_cyt] = 1.0 - this->prm_[f_pde4_part];
    this->prm_[kPDE2] = 20.0;
    this->prm_[kPDE3] = 2.5;
    this->prm_[kPDE4] = 4.0;
    this->prm_[KmPDE2] = 50.0;
    this->prm_[KmPDE3] = 0.8;
    this->prm_[KmPDE4] = 1.4;
    this->prm_[KmIbmxPde2] = 21.58;
    this->prm_[h_ibmx_pde2] = 1.167;
    this->prm_[h_ibmx_pde3] = 0.7629;
    this->prm_[h_ibmx_pde4] = 0.9024;
    this->prm_[KmIbmxPde3] = 2.642;
    this->prm_[KmIbmxPde4] = 11.89;
    this->prm_[KPDEp] = 0.52218;
    this->prm_[delta_k_pde34] = 3.0;
    this->prm_[ff_pde3_cyt] = 0.35;
    this->prm_[kfPDEp] = 0.0196;
    this->prm_[r_pde34_frac] = 3.71;
    this->prm_[r_pde3_cyt] = this->prm_[ff_pde3_cyt] / (1.0 - this->prm_[ff_pde3_cyt]);
    this->prm_[kbPDEp] = this->prm_[KPDEp] * this->prm_[kfPDEp];
    this->prm_[pde_PDE3_tot_alpha] = this->prm_[r_pde3_cyt] * (this->prm_[f_pde4_part] * (1.0 + this->prm_[r_pde34_frac] - this->prm_[r_pde34_frac] * this->prm_[f_pde2_part] - this->prm_[f_pde_part]) + this->prm_[f_pde2_part] * (this->prm_[f_pde_part] - 1.0)) + this->prm_[r_pde34_frac] * this->prm_[f_pde4_part] * (this->prm_[f_pde_part] - this->prm_[f_pde2_part]);
    this->prm_[pde_PDE3_tot_beta] = this->prm_[f_pde4_part] * (1.0 + this->prm_[r_pde34_frac] + this->prm_[f_pde_part] * (this->prm_[r_pde3_cyt] - this->prm_[r_pde34_frac])) - this->prm_[f_pde_part] * (1.0 + this->prm_[r_pde3_cyt]);
    this->prm_[PDE3_tot] = this->prm_[pde_PDE3_tot_alpha] / this->prm_[pde_PDE3_tot_beta] * this->prm_[PDE2_tot];
    this->prm_[PDE4_tot] = ((this->prm_[f_pde_part] - this->prm_[f_pde2_part]) * this->prm_[PDE2_tot] + this->prm_[f_pde_part] * this->prm_[PDE3_tot]) / ((1.0 + this->prm_[r_pde34_frac]) * this->prm_[f_pde4_part] - this->prm_[f_pde_part]);
    this->prm_[ibmx_h2] = std::pow(this->prm_[ibmx], this->prm_[h_ibmx_pde2]);
    this->prm_[ibmx_h3] = std::pow(this->prm_[ibmx], this->prm_[h_ibmx_pde3]);
    this->prm_[ibmx_h4] = std::pow(this->prm_[ibmx], this->prm_[h_ibmx_pde4]);
    this->prm_[ibmx2] = (1.0 - this->prm_[ibmx_h2] / (this->prm_[KmIbmxPde2] + this->prm_[ibmx_h2])) * this->prm_[PDE2_tot];
    this->prm_[ibmx3] = (1.0 - this->prm_[ibmx_h3] / (this->prm_[KmIbmxPde3] + this->prm_[ibmx_h3])) * this->prm_[PDE3_tot];
    this->prm_[ibmx4] = (1.0 - this->prm_[ibmx_h4] / (this->prm_[KmIbmxPde4] + this->prm_[ibmx_h4])) * this->prm_[PDE4_tot];
    this->prm_[f_pde3_cav] = this->prm_[r_pde34_frac] * this->prm_[f_pde4_part] * this->prm_[PDE4_tot] / this->prm_[PDE3_tot];
    this->prm_[f_pde3_cyt] = 1.0 - this->prm_[f_pde3_cav];
    this->prm_[PDE2_cav] = this->prm_[ibmx2] * this->prm_[f_pde2_cav] * this->prm_[vr_cav];
    this->prm_[PDE2_eca] = this->prm_[ibmx2] * this->prm_[f_pde2_eca] * this->prm_[vr_eca];
    this->prm_[PDE2_cyt] = this->prm_[ibmx2] * this->prm_[f_pde2_cyt] * this->prm_[vr_cyt];
    this->prm_[PDE3_tot] = this->prm_[pde_PDE3_tot_alpha] / this->prm_[pde_PDE3_tot_beta] * this->prm_[PDE2_tot];
    this->prm_[PDE3_tot] = this->prm_[pde_PDE3_tot_alpha] / this->prm_[pde_PDE3_tot_beta] * this->prm_[PDE2_tot];
    this->prm_[PDE3_cav] = this->prm_[ibmx3] * this->prm_[f_pde3_cav] * this->prm_[vr_cav];
    this->prm_[PDE3_cyt] = this->prm_[ibmx3] * this->prm_[f_pde3_cyt] * this->prm_[vr_cyt];
    this->prm_[PDE4_cav] = this->prm_[ibmx4] * this->prm_[f_pde4_cav] * this->prm_[vr_cav];
    this->prm_[PDE4_eca] = this->prm_[ibmx4] * this->prm_[f_pde4_eca] * this->prm_[vr_eca];
    this->prm_[PDE4_cyt] = this->prm_[ibmx4] * this->prm_[f_pde4_cyt] * this->prm_[vr_cyt];
    this->prm_[beta_R_b1_tot] = 0.85 * 0.025;
    this->prm_[beta_R_b2_tot] = 0.15 * 0.025;
    this->prm_[Gs_tot] = 224.0 * this->prm_[beta_R_b1_tot];
    this->prm_[Gi_tot] = 0.5;
    this->prm_[f_Gs_eca] = 0.5664;
    this->prm_[f_Gs_cav] = 0.0011071;
    this->prm_[f_Gs_cyt] = 1.0 - this->prm_[f_Gs_cav] - this->prm_[f_Gs_eca];
    this->prm_[f_Gi_cav] = 0.85;
    this->prm_[f_Gi_eca] = 1.0 - this->prm_[f_Gi_cav];
    this->prm_[f_Rb1_cav] = 0.081161;
    this->prm_[f_Rb1_eca] = 0.48744;
    this->prm_[f_Rb1_cyt] = 1.0 - this->prm_[f_Rb1_cav] - this->prm_[f_Rb1_eca];
    this->prm_[f_Rb2_cav] = 0.85;
    this->prm_[f_Rb2_eca] = 1.0 - this->prm_[f_Rb2_cav];
    this->prm_[k_b1_l] = 0.567;
    this->prm_[k_b1_c] = 2.449;
    this->prm_[k_b1_h] = 0.062;
    this->prm_[k_b2_n] = 1.053;
    this->prm_[k_b2_h] = 0.012;
    this->prm_[k_b2_a] = 1.6655;
    this->prm_[k_b2_c] = 1.8463;
    this->prm_[k_b2_l] = 1.053;
    this->prm_[k_b2_f] = 0.1;
    this->prm_[k_act1_Gs] = 4.9054;
    this->prm_[k_act2_Gs] = 0.25945;
    this->prm_[k_act1_Gi] = 4.0;
    this->prm_[k_act2_Gi] = 0.05;
    this->prm_[k_hydr_Gs] = 0.8;
    this->prm_[k_hydr_Gi] = this->prm_[k_hydr_Gs];
    this->prm_[k_reas_Gs] = 1210.0;
    this->prm_[k_reas_Gi] = this->prm_[k_reas_Gs];
    this->prm_[rate_bds] = 0.35;
    this->prm_[k_grk_dp] = this->prm_[rate_bds] * 0.0009833;
    this->prm_[k_grk_p] = this->prm_[rate_bds] * 0.00133;
    this->prm_[k_pka_p] = this->prm_[rate_bds] * 0.0065;
    this->prm_[k_pka_dp] = 0.15629 * this->prm_[k_pka_p];
    this->prm_[beta_cav_GRK] = 1.0;
    this->prm_[beta_cav_R_b1_tot] = this->prm_[f_Rb1_cav] * this->prm_[beta_R_b1_tot] * this->prm_[vr_cav];
    this->prm_[beta_cav_k_GsAct_b2] = 1.0;
    this->prm_[beta_cav_R_b2_tot] = this->prm_[f_Rb2_cav] * this->prm_[beta_R_b2_tot] * this->prm_[vr_cav];
    this->prm_[beta_cav_Rb2_pka_f_a] = (this->prm_[k_b2_f] + this->prm_[iso_L]) * (this->prm_[k_b2_n] + this->prm_[iso_L]) / this->prm_[k_b2_n];
    this->prm_[beta_cav_Gs_f_c22] = this->prm_[k_b2_c] * this->prm_[k_b2_h] * this->prm_[k_b1_l] * (this->prm_[k_b1_h] + this->prm_[iso_L]) * (this->prm_[k_b2_l] + this->prm_[iso_L]);
    this->prm_[beta_cav_Gs_f_a] = this->prm_[k_b1_l] * this->prm_[k_b2_l] * (this->prm_[k_b1_h] + this->prm_[iso_L]) * (this->prm_[k_b2_h] + this->prm_[iso_L]);
    this->prm_[beta_cav_Gs_f_c33] = this->prm_[k_b1_c] * this->prm_[k_b2_c] * this->prm_[k_b1_h] * this->prm_[k_b2_h] * (this->prm_[k_b1_l] + this->prm_[iso_L]) * (this->prm_[k_b2_l] + this->prm_[iso_L]);
    this->prm_[beta_cav_Gs_f_c11] = this->prm_[k_b1_c] * this->prm_[k_b1_h] * this->prm_[k_b2_l] * (this->prm_[k_b2_h] + this->prm_[iso_L]) * (this->prm_[k_b1_l] + this->prm_[iso_L]);
    this->prm_[beta_eca_GRK] = 1.0;
    this->prm_[beta_eca_k_GsAct_b2] = 1.0;
    this->prm_[beta_eca_R_b2_tot] = this->prm_[f_Rb2_eca] * this->prm_[beta_R_b2_tot] * this->prm_[vr_eca];
    this->prm_[beta_eca_R_b1_tot] = this->prm_[f_Rb1_eca] * this->prm_[beta_R_b1_tot] * this->prm_[vr_eca];
    this->prm_[beta_eca_Rb2_pka_f_a] = (this->prm_[k_b2_f] + this->prm_[iso_L]) * (this->prm_[k_b2_n] + this->prm_[iso_L]) / this->prm_[k_b2_n];
    this->prm_[beta_eca_Gs_f_c11] = this->prm_[k_b1_c] * this->prm_[k_b1_h] * this->prm_[k_b2_l] * (this->prm_[k_b2_h] + this->prm_[iso_L]) * (this->prm_[k_b1_l] + this->prm_[iso_L]);
    this->prm_[beta_eca_Gs_f_c33] = this->prm_[k_b1_c] * this->prm_[k_b2_c] * this->prm_[k_b1_h] * this->prm_[k_b2_h] * (this->prm_[k_b1_l] + this->prm_[iso_L]) * (this->prm_[k_b2_l] + this->prm_[iso_L]);
    this->prm_[beta_eca_Gs_f_a] = this->prm_[k_b1_l] * this->prm_[k_b2_l] * (this->prm_[k_b1_h] + this->prm_[iso_L]) * (this->prm_[k_b2_h] + this->prm_[iso_L]);
    this->prm_[beta_eca_Gs_f_c22] = this->prm_[k_b2_c] * this->prm_[k_b2_h] * this->prm_[k_b1_l] * (this->prm_[k_b1_h] + this->prm_[iso_L]) * (this->prm_[k_b2_l] + this->prm_[iso_L]);
    this->prm_[beta_cyt_R_b1_tot] = this->prm_[f_Rb1_cyt] * this->prm_[beta_R_b1_tot] * this->prm_[vr_cyt];
    this->prm_[beta_cyt_GRK] = 1.0;
    this->prm_[beta_cyt_Rb1_np_f_a] = (this->prm_[k_b1_h] + this->prm_[iso_L]) * (this->prm_[k_b1_l] + this->prm_[iso_L]) / this->prm_[k_b1_l];
    this->prm_[ATP] = 5000.0;
    this->prm_[KmATP] = 315.0;
    this->prm_[AC_tot] = 3.0 * this->prm_[beta_R_b1_tot];
    this->prm_[f_AC47_eca] = 0.16479;
    this->prm_[f_AC56_cav] = 0.087459;
    this->prm_[f_AC56_AC47] = 1.0 / (1.0 + 0.35);
    this->prm_[hGsAC47] = 1.0043;
    this->prm_[hGsAC56] = 1.3574;
    this->prm_[hGsGiAC56] = 0.6623;
    this->prm_[KmGsAC47] = 0.031544;
    this->prm_[KmGsAC56] = 0.0852;
    this->prm_[KmGiAC56] = 0.0465;
    this->prm_[KmGsGiAC56] = 0.4824;
    this->prm_[basalAC47] = 0.03135;
    this->prm_[basalAC56] = 0.037696;
    this->prm_[afAC47] = 3.3757;
    this->prm_[afAC56] = 41.32;
    this->prm_[vGsGiAC56] = 0.8569;
    this->prm_[AC47_cyt] = (1.0 - this->prm_[f_AC47_eca]) * (1.0 - this->prm_[f_AC56_AC47]) * this->prm_[AC_tot] * this->prm_[vr_cyt];
    this->prm_[AC56_cav] = this->prm_[f_AC56_cav] * this->prm_[f_AC56_AC47] * this->prm_[AC_tot] * this->prm_[vr_cav];
    this->prm_[fATP] = this->prm_[ATP] / (this->prm_[KmATP] + this->prm_[ATP]);
    this->prm_[AC47_eca] = this->prm_[f_AC47_eca] * (1.0 - this->prm_[f_AC56_AC47]) * this->prm_[AC_tot] * this->prm_[vr_eca];
    this->prm_[AC56_cyt] = (1.0 - this->prm_[f_AC56_cav]) * this->prm_[f_AC56_AC47] * this->prm_[AC_tot] * this->prm_[vr_cyt];
    this->prm_[ka_inak] = 0.015265;
    this->prm_[kp_inak] = 0.092455;
    this->prm_[Ka_inak] = 0.0011001;
    this->prm_[Kp_inak] = 5.7392;
    this->prm_[Ka_ina] = 0.10988;
    this->prm_[Kp_ina] = 7.8605;
    this->prm_[ka_ina] = 0.01368;
    this->prm_[kp_ina] = 0.052811;
    this->prm_[ka_plb] = 0.11348;
    this->prm_[kp_plb] = 0.48302;
    this->prm_[Ka_plb] =  9.88539999999999992e-04;
    this->prm_[Kp_plb] = 0.80737;
    this->prm_[ka_ikur] = 0.069537;
    this->prm_[kp_ikur] = 0.317;
    this->prm_[Ka_ikur] = 0.27623;
    this->prm_[Kp_ikur] = 0.002331;
    this->prm_[ka_iks] = 0.16305;
    this->prm_[kp_iks] = 1.0542;
    this->prm_[Ka_iks] =  9.97940000000000003e-05;
    this->prm_[Kp_iks] =  1.11470000000000002e-04;
    this->prm_[M] = 0.01;
    this->prm_[iks_sig_L] = 0.0001;
    this->prm_[iks_sig_K] = 0.01;
    this->prm_[Yotiao] = 0.025;
    this->prm_[IKs_tot] = 0.025;
    this->prm_[iks_sig_PKAf_sum] = 1.0 + (this->prm_[Yotiao] - this->prm_[PKA_eca]) / this->prm_[M];
    this->prm_[iks_sig_PKAf] = this->prm_[M] / 2.0 * (std::sqrt(std::pow(this->prm_[iks_sig_PKAf_sum], 2.0) + 4.0 * this->prm_[PKA_eca] / this->prm_[M]) - this->prm_[iks_sig_PKAf_sum]);
    this->prm_[iks_sig_PP1f_eca_sum] = 1.0 + (this->prm_[Yotiao] - this->prm_[PP1_eca]) / this->prm_[iks_sig_K];
    this->prm_[PP1f_eca] = this->prm_[iks_sig_K] / 2.0 * (std::sqrt(std::pow(this->prm_[iks_sig_PP1f_eca_sum], 2.0) + 4.0 * this->prm_[PP1_eca] / this->prm_[iks_sig_K]) - this->prm_[iks_sig_PP1f_eca_sum]);
    this->prm_[iks_sig_IKsf_sum] = 1.0 + (this->prm_[Yotiao] - this->prm_[IKs_tot]) / this->prm_[iks_sig_L];
    this->prm_[IKsf] = this->prm_[iks_sig_L] / 2.0 * (std::sqrt(std::pow(this->prm_[iks_sig_IKsf_sum], 2.0) + 4.0 * this->prm_[IKs_tot] / this->prm_[iks_sig_L]) - this->prm_[iks_sig_IKsf_sum]);
    this->prm_[Yotiaof] = (this->prm_[Yotiao] - this->prm_[IKs_tot] + this->prm_[IKsf]) / ((1.0 + this->prm_[PP1f_eca] / this->prm_[iks_sig_K]) * (1.0 + this->prm_[iks_sig_PKAf] / this->prm_[M]));
    this->prm_[IKs_arn] = this->prm_[IKsf] * this->prm_[Yotiaof] * this->prm_[iks_sig_PKAf] / (this->prm_[iks_sig_L] * this->prm_[M]);
    this->prm_[IKs_arp] = this->prm_[IKs_arn] * this->prm_[PP1f_eca] / this->prm_[iks_sig_K];
    this->prm_[ka_tni] = 0.10408;
    this->prm_[kp_tni] = 0.052633;
    this->prm_[Ka_tni] =  2.71430000000000008e-05;
    this->prm_[Kp_tni] = 0.26714;
    this->prm_[RyR_tot] = 0.125;
    this->prm_[RyR_akap] = 0.125;
    this->prm_[ka_ryr] = 0.0025548;
    this->prm_[kp_ryr] = 0.0038257;
    this->prm_[Ka_ryr] =  6.62979999999999944e-05;
    this->prm_[Kp_ryr] = 0.043003;
    this->prm_[Mr] = 0.01;
    this->prm_[Lr] = 0.0001;
    this->prm_[Kr] = 0.01;
    this->prm_[akap_sig_RyRf_sum] = 1.0 + (this->prm_[RyR_akap] - this->prm_[RyR_tot]) / this->prm_[Lr];
    this->prm_[RyRf] = this->prm_[Lr] / 2.0 * (std::sqrt(std::pow(this->prm_[akap_sig_RyRf_sum], 2.0) + 4.0 * this->prm_[RyR_tot] / this->prm_[Lr]) - this->prm_[akap_sig_RyRf_sum]);
    this->prm_[ICaL_tot] = 0.025;
    this->prm_[ICaL_akap] = 0.025;
    this->prm_[ka_ical] =  5.10090000000000044e-04;
    this->prm_[kp_ical] = 0.0006903;
    this->prm_[Ka_ical] =  1.27019999999999993e-06;
    this->prm_[Kp_ical] = 0.0063064;
    this->prm_[Mi] = 0.01;
    this->prm_[Li] = 0.0001;
    this->prm_[Ki] = 0.01;
    this->prm_[akap_sig_ICaLf_sum] = 1.0 + (this->prm_[ICaL_akap] - this->prm_[ICaL_tot]) / this->prm_[Li];
    this->prm_[ICaLf] = this->prm_[Li] / 2.0 * (std::sqrt(std::pow(this->prm_[akap_sig_ICaLf_sum], 2.0) + 4.0 * this->prm_[ICaL_tot] / this->prm_[Li]) - this->prm_[akap_sig_ICaLf_sum]);
    this->prm_[akap_sig_PP1f_cav_b] = this->prm_[ICaL_akap] + this->prm_[RyR_akap] + this->prm_[Ki] + this->prm_[Kr] - this->prm_[PP1_cav];
    this->prm_[akap_sig_PP1f_cav_c] = this->prm_[ICaL_akap] * this->prm_[Kr] + this->prm_[RyR_akap] * this->prm_[Ki] + this->prm_[Ki] * this->prm_[Kr] - this->prm_[PP1_cav] * (this->prm_[Ki] + this->prm_[Kr]);
    this->prm_[akap_sig_PP1f_cav_d] = this->prm_[PP1_cav] * this->prm_[Ki] * this->prm_[Kr];
    this->prm_[akap_sig_PP1f_cav_rr] = -this->prm_[akap_sig_PP1f_cav_d] / 27.0 * std::pow(this->prm_[akap_sig_PP1f_cav_b], 3.0) - this->prm_[akap_sig_PP1f_cav_b] * this->prm_[akap_sig_PP1f_cav_b] * this->prm_[akap_sig_PP1f_cav_c] * this->prm_[akap_sig_PP1f_cav_c] / 108.0 + this->prm_[akap_sig_PP1f_cav_b] * this->prm_[akap_sig_PP1f_cav_c] * this->prm_[akap_sig_PP1f_cav_d] / 6.0 + std::pow(this->prm_[akap_sig_PP1f_cav_c], 3.0) / 27.0 + this->prm_[akap_sig_PP1f_cav_d] * this->prm_[akap_sig_PP1f_cav_d] / 4.0;

    this->prm_[akap_sig_PP1f_cav_yi] = 0.;
    if (this->prm_[akap_sig_PP1f_cav_rr] < 0.0) this->prm_[akap_sig_PP1f_cav_yi] = std::sqrt(-this->prm_[akap_sig_PP1f_cav_rr]);

    this->prm_[akap_sig_PP1f_cav_yr] = this->prm_[akap_sig_PP1f_cav_d] / 2.0 + this->prm_[akap_sig_PP1f_cav_b] * this->prm_[akap_sig_PP1f_cav_c] / 6.0 - std::pow(this->prm_[akap_sig_PP1f_cav_b], 3.0) / 27.0;
    if (this->prm_[akap_sig_PP1f_cav_rr] > 0.0) this->prm_[akap_sig_PP1f_cav_yr] += std::sqrt(this->prm_[akap_sig_PP1f_cav_rr]); 

    this->prm_[akap_sig_PP1f_cav_mag] = (this->prm_[akap_sig_PP1f_cav_yr] * this->prm_[akap_sig_PP1f_cav_yr] + this->prm_[akap_sig_PP1f_cav_yi] * std::pow(this->prm_[akap_sig_PP1f_cav_yi], (1.0 / 6.0)));
    this->prm_[akap_sig_PP1f_cav_arg] = std::atan(this->prm_[akap_sig_PP1f_cav_yi] / this->prm_[akap_sig_PP1f_cav_yr]) / 3.0;
    this->prm_[akap_sig_PP1f_cav_x] = (this->prm_[akap_sig_PP1f_cav_c] / 3.0 - this->prm_[akap_sig_PP1f_cav_b] * this->prm_[akap_sig_PP1f_cav_b] / 9.0) / (this->prm_[akap_sig_PP1f_cav_mag] * this->prm_[akap_sig_PP1f_cav_mag]);
    this->prm_[PP1f_cav] = this->prm_[akap_sig_PP1f_cav_mag] * cos(this->prm_[akap_sig_PP1f_cav_arg]) * (1.0 - this->prm_[akap_sig_PP1f_cav_x]) - this->prm_[akap_sig_PP1f_cav_b] / 3.0;
    this->prm_[akap_sig_PKAf_d] = this->prm_[PKA_cav] * this->prm_[Mi] * this->prm_[Mr];
    this->prm_[akap_sig_PKAf_b] = this->prm_[ICaL_akap] + this->prm_[RyR_akap] + this->prm_[Mi] + this->prm_[Mr] - this->prm_[PKA_cav];
    this->prm_[akap_sig_PKAf_c] = this->prm_[ICaL_akap] * this->prm_[Mr] + this->prm_[RyR_akap] * this->prm_[Mi] + this->prm_[Mi] * this->prm_[Mr] - this->prm_[PKA_cav] * (this->prm_[Mi] + this->prm_[Mr]);
    this->prm_[akap_sig_PKAf_rr] = -this->prm_[akap_sig_PKAf_d] / 27.0 * std::pow(this->prm_[akap_sig_PKAf_b], 3.0) - this->prm_[akap_sig_PKAf_b] * this->prm_[akap_sig_PKAf_b] * this->prm_[akap_sig_PKAf_c] * this->prm_[akap_sig_PKAf_c] / 108.0 + this->prm_[akap_sig_PKAf_b] * this->prm_[akap_sig_PKAf_c] * this->prm_[akap_sig_PKAf_d] / 6.0 + std::pow(this->prm_[akap_sig_PKAf_c], 3.0) / 27.0 + this->prm_[akap_sig_PKAf_d] * this->prm_[akap_sig_PKAf_d] / 4.0;

    this->prm_[akap_sig_PKAf_yr] = this->prm_[akap_sig_PKAf_d] / 2.0 + this->prm_[akap_sig_PKAf_b] * this->prm_[akap_sig_PKAf_c] / 6.0 - std::pow(this->prm_[akap_sig_PKAf_b], 3.0) / 27.0;
    if (this->prm_[akap_sig_PKAf_rr] > 0.0) this->prm_[akap_sig_PKAf_yr] += std::sqrt(this->prm_[akap_sig_PKAf_rr]); 

    this->prm_[akap_sig_PKAf_yi] = 0.;
    if (this->prm_[akap_sig_PKAf_rr] < 0.0) this->prm_[akap_sig_PKAf_yi] = std::sqrt(-this->prm_[akap_sig_PKAf_rr]);

    this->prm_[akap_sig_PKAf_mag] = std::pow((this->prm_[akap_sig_PKAf_yr] * this->prm_[akap_sig_PKAf_yr] + this->prm_[akap_sig_PKAf_yi] * this->prm_[akap_sig_PKAf_yi]), (1.0 / 6.0));
    this->prm_[akap_sig_PKAf_arg] = std::atan(this->prm_[akap_sig_PKAf_yi] / this->prm_[akap_sig_PKAf_yr]) / 3.0;
    this->prm_[akap_sig_PKAf_x] = (this->prm_[akap_sig_PKAf_c] / 3.0 - this->prm_[akap_sig_PKAf_b] * this->prm_[akap_sig_PKAf_b] / 9.0) / (this->prm_[akap_sig_PKAf_mag] * this->prm_[akap_sig_PKAf_mag]);
    this->prm_[akap_sig_PKAf] = this->prm_[akap_sig_PKAf_mag] * std::cos(this->prm_[akap_sig_PKAf_arg]) * (1.0 - this->prm_[akap_sig_PKAf_x]) - this->prm_[akap_sig_PKAf_b] / 3.0;
    this->prm_[RyR_akapf] = (this->prm_[RyR_akap] - this->prm_[RyR_tot] + this->prm_[RyRf]) / ((this->prm_[PP1f_cav] / this->prm_[Kr] + 1.0) * (this->prm_[akap_sig_PKAf] / this->prm_[Mr] + 1.0));
    this->prm_[RyR_arn] = this->prm_[RyRf] * this->prm_[RyR_akapf] * this->prm_[akap_sig_PKAf] / (this->prm_[Lr] * this->prm_[Mr]);
    this->prm_[RyR_arp] = this->prm_[RyR_arn] * this->prm_[PP1f_cav] / this->prm_[Kr];
    this->prm_[ICaL_akapf] = (this->prm_[ICaL_akap] - this->prm_[ICaL_tot] + this->prm_[ICaLf]) / ((this->prm_[PP1f_cav] / this->prm_[Ki] + 1.0) * (this->prm_[akap_sig_PKAf] / this->prm_[Mi] + 1.0));
    this->prm_[ICaL_arn] = this->prm_[ICaLf] * this->prm_[ICaL_akapf] * this->prm_[akap_sig_PKAf] / (this->prm_[Li] * this->prm_[Mi]);
    this->prm_[ICaL_arp] = this->prm_[ICaL_arn] * this->prm_[PP1f_cav] / this->prm_[Ki];
    this->prm_[beta_0] = 0.6667 * 4.75;
    this->prm_[irel_fhat_ratio] = 0.0329 + this->prm_[RyR_arn] / this->prm_[RyR_tot];
    this->prm_[ical_f_hat_ratio] = 0.0269 + this->prm_[ICaL_arn] / this->prm_[ICaL_tot];
    this->prm_[iks_f_hat_ratio] = 0.0306 + this->prm_[IKs_arn] / this->prm_[IKs_tot];

}


void Gong2020::Compute(double v_new, double dt, double stim_current)
{
    using namespace Gng20Var;
    using namespace Gng20Prm;
    using namespace Gng20Cur;

    double inhib1_p_old = this->var_[inhib1_p];

    // Compute b-adrenergic signaling
    if (static_cast<int>(this->prm_[signaling]) == 1) this->ComputeBetaAdrenergicSignaling(dt);
    this->ComputeEffectiveFraction();

    // Concentration of uninhibited PP1 in the cytosolic compartment to calculate the fraction of CaMK that has trapped calmodulin or dCaMKt.
    if (static_cast<int>(this->prm_[electrophys]) == 1) {
        this->prm_[Whole_cell_PP1] = 0.1371;
        if (static_cast<int>(this->prm_[signaling]) == 1) {
            double pp1_PP1f_cyt_sum = this->prm_[pp1_K] - this->prm_[PP1_cyt] + inhib1_p_old;
            double PP1f_cyt = 0.5 * (std::sqrt(pp1_PP1f_cyt_sum*pp1_PP1f_cyt_sum + 4. * this->prm_[pp1_K] * this->prm_[PP1_cyt]) - pp1_PP1f_cyt_sum);
            this->prm_[Whole_cell_PP1] = this->prm_[PP1_cav] / this->prm_[vr_cav] + this->prm_[PP1_eca] / this->prm_[vr_eca] + PP1f_cyt / this->prm_[vr_cyt];
        }

        // Compute electrophysiology model.
        this->ComputeElectrophysiology(v_new, dt, stim_current);
    }

}


void Gong2020::ComputeBetaAdrenergicSignaling(double dt)
{
    using namespace Gng20Var;
    using namespace Gng20Prm;

    // iup_f_plb
    if (this->var_[iup_f_plb] < 0.0) this->var_[iup_f_plb] = 0.0001;
    if (this->var_[iup_f_plb] > 1.0) this->var_[iup_f_plb] = 0.9999;

    // f_tni
    if (this->var_[f_tni] < 0.0) this->var_[f_tni] = 0.0001;
    if (this->var_[f_tni] > 1.0) this->var_[f_tni] = 0.9999;

    // ina_f_ina
    if (this->var_[ina_f_ina] < 0.0) this->var_[ina_f_ina] = 0.0001;
    if (this->var_[ina_f_ina] > 1.0) this->var_[ina_f_ina] = 0.9999;

    // f_inak
    if (this->var_[f_inak] < 0.0) this->var_[f_inak] = 0.0001;
    if (this->var_[f_inak] > 1.0) this->var_[f_inak] = 0.9999;

    // f_ikur
    if (this->var_[f_ikur] < 0.0) this->var_[f_ikur] = 0.0001;
    if (this->var_[f_ikur] > 1.0) this->var_[f_ikur] = 0.9999;

    // Total concentration of non-phosphorylated B1AR in the caveolar subspace
    double beta_cav_Rb1_np_tot = this->prm_[beta_cav_R_b1_tot] - this->var_[beta_cav_Rb1_pka_tot] - this->var_[beta_cav_Rb1_grk_tot];

    // Total concentration of non-phosphorylated B2AR in the caveolar subspace
    double beta_cav_Rb2_np_tot = this->prm_[beta_cav_R_b2_tot] - this->var_[beta_cav_Rb2_pka_tot] - this->var_[beta_cav_Rb2_grk_tot];

    // Concentration of Gi holoenzyme in the caveolar subspace
    double beta_cav_Gi_abg = this->prm_[f_Gi_cav] * this->prm_[Gi_tot] * this->prm_[vr_cav] - this->var_[beta_cav_Gi_aGTP] - this->var_[beta_cav_Gi_aGDP];

    // Concentration of Gs holoenzyme in the caveolar subspace
    double beta_cav_Gs_abg = this->prm_[f_Gs_cav] * this->prm_[Gs_tot] * this->prm_[vr_cav] - this->var_[beta_cav_Gs_aGTP] - this->var_[beta_cav_Gs_aGDP];
    double beta_cav_Gs_f_d = beta_cav_Gs_abg * this->prm_[beta_cav_Gs_f_c33] / this->prm_[beta_cav_Gs_f_a];
    double beta_cav_Gs_f_b = (this->prm_[beta_cav_Gs_f_c11] + this->prm_[beta_cav_Gs_f_c22]) / this->prm_[beta_cav_Gs_f_a] + beta_cav_Rb1_np_tot + beta_cav_Rb2_np_tot - beta_cav_Gs_abg;
    double beta_cav_Gs_f_c = (this->prm_[beta_cav_Gs_f_c22] * (beta_cav_Rb1_np_tot - beta_cav_Gs_abg) + this->prm_[beta_cav_Gs_f_c11] * (beta_cav_Rb2_np_tot - beta_cav_Gs_abg) + this->prm_[beta_cav_Gs_f_c33]) / this->prm_[beta_cav_Gs_f_a];
    double beta_cav_Gs_f_rr = -beta_cav_Gs_f_d / 27.0 * std::pow(beta_cav_Gs_f_b, 3.) - beta_cav_Gs_f_b * beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_c / 108.0 + beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_d / 6.0 + std::pow(beta_cav_Gs_f_c, 3.) / 27.0 + beta_cav_Gs_f_d * beta_cav_Gs_f_d / 4.0;
    
    double beta_cav_Gs_f_yr = 0.5*beta_cav_Gs_f_d + beta_cav_Gs_f_b * beta_cav_Gs_f_c / 6.0 - std::pow(beta_cav_Gs_f_b, 3.0) / 27.0;
    if (beta_cav_Gs_f_rr > 0.0) beta_cav_Gs_f_yr += std::sqrt(beta_cav_Gs_f_rr); 
        
    double beta_cav_Gs_f_yi = 0.;
    if (beta_cav_Gs_f_rr < 0.0) beta_cav_Gs_f_yi += std::sqrt(-beta_cav_Gs_f_rr);
    
    double beta_cav_Gs_f_mag = std::pow((beta_cav_Gs_f_yr * beta_cav_Gs_f_yr + beta_cav_Gs_f_yi * beta_cav_Gs_f_yi), (1.0 / 6.0));
    double beta_cav_Gs_f_arg = std::atan(beta_cav_Gs_f_yi / beta_cav_Gs_f_yr) / 3.0;
    double beta_cav_Gs_f_x = (beta_cav_Gs_f_c / 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b / 9.0) / (beta_cav_Gs_f_mag * beta_cav_Gs_f_mag);
    double beta_cav_Gs_f_r = beta_cav_Gs_f_mag * std::cos(beta_cav_Gs_f_arg) * (1.0 - beta_cav_Gs_f_x) - beta_cav_Gs_f_b / 3.0;
    double beta_cav_Gs_f_i = beta_cav_Gs_f_mag * std::sin(beta_cav_Gs_f_arg) * (1.0 + beta_cav_Gs_f_x);

    // Concentration of free Gs in the caveolar subspace
    double beta_cav_Gs_f = std::sqrt(beta_cav_Gs_f_r * beta_cav_Gs_f_r + beta_cav_Gs_f_i * beta_cav_Gs_f_i);

    // Concentration of free non-phosphorylated beta1AR in caveolar subspace
    double beta_cav_Rb1_f = beta_cav_Rb1_np_tot / (1.0 + this->prm_[iso_L] / this->prm_[k_b1_l] + beta_cav_Gs_f * (this->prm_[k_b1_h] + this->prm_[iso_L]) / (this->prm_[k_b1_c] * this->prm_[k_b1_h]));

    // Concentration of non-phosphorylated Ligand / Receptor1 complexes in caveolar subspace
    double beta_cav_LRb1 = this->prm_[iso_L] * beta_cav_Rb1_f / this->prm_[k_b1_l];

    // Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in caveolar subspace
    double beta_cav_LRb1Gs = this->prm_[iso_L] * beta_cav_Rb1_f * beta_cav_Gs_f / (this->prm_[k_b1_c] * this->prm_[k_b1_h]);

    // Concentration of free non-phosphorylated beta2AR in caveolar subspace
    double beta_cav_Rb2_f = beta_cav_Rb2_np_tot / (1.0 + this->prm_[iso_L] / this->prm_[k_b2_l] + beta_cav_Gs_f * (this->prm_[k_b2_h] + this->prm_[iso_L]) / (this->prm_[k_b2_c] * this->prm_[k_b2_h]));

    // Concentration of non-phosphorylated Ligand / Receptor2 complexes in caveolar subspace
    double beta_cav_LRb2 = this->prm_[iso_L] * beta_cav_Rb2_f / this->prm_[k_b2_l];

    // Concentration of non-phosphorylated Receptor2 / G-protein complexes in caveolar subspace
    double beta_cav_Rb2Gs = beta_cav_Rb2_f * beta_cav_Gs_f / this->prm_[k_b2_c];

    // Concentration of non-phosphorylated Receptor1 / G-protein complexes in caveolar subspace
    double beta_cav_Rb1Gs = beta_cav_Rb1_f * beta_cav_Gs_f / this->prm_[k_b1_c];

    // Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in caveolar subspace
    double beta_cav_LRb2Gs = this->prm_[iso_L] * beta_cav_Rb2_f * beta_cav_Gs_f / (this->prm_[k_b2_c] * this->prm_[k_b2_h]);

    // Concentration of total PKA-phosphorylated beta1 receptors
    double dbeta_cav_Rb1_pka_tot = 0.001 * (this->prm_[k_pka_p] * this->var_[pka_cav_C] * beta_cav_Rb1_np_tot - this->prm_[k_pka_dp] * this->var_[beta_cav_Rb1_pka_tot]);
    this->var_[beta_cav_Rb1_pka_tot] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Rb1_pka_tot], dt, dbeta_cav_Rb1_pka_tot);

    // Concentration of total GRK-phosphorylated beta1 receptors
    double dbeta_cav_Rb1_grk_tot = 0.001 * (this->prm_[k_grk_p]* this->prm_[beta_cav_GRK] * (beta_cav_LRb1 + beta_cav_LRb1Gs) - this->prm_[k_grk_dp] * this->var_[beta_cav_Rb1_grk_tot]);
    this->var_[beta_cav_Rb1_grk_tot] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Rb1_grk_tot], dt, dbeta_cav_Rb1_grk_tot);

    // Concentration of total PKA-phosphorylated beta2 receptors
    double dbeta_cav_Rb2_pka_tot = 0.001 * (this->prm_[k_pka_p] * this->var_[pka_cav_C] * beta_cav_Rb2_np_tot - this->prm_[k_pka_dp] * this->var_[beta_cav_Rb2_pka_tot]);
    this->var_[beta_cav_Rb2_pka_tot] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Rb2_pka_tot], dt, dbeta_cav_Rb2_pka_tot);
    
    // Concentration of total GRK-phosphorylated beta2 receptors
    double dbeta_cav_Rb2_grk_tot = 0.001 * (this->prm_[k_grk_p] * this->prm_[beta_cav_GRK] * (beta_cav_LRb2 + beta_cav_LRb2Gs) - this->prm_[k_grk_dp] * this->var_[beta_cav_Rb2_grk_tot]);
    this->var_[beta_cav_Rb2_grk_tot] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Rb2_grk_tot], dt, dbeta_cav_Rb2_grk_tot);

    // Concentration of Gs holoenzyme in the extracaveolar space
    double beta_eca_Gs_abg = this->prm_[f_Gs_eca] * this->prm_[Gs_tot] * this->prm_[vr_eca] - this->var_[beta_eca_Gs_aGTP] - this->var_[beta_eca_Gs_aGDP];

    // Total concentration of non-phosphorylated B2AR in the extracaveolar space
    double beta_eca_Rb2_np_tot = this->prm_[beta_eca_R_b2_tot] - this->var_[beta_eca_Rb2_pka_tot] - this->var_[beta_eca_Rb2_grk_tot];

    // Total concentration of non-phosphorylated B1AR in the extracaveolar space
    double beta_eca_Rb1_np_tot = this->prm_[beta_eca_R_b1_tot] - this->var_[beta_eca_Rb1_pka_tot] - this->var_[beta_eca_Rb1_grk_tot];
    double beta_eca_Gs_f_d = beta_eca_Gs_abg * this->prm_[beta_eca_Gs_f_c33] / this->prm_[beta_eca_Gs_f_a];
    double beta_eca_Gs_f_c = (this->prm_[beta_eca_Gs_f_c22] * (beta_eca_Rb1_np_tot - beta_eca_Gs_abg) + this->prm_[beta_eca_Gs_f_c11] * (beta_eca_Rb2_np_tot - beta_eca_Gs_abg) + this->prm_[beta_eca_Gs_f_c33]) / this->prm_[beta_eca_Gs_f_a];
    double beta_eca_Gs_f_b = (this->prm_[beta_eca_Gs_f_c11] + this->prm_[beta_eca_Gs_f_c22]) / this->prm_[beta_eca_Gs_f_a] + beta_eca_Rb1_np_tot + beta_eca_Rb2_np_tot - beta_eca_Gs_abg;
    double beta_eca_Gs_f_rr = -beta_eca_Gs_f_d / 27.0 * std::pow(beta_eca_Gs_f_b, 3.0) - beta_eca_Gs_f_b * beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_c / 108.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_d / 6.0 + std::pow(beta_eca_Gs_f_c, 3.0) / 27.0 + beta_eca_Gs_f_d * beta_eca_Gs_f_d / 4.0;
    
    double beta_eca_Gs_f_yi = 0.;
    if (beta_eca_Gs_f_rr < 0.0) beta_eca_Gs_f_yi += std::sqrt(-beta_eca_Gs_f_rr);
    
    double beta_eca_Gs_f_yr = beta_eca_Gs_f_d / 2.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c / 6.0 - std::pow(beta_eca_Gs_f_b, 3.0) / 27.0;
    if (beta_eca_Gs_f_rr > 0.0) beta_eca_Gs_f_yr += std::sqrt(beta_eca_Gs_f_rr); 
    
    double beta_eca_Gs_f_mag = std::pow((beta_eca_Gs_f_yr * beta_eca_Gs_f_yr + beta_eca_Gs_f_yi * beta_eca_Gs_f_yi), (1.0 / 6.0));
    double beta_eca_Gs_f_arg = std::atan(beta_eca_Gs_f_yi / beta_eca_Gs_f_yr) / 3.0;
    double beta_eca_Gs_f_x = (beta_eca_Gs_f_c / 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b / 9.0) / (beta_eca_Gs_f_mag * beta_eca_Gs_f_mag);
    double beta_eca_Gs_f_i = beta_eca_Gs_f_mag * std::sin(beta_eca_Gs_f_arg) * (1.0 + beta_eca_Gs_f_x);
    double beta_eca_Gs_f_r = beta_eca_Gs_f_mag * std::cos(beta_eca_Gs_f_arg) * (1.0 - beta_eca_Gs_f_x) - beta_eca_Gs_f_b / 3.0;

    // Concentration of free Gs in the caveolar subspace
    double beta_eca_Gs_f = std::sqrt(beta_eca_Gs_f_r * beta_eca_Gs_f_r + beta_eca_Gs_f_i * beta_eca_Gs_f_i);

    // Concentration of free non-phosphorylated beta1AR in the extracaveolar space
    double beta_eca_Rb1_f = beta_eca_Rb1_np_tot / (1.0 + this->prm_[iso_L] / this->prm_[k_b1_l] + beta_eca_Gs_f * (this->prm_[k_b1_h] + this->prm_[iso_L]) / (this->prm_[k_b1_c] * this->prm_[k_b1_h]));

    // Concentration of free non-phosphorylated beta2AR in the extracaveolar space
    double beta_eca_Rb2_f = beta_eca_Rb2_np_tot / (1.0 + this->prm_[iso_L] / this->prm_[k_b2_l] + beta_eca_Gs_f * (this->prm_[k_b2_h] + this->prm_[iso_L]) / (this->prm_[k_b2_c] * this->prm_[k_b2_h]));

    // Concentration of non-phosphorylated Ligand / Receptor1 complexes in the extracaveolar space
    double beta_eca_LRb1 = this->prm_[iso_L] * beta_eca_Rb1_f / this->prm_[k_b1_l];

    // Concentration of non-phosphorylated Ligand / Receptor2 complexes in the extracaveolar space
    double beta_eca_LRb2 = this->prm_[iso_L] * beta_eca_Rb2_f / this->prm_[k_b2_l];

    // Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in the extracaveolar space
    double beta_eca_LRb2Gs = this->prm_[iso_L] * beta_eca_Rb2_f * beta_eca_Gs_f / (this->prm_[k_b2_c] * this->prm_[k_b2_h]);

    // Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in the extracaveolar space
    double beta_eca_LRb1Gs = this->prm_[iso_L] * beta_eca_Rb1_f * beta_eca_Gs_f / (this->prm_[k_b1_c] * this->prm_[k_b1_h]);

    // Concentration of non-phosphorylated Receptor2 / G-protein complexes in the extracaveolar space
    double beta_eca_Rb2Gs = beta_eca_Rb2_f * beta_eca_Gs_f / this->prm_[k_b2_c];

    // Concentration of non-phosphorylated Receptor1 / G-protein complexes in the extracaveolar space
    double beta_eca_Rb1Gs = beta_eca_Rb1_f * beta_eca_Gs_f / this->prm_[k_b1_c];
    double beta_eca_RGs_tot = beta_eca_Rb1Gs + this->prm_[beta_eca_k_GsAct_b2] * beta_eca_Rb2Gs;
    double beta_eca_LRGs_tot = beta_eca_LRb1Gs + this->prm_[beta_eca_k_GsAct_b2] * beta_eca_LRb2Gs;

    // Concentration of Gi holoenzyme in the extracaveolar space
    double beta_eca_Gi_abg = this->prm_[f_Gi_eca] * this->prm_[Gi_tot] * this->prm_[vr_eca] - this->var_[beta_eca_Gi_aGTP] - this->var_[beta_eca_Gi_aGDP];
    double beta_eca_Rb2_pka_f_c = -this->var_[beta_eca_Rb2_pka_tot] * this->prm_[k_b2_a] * this->prm_[k_b2_f];
    double beta_eca_Rb2_pka_f_b = beta_eca_Gi_abg * (this->prm_[iso_L] + this->prm_[k_b2_f]) - this->var_[beta_eca_Rb2_pka_tot] * (this->prm_[k_b2_f] + this->prm_[iso_L]) + this->prm_[k_b2_a] * this->prm_[k_b2_f] * (1.0 + this->prm_[iso_L] / this->prm_[k_b2_n]);

    // Concentration of free PKA-phosphorylated beta2AR in the extracaveolar space
    double beta_eca_Rb2_pka_f = (-beta_eca_Rb2_pka_f_b + std::sqrt(beta_eca_Rb2_pka_f_b * beta_eca_Rb2_pka_f_b - 4.0 * this->prm_[beta_eca_Rb2_pka_f_a] * beta_eca_Rb2_pka_f_c)) / (2.0 * this->prm_[beta_eca_Rb2_pka_f_a]);

    // Concentration of total PKA-phosphorylated beta1 receptors in the extracaveolar space
    double dbeta_eca_Rb1_pka_tot = 0.001 * (this->prm_[k_pka_p] * this->var_[pka_eca_C] * beta_eca_Rb1_np_tot - this->prm_[k_pka_dp] * this->var_[beta_eca_Rb1_pka_tot]);
    this->var_[beta_eca_Rb1_pka_tot] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Rb1_pka_tot], dt, dbeta_eca_Rb1_pka_tot);

    // Concentration of total GRK-phosphorylated beta1 receptors in the extracaveolar space
    double dbeta_eca_Rb1_grk_tot = 0.001 * (this->prm_[k_grk_p] * this->prm_[beta_eca_GRK] * (beta_eca_LRb1 + beta_eca_LRb1Gs) - this->prm_[k_grk_dp] * this->var_[beta_eca_Rb1_grk_tot]);
    this->var_[beta_eca_Rb1_grk_tot] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Rb1_grk_tot], dt, dbeta_eca_Rb1_grk_tot);

    // Concentration of total PKA-phosphorylated beta2 receptors in the extracaveolar space
    double dbeta_eca_Rb2_pka_tot = 0.001 * (this->prm_[k_pka_p] * this->var_[pka_eca_C] * beta_eca_Rb2_np_tot - this->prm_[k_pka_dp] * this->var_[beta_eca_Rb2_pka_tot]);
    this->var_[beta_eca_Rb2_pka_tot] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Rb2_pka_tot], dt, dbeta_eca_Rb2_pka_tot);

    // Concentration of total GRK-phosphorylated beta2 receptors in the extracaveolar space
    double dbeta_eca_Rb2_grk_tot = 0.001 * (this->prm_[k_grk_p] * this->prm_[beta_eca_GRK] * (beta_eca_LRb2 + beta_eca_LRb2Gs) - this->prm_[k_grk_dp] * this->var_[beta_eca_Rb2_grk_tot]);
    this->var_[beta_eca_Rb2_grk_tot] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Rb2_grk_tot], dt, dbeta_eca_Rb2_grk_tot);

    // Concentration of Gs holoenzyme in the cytoplasm
    double beta_cyt_Gs_abg = this->prm_[f_Gs_cyt] * this->prm_[Gs_tot] * this->prm_[vr_cyt] - this->var_[beta_cyt_Gs_aGTP] - this->var_[beta_cyt_Gs_aGDP];

    // Total concentration of non-phosphorylated beta-1 AR in the cytoplasm
    double beta_cyt_Rb1_np_tot = this->prm_[beta_cyt_R_b1_tot] - this->var_[beta_cyt_Rb1_pka_tot] - this->var_[beta_cyt_Rb1_grk_tot];
    double beta_cyt_Rb1_np_f_b = beta_cyt_Gs_abg * (this->prm_[k_b1_h] + this->prm_[iso_L]) - beta_cyt_Rb1_np_tot * (this->prm_[k_b1_h] + this->prm_[iso_L]) + this->prm_[k_b1_c] * this->prm_[k_b1_h] * (1.0 + this->prm_[iso_L] / this->prm_[k_b1_l]);
    double beta_cyt_Rb1_np_f_c = -beta_cyt_Rb1_np_tot * this->prm_[k_b1_h] * this->prm_[k_b1_c];

    // Concentration of free non-phosphorylated beta-1 AR in the cytoplasm
    double Rb1_np_f = (-beta_cyt_Rb1_np_f_b + std::sqrt(beta_cyt_Rb1_np_f_b * beta_cyt_Rb1_np_f_b - 4.0 * this->prm_[beta_cyt_Rb1_np_f_a] * beta_cyt_Rb1_np_f_c)) / (2.0 * this->prm_[beta_cyt_Rb1_np_f_a]);

    // Concentration of free (non-complexed) Gi in the cytoplasm
    double beta_cyt_Gs_f = beta_cyt_Gs_abg / (1.0 + Rb1_np_f / this->prm_[k_b1_c] * (1.0 + this->prm_[iso_L] / this->prm_[k_b1_h]));

    // Concentration of non-phosphorylated ligand / beta-1 AR complexes in the cytoplasm
    double LRb1_np = this->prm_[iso_L] * Rb1_np_f / this->prm_[k_b1_l];

    // Concentration of non-phosphorylated ligand / beta-1 AR / Gs complexes in the cytoplasm
    double LRb1Gs_np = this->prm_[iso_L] * Rb1_np_f * beta_cyt_Gs_f / (this->prm_[k_b1_c] * this->prm_[k_b1_h]);

    // Concentration of non-phosphorylated beta-1 AR / Gs complexes in the cytoplasm
    double Rb1Gs_np = beta_cyt_Gs_f * Rb1_np_f / this->prm_[k_b1_c];

    // Concentration of total PKA-phosphorylated receptors in the cytoplasm
    double dbeta_cyt_Rb1_pka_tot = 0.001 * (this->prm_[k_pka_p] * this->var_[pka_cyt_C] * beta_cyt_Rb1_np_tot - this->prm_[k_pka_dp] * this->var_[beta_cyt_Rb1_pka_tot]);
    this->var_[beta_cyt_Rb1_pka_tot] = ALGORITHM::ForwardEuler(this->var_[beta_cyt_Rb1_pka_tot], dt, dbeta_cyt_Rb1_pka_tot);

    // Concentration of total GRK-phosphorylated receptors in the cytoplasm
    double dbeta_cyt_Rb1_grk_tot = 0.001 * (this->prm_[k_grk_p] * this->prm_[beta_cyt_GRK] * (LRb1_np + LRb1Gs_np) - this->prm_[k_grk_dp] * this->var_[beta_cyt_Rb1_grk_tot]);
    this->var_[beta_cyt_Rb1_grk_tot] = ALGORITHM::ForwardEuler(this->var_[beta_cyt_Rb1_grk_tot], dt, dbeta_cyt_Rb1_grk_tot);

    //// function Mod_GprotAct() in the other version //////
    double beta_cav_RGs_tot = beta_cav_Rb1Gs + this->prm_[beta_cav_k_GsAct_b2] * beta_cav_Rb2Gs;
    double beta_cav_LRGs_tot = beta_cav_LRb1Gs + this->prm_[beta_cav_k_GsAct_b2] * beta_cav_LRb2Gs;
    double beta_cav_Rb2_pka_f_c = -this->var_[beta_cav_Rb2_pka_tot] * this->prm_[k_b2_a] * this->prm_[k_b2_f];
    double beta_cav_Rb2_pka_f_b = beta_cav_Gi_abg * (this->prm_[iso_L] + this->prm_[k_b2_f]) - this->var_[beta_cav_Rb2_pka_tot] * (this->prm_[k_b2_f] + this->prm_[iso_L]) + this->prm_[k_b2_a] * this->prm_[k_b2_f] * (1.0 + this->prm_[iso_L] / this->prm_[k_b2_n]);

    // Concentration of free PKA-phosphorylated beta2AR in the caveolar subspace
    double beta_cav_Rb2_pka_f = (-beta_cav_Rb2_pka_f_b + std::sqrt(beta_cav_Rb2_pka_f_b * beta_cav_Rb2_pka_f_b - 4.0 * this->prm_[beta_cav_Rb2_pka_f_a] * beta_cav_Rb2_pka_f_c)) / (2.0 * this->prm_[beta_cav_Rb2_pka_f_a]);

    // Concentration of free (non-complexed) Gi in the caveolar subspace
    double beta_cav_Gi_f = beta_cav_Gi_abg / (1.0 + beta_cav_Rb2_pka_f / this->prm_[k_b2_a] * (1.0 + this->prm_[iso_L] / this->prm_[k_b2_f]));

    // Concentration of PKA phosphorylated b2AR/Gi complexes in caveolar subspace
    double beta_cav_Rb2Gi = beta_cav_Rb2_pka_f * beta_cav_Gi_f / this->prm_[k_b2_a];

    // Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in caveolar subspace
    double beta_cav_LRb2Gi = beta_cav_Rb2Gi * this->prm_[iso_L] / this->prm_[k_b2_f];

    // Concentration of free (non-complexed) Gi in the extracaveolar space
    double beta_eca_Gi_f = beta_eca_Gi_abg / (1.0 + beta_eca_Rb2_pka_f / this->prm_[k_b2_a] * (1.0 + this->prm_[iso_L] / this->prm_[k_b2_f]));

    // Concentration of PKA phosphorylated b2AR/Gi complexes in extracaveolar space
    double beta_eca_Rb2Gi = beta_eca_Rb2_pka_f * beta_eca_Gi_f / this->prm_[k_b2_a];

    // Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in extracaveolar space
    double beta_eca_LRb2Gi = this->prm_[iso_L] / this->prm_[k_b2_f] * beta_eca_Rb2Gi;

    // Concentration of active Gs alpha subunit in caveolar subspace (tri-phosphate)
    double dbeta_cav_Gs_aGTP = 0.001 * (this->prm_[k_act2_Gs] * beta_cav_RGs_tot + this->prm_[k_act1_Gs] * beta_cav_LRGs_tot - this->prm_[k_hydr_Gs] * this->var_[beta_cav_Gs_aGTP]);
    this->var_[beta_cav_Gs_aGTP] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Gs_aGTP], dt, dbeta_cav_Gs_aGTP);
    
    // Concentration of active Gi alpha subunit in caveolar subspace
    double dbeta_cav_Gi_aGTP = 0.001 * (this->prm_[k_act2_Gi] * beta_cav_Rb2Gi + this->prm_[k_act1_Gi] * beta_cav_LRb2Gi - this->prm_[k_hydr_Gi] * this->var_[beta_cav_Gi_aGTP]);
    this->var_[beta_cav_Gi_aGTP] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Gi_aGTP], dt, dbeta_cav_Gi_aGTP);

    // Concentration of active Gs alpha subunit in the extracaveolar space (tri-phosphate)
    double dbeta_eca_Gs_aGTP = 0.001 * (this->prm_[k_act2_Gs] * beta_eca_RGs_tot + this->prm_[k_act1_Gs] * beta_eca_LRGs_tot - this->prm_[k_hydr_Gs] * this->var_[beta_eca_Gs_aGTP]);
    this->var_[beta_eca_Gs_aGTP] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Gs_aGTP], dt, dbeta_eca_Gs_aGTP);

    // Concentration of active Gi alpha subunit in the extracaveolar space
    double dbeta_eca_Gi_aGTP = 0.001 * (this->prm_[k_act2_Gi] * beta_eca_Rb2Gi + this->prm_[k_act1_Gi] * beta_eca_LRb2Gi - this->prm_[k_hydr_Gi] * this->var_[beta_eca_Gi_aGTP]);
    this->var_[beta_eca_Gi_aGTP] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Gi_aGTP], dt, dbeta_eca_Gi_aGTP);

    // Concentration of active Gs alpha subunit in cytoplasm
    double dbeta_cyt_Gs_aGTP = 0.001 * (this->prm_[k_act2_Gs] * Rb1Gs_np + this->prm_[k_act1_Gs] * LRb1Gs_np - this->prm_[k_hydr_Gs] * this->var_[beta_cyt_Gs_aGTP]);
    this->var_[beta_cyt_Gs_aGTP] = ALGORITHM::ForwardEuler(this->var_[beta_cyt_Gs_aGTP], dt, dbeta_cyt_Gs_aGTP);

    // Concentration of active Gs beta-gamma subunit in caveolar subspace
    double dbeta_cav_Gs_bg = 0.001 * (this->prm_[k_act2_Gs] * beta_cav_RGs_tot + this->prm_[k_act1_Gs] * beta_cav_LRGs_tot - this->prm_[k_reas_Gs] * this->var_[beta_cav_Gs_bg] * this->var_[beta_cav_Gs_aGDP]);
    this->var_[beta_cav_Gs_bg] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Gs_bg], dt, dbeta_cav_Gs_bg);

    // Concentration of active Gi beta-gamma subunit in caveolar subspace
    double dbeta_cav_Gi_bg = 0.001 * (this->prm_[k_act2_Gi] * beta_cav_Rb2Gi + this->prm_[k_act1_Gi] * beta_cav_LRb2Gi - this->prm_[k_reas_Gi] * this->var_[beta_cav_Gi_bg] * this->var_[beta_cav_Gi_aGDP]);
    this->var_[beta_cav_Gi_bg] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Gi_bg], dt, dbeta_cav_Gi_bg);

    // Concentration of active Gs beta-gamma subunit in the extracaveolar space
    double dbeta_eca_Gs_bg = 0.001 * (this->prm_[k_act2_Gs] * beta_eca_RGs_tot + this->prm_[k_act1_Gs] * beta_eca_LRGs_tot - this->prm_[k_reas_Gs] * this->var_[beta_eca_Gs_bg] * this->var_[beta_eca_Gs_aGDP]);
    this->var_[beta_eca_Gs_bg] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Gs_bg], dt, dbeta_eca_Gs_bg);

    // Concentration of active Gi beta-gamma subunit in the extracaveolar space
    double dbeta_eca_Gi_bg = 0.001 * (this->prm_[k_act2_Gi] * beta_eca_Rb2Gi + this->prm_[k_act1_Gi] * beta_eca_LRb2Gi - this->prm_[k_reas_Gi] * this->var_[beta_eca_Gi_bg] * this->var_[beta_eca_Gi_aGDP]);
    this->var_[beta_eca_Gi_bg] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Gi_bg], dt, dbeta_eca_Gi_bg);

    // Concentration of active Gs beta-gamma subunit in cytoplasm
    double dbeta_cyt_Gs_bg = 0.001 * (this->prm_[k_act2_Gs] * Rb1Gs_np + this->prm_[k_act1_Gs] * LRb1Gs_np - this->prm_[k_reas_Gs] * this->var_[beta_cyt_Gs_bg] * this->var_[beta_cyt_Gs_aGDP]);
    this->var_[beta_cyt_Gs_bg] = ALGORITHM::ForwardEuler(this->var_[beta_cyt_Gs_bg], dt, dbeta_cyt_Gs_bg);

    // Concentration of inactive Gs alpha subunit in caveolar subspace (di-phosphate)
    double dbeta_cav_Gs_aGDP = 0.001 * (this->prm_[k_hydr_Gs] * this->var_[beta_cav_Gs_aGTP] - this->prm_[k_reas_Gs] * this->var_[beta_cav_Gs_bg] * this->var_[beta_cav_Gs_aGDP]);
    this->var_[beta_cav_Gs_aGDP] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Gs_aGDP], dt, dbeta_cav_Gs_aGDP);

    // Concentration of inactive Gi alpha subunit in caveolar subspace
    double dbeta_cav_Gi_aGDP = 0.001 * (this->prm_[k_hydr_Gi] * this->var_[beta_cav_Gi_aGTP] - this->prm_[k_reas_Gi] * this->var_[beta_cav_Gi_bg] * this->var_[beta_cav_Gi_aGDP]);
    this->var_[beta_cav_Gi_aGDP] = ALGORITHM::ForwardEuler(this->var_[beta_cav_Gi_aGDP], dt, dbeta_cav_Gi_aGDP);

    // Concentration of inactive Gs alpha subunit in the extracaveolar space (di-phosphate)
    double dbeta_eca_Gs_aGDP = 0.001 * (this->prm_[k_hydr_Gs] * this->var_[beta_eca_Gs_aGTP] - this->prm_[k_reas_Gs] * this->var_[beta_eca_Gs_bg] * this->var_[beta_eca_Gs_aGDP]);
    this->var_[beta_eca_Gs_aGDP] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Gs_aGDP], dt, dbeta_eca_Gs_aGDP);

    // Concentration of inactive Gi alpha subunit in the extracaveolar space
    double dbeta_eca_Gi_aGDP = 0.001 * (this->prm_[k_hydr_Gi] * this->var_[beta_eca_Gi_aGTP] - this->prm_[k_reas_Gi] * this->var_[beta_eca_Gi_bg] * this->var_[beta_eca_Gi_aGDP]);
    this->var_[beta_eca_Gi_aGDP] = ALGORITHM::ForwardEuler(this->var_[beta_eca_Gi_aGDP], dt, dbeta_eca_Gi_aGDP);
    
    // Concentration of inactive Gs alpha subunit in cytoplasm
    double dbeta_cyt_Gs_aGDP = 0.001 * (this->prm_[k_hydr_Gs] * this->var_[beta_cyt_Gs_aGTP] - this->prm_[k_reas_Gs] * this->var_[beta_cyt_Gs_bg] * this->var_[beta_cyt_Gs_aGDP]);
    this->var_[beta_cyt_Gs_aGDP] = ALGORITHM::ForwardEuler(this->var_[beta_cyt_Gs_aGDP], dt, dbeta_cyt_Gs_aGDP);

    // Concentration of free PKA RC subunits in the caveolar compartment
    double pka_cav_RCf = this->prm_[PKA_cav] - this->var_[pka_cav_ARC] - this->var_[pka_cav_A2RC] - this->var_[pka_cav_A2R];

    // Caveolar concentration of PKA RC dimer with 1 cAMP molecule bound
    double dpka_cav_ARC = 0.001 * (this->prm_[pka_cav_f1] * pka_cav_RCf * this->var_[cAMP_cav] - this->prm_[pka_cav_b1] * this->var_[pka_cav_ARC] - this->prm_[pka_cav_f2] * this->var_[pka_cav_ARC] * this->var_[cAMP_cav] + this->prm_[pka_cav_b2] * this->var_[pka_cav_A2RC]);
    this->var_[pka_cav_ARC] = ALGORITHM::ForwardEuler(this->var_[pka_cav_ARC], dt, dpka_cav_ARC);

    // Caveolar concentration of PKA RC dimer with 2 cAMP molecules bound
    double dpka_cav_A2RC = 0.001 * (this->prm_[pka_cav_f2] * this->var_[pka_cav_ARC] * this->var_[cAMP_cav] - (this->prm_[pka_cav_b2] + this->prm_[pka_cav_f3]) * this->var_[pka_cav_A2RC] + this->prm_[pka_cav_b3] * this->var_[pka_cav_A2R] * this->var_[pka_cav_C]);
    this->var_[pka_cav_A2RC] = ALGORITHM::ForwardEuler(this->var_[pka_cav_A2RC], dt, dpka_cav_A2RC);

    // Caveolar concentration of PKA this->prm_[R] subunit with 2 cAMP molecules bound
    double dpka_cav_A2R = 0.001 * (this->prm_[pka_cav_f3] * this->var_[pka_cav_A2RC] - this->prm_[pka_cav_b3] * this->var_[pka_cav_A2R] * this->var_[pka_cav_C]);
    this->var_[pka_cav_A2R] = ALGORITHM::ForwardEuler(this->var_[pka_cav_A2R], dt, dpka_cav_A2R);

    // Caveolar concentration of free PKA catalytic subunit
    double dpka_cav_C = 0.001 * (this->prm_[pka_cav_f3] * this->var_[pka_cav_A2RC] - this->prm_[pka_cav_b3] * this->var_[pka_cav_A2R] * this->var_[pka_cav_C] + this->prm_[b_pki] * this->var_[pka_cav_PKIC] - this->prm_[f_pki] * (this->prm_[PKI_cav] - this->var_[pka_cav_PKIC]) * this->var_[pka_cav_C]);
    this->var_[pka_cav_C] = ALGORITHM::ForwardEuler(this->var_[pka_cav_C], dt, dpka_cav_C);

    // Caveolar concentration of free PKI inactivated PKA C subunit
    double dpka_cav_PKIC = 0.001 * (this->prm_[f_pki] * (this->prm_[PKI_cav] - this->var_[pka_cav_PKIC]) * this->var_[pka_cav_C] - this->prm_[b_pki] * this->var_[pka_cav_PKIC]);
    this->var_[pka_cav_PKIC] = ALGORITHM::ForwardEuler(this->var_[pka_cav_PKIC], dt, dpka_cav_PKIC);

    // Concentration of free PKA RC subunits in the Extracaveolar compartment
    double pka_eca_RCf = this->prm_[PKA_eca] - this->var_[pka_eca_ARC] - this->var_[pka_eca_A2RC] - this->var_[pka_eca_A2R];

    // Extracaveolar rate of change in free cAMP through binding by PKA
    double pka_eca_dcAMP = -this->prm_[pka_eca_f1] * pka_eca_RCf * this->var_[cAMP_eca] + this->prm_[pka_eca_b1] * this->var_[pka_eca_ARC] - this->prm_[pka_eca_f2] * this->var_[pka_eca_ARC] * this->var_[cAMP_eca] + this->prm_[pka_eca_b2] * this->var_[pka_eca_A2RC];

    // Extracaveolar concentration of PKA RC dimer with 1 cAMP molecule bound
    double dpka_eca_ARC = 0.001 * (this->prm_[pka_eca_f1] * pka_eca_RCf * this->var_[cAMP_eca] - this->prm_[pka_eca_b1] * this->var_[pka_eca_ARC] - this->prm_[pka_eca_f2] * this->var_[pka_eca_ARC] * this->var_[cAMP_eca] + this->prm_[pka_eca_b2] * this->var_[pka_eca_A2RC]);
    this->var_[pka_eca_ARC] = ALGORITHM::ForwardEuler(this->var_[pka_eca_ARC], dt, dpka_eca_ARC);

    // Extracaveolar concentration of PKA RC dimer with 2 cAMP molecules bound
    double dpka_eca_A2RC = 0.001 * (this->prm_[pka_eca_f2] * this->var_[pka_eca_ARC] * this->var_[cAMP_eca] - (this->prm_[pka_eca_b2] + this->prm_[pka_eca_f3]) * this->var_[pka_eca_A2RC] + this->prm_[pka_eca_b3] * this->var_[pka_eca_A2R] * this->var_[pka_eca_C]);
    this->var_[pka_eca_A2RC] = ALGORITHM::ForwardEuler(this->var_[pka_eca_A2RC], dt, dpka_eca_A2RC);

    // Extracaveolar concentration of PKA this->prm_[R] subunit with 2 cAMP molecules bound
    double dpka_eca_A2R = 0.001 * (this->prm_[pka_eca_f3] * this->var_[pka_eca_A2RC] - this->prm_[pka_eca_b3] * this->var_[pka_eca_A2R] * this->var_[pka_eca_C]);
    this->var_[pka_eca_A2R] = ALGORITHM::ForwardEuler(this->var_[pka_eca_A2R], dt, dpka_eca_A2R);

    // Extracaveolar concentration of free PKA catalytic subunit
    double dpka_eca_C = 0.001 * (this->prm_[pka_eca_f3] * this->var_[pka_eca_A2RC] - this->prm_[pka_eca_b3] * this->var_[pka_eca_A2R] * this->var_[pka_eca_C] + this->prm_[b_pki] * this->var_[pka_eca_PKIC] - this->prm_[f_pki] * (this->prm_[PKI_eca] - this->var_[pka_eca_PKIC]) * this->var_[pka_eca_C]);
    this->var_[pka_eca_C] = ALGORITHM::ForwardEuler(this->var_[pka_eca_C], dt, dpka_eca_C);

    // Extracaveolar concentration of free PKI inactivated PKA C subunit
    double dpka_eca_PKIC = 0.001 * (this->prm_[f_pki] * (this->prm_[PKI_eca] - this->var_[pka_eca_PKIC]) * this->var_[pka_eca_C] - this->prm_[b_pki] * this->var_[pka_eca_PKIC]);
    this->var_[pka_eca_PKIC] = ALGORITHM::ForwardEuler(this->var_[pka_eca_PKIC], dt, dpka_eca_PKIC);

    // Concentration of free PKA RC subunits in the Cytosolic compartment
    double pka_cyt_RCf = this->prm_[PKA_cyt] - this->var_[pka_cyt_ARC] - this->var_[pka_cyt_A2RC] - this->var_[pka_cyt_A2R];

    // Cytosolic concentration of PKA RC dimer with 1 cAMP molecule bound
    double dpka_cyt_ARC = 0.001 * (this->prm_[pka_cyt_f1] * pka_cyt_RCf * this->var_[cAMP_cyt] - this->prm_[pka_cyt_b1] * this->var_[pka_cyt_ARC] - this->prm_[pka_cyt_f2] * this->var_[pka_cyt_ARC] * this->var_[cAMP_cyt] + this->prm_[pka_cyt_b2] * this->var_[pka_cyt_A2RC]);
    this->var_[pka_cyt_ARC] = ALGORITHM::ForwardEuler(this->var_[pka_cyt_ARC], dt, dpka_cyt_ARC);

    // Cytosolic concentration of PKA RC dimer with 2 cAMP molecules bound
    double dpka_cyt_A2RC = 0.001 * (this->prm_[pka_cyt_f2] * this->var_[pka_cyt_ARC] * this->var_[cAMP_cyt] - (this->prm_[pka_cyt_b2] + this->prm_[pka_cyt_f3]) * this->var_[pka_cyt_A2RC] + this->prm_[pka_cyt_b3] * this->var_[pka_cyt_A2R] * this->var_[pka_cyt_C]);
    this->var_[pka_cyt_A2RC] = ALGORITHM::ForwardEuler(this->var_[pka_cyt_A2RC], dt, dpka_cyt_A2RC);

    // Cytosolic concentration of PKA this->prm_[R] subunit with 2 cAMP molecules bound
    double dpka_cyt_A2R = 0.001 * (this->prm_[pka_cyt_f3] * this->var_[pka_cyt_A2RC] - this->prm_[pka_cyt_b3] * this->var_[pka_cyt_A2R] * this->var_[pka_cyt_C]);
    this->var_[pka_cyt_A2R] = ALGORITHM::ForwardEuler(this->var_[pka_cyt_A2R], dt, dpka_cyt_A2R);

    // Cytosolic concentration of free PKA catalytic subunit
    double dpka_cyt_C = 0.001 * (this->prm_[pka_cyt_f3] * this->var_[pka_cyt_A2RC] - this->prm_[pka_cyt_b3] * this->var_[pka_cyt_A2R] * this->var_[pka_cyt_C] + this->prm_[b_pki] * this->var_[pka_cyt_PKIC] - this->prm_[f_pki] * (this->prm_[PKI_cyt] - this->var_[pka_cyt_PKIC]) * this->var_[pka_cyt_C]);
    this->var_[pka_cyt_C] = ALGORITHM::ForwardEuler(this->var_[pka_cyt_C], dt, dpka_cyt_C);

    // Cytosolic concentration of free PKI inactivated PKA C subunit
    double dpka_cyt_PKIC = 0.001 * (this->prm_[f_pki] * (this->prm_[PKI_cyt] - this->var_[pka_cyt_PKIC]) * this->var_[pka_cyt_C] - this->prm_[b_pki] * this->var_[pka_cyt_PKIC]);
    this->var_[pka_cyt_PKIC] = ALGORITHM::ForwardEuler(this->var_[pka_cyt_PKIC], dt, dpka_cyt_PKIC);

    // Caveolar rate of change in free cAMP through binding by PKA
    double pka_cav_dcAMP = -this->prm_[pka_cav_f1] * pka_cav_RCf * this->var_[cAMP_cav] + this->prm_[pka_cav_b1] * this->var_[pka_cav_ARC] - this->prm_[pka_cav_f2] * this->var_[pka_cav_ARC] * this->var_[cAMP_cav] + this->prm_[pka_cav_b2] * this->var_[pka_cav_A2RC];

    // Cytosolic rate of change in free cAMP through binding by PKA
    double pka_cyt_dcAMP = -this->prm_[pka_cyt_f1] * pka_cyt_RCf * this->var_[cAMP_cyt] + this->prm_[pka_cyt_b1] * this->var_[pka_cyt_ARC] - this->prm_[pka_cyt_f2] * this->var_[pka_cyt_ARC] * this->var_[cAMP_cyt] + this->prm_[pka_cyt_b2] * this->var_[pka_cyt_A2RC];

    //PDE
    // Rate of cAMP degradation by PDE2 in cytosolic subspace
    double dcAMP_PDE2_cyt = this->prm_[PDE2_cyt] * this->prm_[kPDE2] / (1.0 + this->prm_[KmPDE2] / this->var_[cAMP_cyt]);

    // Rate of cAMP degradation by PDE2 in extracaveolar subspace
    double dcAMP_PDE2_eca = this->prm_[PDE2_eca] * this->prm_[kPDE2] / (1.0 + this->prm_[KmPDE2] / this->var_[cAMP_eca]);

    // Rate of cAMP degradation by PDE2 in caveolar subspace
    double dcAMP_PDE2_cav = this->prm_[PDE2_cav] * this->prm_[kPDE2] / (1.0 + this->prm_[KmPDE2] / this->var_[cAMP_cav]);

    // Rate of cAMP degradation by PDE3 in caveolar subspace
    double dcAMP_PDE3_cav = (this->prm_[PDE3_cav] + (this->prm_[delta_k_pde34] - 1.0) * this->var_[PDE3_P_cav]) * this->prm_[kPDE3] / (1.0 + this->prm_[KmPDE3] / this->var_[cAMP_cav]);

    // Rate of cAMP degradation by PDE4 in cytosolic subspace
    double dcAMP_PDE4_cyt = (this->prm_[PDE4_cyt] + (this->prm_[delta_k_pde34] - 1.0) * this->var_[PDE4_P_cyt]) * this->prm_[kPDE4] / (1.0 + this->prm_[KmPDE4] / this->var_[cAMP_cyt]);

    // Rate of cAMP degradation by PDE4 in extracaveolar subspace
    double dcAMP_PDE4_eca = (this->prm_[PDE4_eca] + (this->prm_[delta_k_pde34] - 1.0) * this->var_[PDE4_P_eca]) * this->prm_[kPDE4] / (1.0 + this->prm_[KmPDE4] / this->var_[cAMP_eca]);

    // Rate of cAMP degradation by PDE4 in caveolar subspace
    double dcAMP_PDE4_cav = (this->prm_[PDE4_cav] + (this->prm_[delta_k_pde34] - 1.0) * this->var_[PDE4_P_cav]) * this->prm_[kPDE4] / (1.0 + this->prm_[KmPDE4] / this->var_[cAMP_cav]);

    // Rate of cAMP degradation by PDE3 in cytosolic subspace
    double dcAMP_PDE3_cyt = (this->prm_[PDE3_cyt] + (this->prm_[delta_k_pde34] - 1.0) * this->var_[PDE3_P_cyt]) * this->prm_[kPDE3] / (1.0 + this->prm_[KmPDE3] / this->var_[cAMP_cyt]);

    double camp_cAMP_cyt_pde = dcAMP_PDE2_cyt + dcAMP_PDE3_cyt + dcAMP_PDE4_cyt;
    double camp_cAMP_cyt_j1 = this->prm_[j_cav_cyt] * (this->var_[cAMP_cav] - this->var_[cAMP_cyt]) / this->prm_[v_cyt];
    double camp_cAMP_cyt_j2 = this->prm_[j_eca_cyt] * (this->var_[cAMP_eca] - this->var_[cAMP_cyt]) / this->prm_[v_cyt];
    double camp_cAMP_eca_pde = dcAMP_PDE2_eca + dcAMP_PDE4_eca;
    double camp_cAMP_eca_j2 = this->prm_[j_eca_cyt] * (this->var_[cAMP_eca] - this->var_[cAMP_cyt]) / this->prm_[v_eca];
    double camp_cAMP_eca_j1 = this->prm_[j_cav_eca] * (this->var_[cAMP_cav] - this->var_[cAMP_eca]) / this->prm_[v_eca];
    double camp_cAMP_cav_pde = dcAMP_PDE2_cav + dcAMP_PDE3_cav + dcAMP_PDE4_cav;
    double camp_cAMP_cav_j2 = this->prm_[j_cav_cyt] * (this->var_[cAMP_cav] - this->var_[cAMP_cyt]) / this->prm_[v_cav];
    double camp_cAMP_cav_j1 = this->prm_[j_cav_eca] * (this->var_[cAMP_cav] - this->var_[cAMP_eca]) / this->prm_[v_cav];

    double ac_kAC47_cyt_gsa = std::pow(this->var_[beta_cyt_Gs_aGTP], this->prm_[hGsAC47]);
    double kAC47_cyt = this->prm_[afAC47] * (this->prm_[basalAC47] + ac_kAC47_cyt_gsa / (this->prm_[KmGsAC47] + ac_kAC47_cyt_gsa));
    double ac_kAC56_cav_gsa = std::pow(this->var_[beta_cav_Gs_aGTP], this->prm_[hGsAC56]);
    double gsi = std::pow(this->var_[beta_cav_Gs_aGTP], this->prm_[hGsGiAC56]);
    double kAC56_cav = this->prm_[afAC56] * (this->prm_[basalAC56] + ac_kAC56_cav_gsa / (this->prm_[KmGsAC56] + ac_kAC56_cav_gsa)) * (1.0 - (1.0 - this->prm_[vGsGiAC56] * gsi / (this->prm_[KmGsGiAC56] + gsi)) * this->var_[beta_cav_Gi_bg] / (this->prm_[KmGiAC56] + this->var_[beta_cav_Gi_bg]));
    double ac_kAC47_eca_gsa = std::pow(this->var_[beta_eca_Gs_aGTP], this->prm_[hGsAC47]);
    double kAC47_eca = this->prm_[afAC47] * (this->prm_[basalAC47] + ac_kAC47_eca_gsa / (this->prm_[KmGsAC47] + ac_kAC47_eca_gsa));
    double ac_kAC56_cyt_gsa = std::pow(this->var_[beta_cyt_Gs_aGTP], this->prm_[hGsAC56]);
    double kAC56_cyt = this->prm_[afAC56] * (this->prm_[basalAC56] + ac_kAC56_cyt_gsa / (this->prm_[KmGsAC56] + ac_kAC56_cyt_gsa));

    // Rate of cAMP production by AC type 4/7 in cytoplasm
    double dcAMP_AC47_cyt = kAC47_cyt * this->prm_[AC47_cyt] * this->prm_[fATP];

    // Rate of cAMP production by AC type 5/6 in cytoplasm
    double dcAMP_AC56_cyt = kAC56_cyt * this->prm_[AC56_cyt] * this->prm_[fATP];

    // Rate of cAMP production by AC type 5/6 in caveolar subspace
    double dcAMP_AC56_cav = kAC56_cav * this->prm_[AC56_cav] * this->prm_[fATP];

    // Rate of cAMP production by AC type 4/7 in extracaveolar subspace
    double dcAMP_AC47_eca = kAC47_eca * this->prm_[AC47_eca] * this->prm_[fATP];

    // Caveolar concentration of cAMP
    double dcAMP_cav = 0.001 * (pka_cav_dcAMP + dcAMP_AC56_cav - camp_cAMP_cav_pde - camp_cAMP_cav_j1 - camp_cAMP_cav_j2);
    this->var_[cAMP_cav] = ALGORITHM::ForwardEuler(this->var_[cAMP_cav], dt, dcAMP_cav);

    // Extracaveolar concentration of cAMP
    double dcAMP_eca = 0.001 * (pka_eca_dcAMP + dcAMP_AC47_eca - camp_cAMP_eca_pde + camp_cAMP_eca_j1 - camp_cAMP_eca_j2);
    this->var_[cAMP_eca] = ALGORITHM::ForwardEuler(this->var_[cAMP_eca], dt, dcAMP_eca);

    // Cytosolic concentration of cAMP
    double dcAMP_cyt = 0.001 * (pka_cyt_dcAMP + dcAMP_AC47_cyt + dcAMP_AC56_cyt - camp_cAMP_cyt_pde + camp_cAMP_cyt_j1 + camp_cAMP_cyt_j2);
    this->var_[cAMP_cyt] = ALGORITHM::ForwardEuler(this->var_[cAMP_cyt], dt, dcAMP_cyt);

    // Concentration of phosphorylated PDE3 in the caveolar subspace
    double dPDE3_P_cav = 0.001 * (this->prm_[kfPDEp] * this->var_[pka_cav_C] * (this->prm_[PDE3_cav] - this->var_[PDE3_P_cav]) - this->prm_[kbPDEp] * this->var_[PDE3_P_cav]);
    this->var_[PDE3_P_cav] = ALGORITHM::ForwardEuler(this->var_[PDE3_P_cav], dt, dPDE3_P_cav);

    // Concentration of phosphorylated PDE3 in the cytosolic subspace
    double dPDE3_P_cyt = 0.001 * (this->prm_[kfPDEp] * this->var_[pka_cyt_C] * (this->prm_[PDE3_cyt] - this->var_[PDE3_P_cyt]) - this->prm_[kbPDEp] * this->var_[PDE3_P_cyt]);
    this->var_[PDE3_P_cyt] = ALGORITHM::ForwardEuler(this->var_[PDE3_P_cyt], dt, dPDE3_P_cyt);

    // Concentration of phosphorylated PDE4 in the caveolar subspace
    double dPDE4_P_cav = 0.001 * (this->prm_[kfPDEp] * this->var_[pka_cav_C] * (this->prm_[PDE4_cav] - this->var_[PDE4_P_cav]) - this->prm_[kbPDEp] * this->var_[PDE4_P_cav]);
    this->var_[PDE4_P_cav] = ALGORITHM::ForwardEuler(this->var_[PDE4_P_cav], dt, dPDE4_P_cav);

    // Concentration of phosphorylated PDE4 in the extracaveolar subspace
    double dPDE4_P_eca = 0.001 * (this->prm_[kfPDEp] * this->var_[pka_eca_C] * (this->prm_[PDE4_eca] - this->var_[PDE4_P_eca]) - this->prm_[kbPDEp] * this->var_[PDE4_P_eca]);
    this->var_[PDE4_P_eca] = ALGORITHM::ForwardEuler(this->var_[PDE4_P_eca], dt, dPDE4_P_eca);

    // Concentration of phosphorylated PDE4 in the cytosolic subspace
    double dPDE4_P_cyt = 0.001 * (this->prm_[kfPDEp] * this->var_[pka_cyt_C] * (this->prm_[PDE4_cyt] - this->var_[PDE4_P_cyt]) - this->prm_[kbPDEp] * this->var_[PDE4_P_cyt]);
    this->var_[PDE4_P_cyt] = ALGORITHM::ForwardEuler(this->var_[PDE4_P_cyt], dt, dPDE4_P_cyt);

    //// Mod_PP1_Inhibition()
    double pp1_PP1f_cyt_sum = this->prm_[pp1_K] - this->prm_[PP1_cyt] + this->var_[inhib1_p];

    // Concentration of uninhibited PP1 in the cytosolic compartment
    double PP1f_cyt = 0.5 * (std::sqrt(std::pow(pp1_PP1f_cyt_sum, 2.0) + 4.0 * this->prm_[pp1_K] * this->prm_[PP1_cyt]) - pp1_PP1f_cyt_sum);
    double di = this->prm_[inhib1_tot] - this->var_[inhib1_p];

    // Concentration of phosphorylated PP1 inhibitor 1 (cytoplasmic)
    double dinhib1_p = 0.001 * (this->prm_[kp] * this->var_[pka_cyt_C] * di / (this->prm_[Kp] + di) - this->prm_[kdp] * this->prm_[PP2A] * this->var_[inhib1_p] / (this->prm_[Kdp] + this->var_[inhib1_p]));
    this->var_[inhib1_p] = ALGORITHM::ForwardEuler(this->var_[inhib1_p], dt, dinhib1_p);

    // Substrates without AKAP
    // Fraction of phosphorylated PLB
    double diup_f_plb = 0.001 * (this->prm_[ka_plb] * this->var_[pka_cyt_C] * (1.0 - this->var_[iup_f_plb]) / (this->prm_[Ka_plb] + 1.0 - this->var_[iup_f_plb]) - this->prm_[kp_plb] * PP1f_cyt * this->var_[iup_f_plb] / (this->prm_[Kp_plb] + this->var_[iup_f_plb]));
    this->var_[iup_f_plb] = ALGORITHM::ForwardEuler(this->var_[iup_f_plb], dt, diup_f_plb);

    // Fraction of phosphorylated Troponin
    double df_tni = 0.001 * (this->prm_[ka_tni] * this->var_[pka_cyt_C] * (1.0 - this->var_[f_tni]) / (this->prm_[Ka_tni] + 1.0 - this->var_[f_tni]) - this->prm_[kp_tni] * this->prm_[PP2A] * this->var_[f_tni] / (this->prm_[Kp_tni] + this->var_[f_tni]));
    this->var_[f_tni] = ALGORITHM::ForwardEuler(this->var_[f_tni], dt, df_tni);

    // Fraction of phosphorylated INa channels
    double dina_f_ina = 0.001 * (this->prm_[ka_ina] * this->var_[pka_cav_C] * (1.0 - this->var_[ina_f_ina]) / (this->prm_[Ka_ina] + 1.0 - this->var_[ina_f_ina]) - this->prm_[kp_ina] * this->prm_[PP1_cav] * this->var_[ina_f_ina] / (this->prm_[Kp_ina] + this->var_[ina_f_ina]));
    this->var_[ina_f_ina] = ALGORITHM::ForwardEuler(this->var_[ina_f_ina], dt, dina_f_ina);

    // Fraction of phosphorylated INaK
    double df_inak = 0.001 * (this->prm_[ka_inak] * this->var_[pka_cav_C] * (1.0 - this->var_[f_inak]) / (this->prm_[Ka_inak] + 1.0 - this->var_[f_inak]) - this->prm_[kp_inak] * this->prm_[PP1_cav] * this->var_[f_inak] / (this->prm_[Kp_inak] + this->var_[f_inak]));
    this->var_[f_inak] = ALGORITHM::ForwardEuler(this->var_[f_inak], dt, df_inak);

    // Fraction of phosphorylated IKur channels
    double df_ikur = 0.001 * (this->prm_[ka_ikur] * this->var_[pka_eca_C] * (1.0 - this->var_[f_ikur]) / (this->prm_[Ka_ikur] + 1.0 - this->var_[f_ikur]) - this->prm_[kp_ikur] * this->prm_[PP1_eca] * this->var_[f_ikur] / (this->prm_[Kp_ikur] + this->var_[f_ikur]));
    this->var_[f_ikur] = ALGORITHM::ForwardEuler(this->var_[f_ikur], dt, df_ikur);

    //Substrates with AKAP
    double iks_sig_IKsp_dif = this->prm_[IKs_arp] - this->var_[IKsp];

    // Concentration of phosphorylated IKs channels
    double dIKsp = 0.001 * (this->prm_[ka_iks] * this->var_[pka_eca_C] * iks_sig_IKsp_dif / (this->prm_[Ka_iks] + iks_sig_IKsp_dif) - this->prm_[kp_iks] * this->prm_[PP1_eca] * this->var_[IKsp] / (this->prm_[Kp_iks] + this->var_[IKsp]));
    this->var_[IKsp] = ALGORITHM::ForwardEuler(this->var_[IKsp], dt, dIKsp);

    double akap_sig_RyRp_dif = this->prm_[RyR_arp] - this->var_[RyRp];
    double dRyRp = 0.001 * (this->prm_[ka_ryr] * this->var_[pka_cav_C] * akap_sig_RyRp_dif / (this->prm_[Ka_ryr] + akap_sig_RyRp_dif) - this->prm_[kp_ryr] * this->prm_[PP1_cav] * this->var_[RyRp] / (this->prm_[Kp_ryr] + this->var_[RyRp]));
    this->var_[RyRp] = ALGORITHM::ForwardEuler(this->var_[RyRp], dt, dRyRp);

    double akap_sig_ICaLp_dif = this->prm_[ICaL_arp] - this->var_[ICaLp];

    // Concentration of phosphorylated L-type Calcium channels
    double dICaLp = 0.001 * (this->prm_[ka_ical] * this->var_[pka_cav_C] * akap_sig_ICaLp_dif / (this->prm_[Ka_ical] + akap_sig_ICaLp_dif) - this->prm_[kp_ical] * this->prm_[PP1_cav] * this->var_[ICaLp] / (this->prm_[Kp_ical] + this->var_[ICaLp]));
    this->var_[ICaLp] = ALGORITHM::ForwardEuler(this->var_[ICaLp], dt, dICaLp);

}


void Gong2020::ComputeEffectiveFraction() 
{
    using namespace Gng20Var;
    using namespace Gng20Prm;

    //// ICaL
    double fp_ICaL = (this->var_[ICaLp] + this->prm_[ICaL_arn]) / this->prm_[ICaL_tot]; 
    if (fp_ICaL < 0.0) fp_ICaL = 0.0001;
    if (fp_ICaL > 1.0) fp_ICaL = 0.9999;
    
    double ical_f_hat = (fp_ICaL - this->prm_[ical_f_hat_ratio]) / (0.9273 - this->prm_[ical_f_hat_ratio]);
    if (ical_f_hat < 0.0) ical_f_hat = 0.0;
    if (ical_f_hat > 1.0) ical_f_hat = 1.0;

    //// IKs
    double fp_iks = (this->var_[IKsp] + this->prm_[IKs_arn]) / this->prm_[IKs_tot];
    if (fp_iks < 0.0) fp_iks = 0.0001;
    if (fp_iks > 1.0) fp_iks = 0.9999;

    double iks_f_hat = (fp_iks - this->prm_[iks_f_hat_ratio]) / (0.785 - this->prm_[iks_f_hat_ratio]);
    if (iks_f_hat < 0.0) iks_f_hat = 0.0;
    if (iks_f_hat > 1.0) iks_f_hat = 1.0;

    //// Iup (PLB)
    double iup_f_pka = (this->var_[iup_f_plb] - 0.6591) / (0.9945 - 0.6591);
    if (iup_f_pka < 0.0) iup_f_pka = 0.0;
    if (iup_f_pka > 1.0) iup_f_pka = 1.0;

    //// Tni
    double calcium_fhat = (this->var_[f_tni] - 0.6735188) / (0.9991797 - 0.6735188);
    if (calcium_fhat < 0.0) calcium_fhat = 0.0;
    if (calcium_fhat > 1.0) calcium_fhat = 1.0;

    //// INa
    double ina_f_pka = (this->var_[ina_f_ina] - 0.2394795) / (0.9501431 - 0.2394795);
    if (ina_f_pka < 0.0) ina_f_pka = 0.0;
    if (ina_f_pka > 1.0) ina_f_pka = 1.0;

    //// INaK
    double inak_fhat = (this->var_[f_inak] - 0.1263453) / (0.9980137 - 0.1263453);
    if (inak_fhat < 0.0) inak_fhat = 0.0;
    if (inak_fhat > 1.0) inak_fhat = 1.0;

    //// RyR
    double fp_RyR = (this->var_[RyRp] + this->prm_[RyR_arn]) / this->prm_[RyR_tot];
    if (fp_RyR < 0.0) fp_RyR = 0.0001;
    if (fp_RyR > 1.0) fp_RyR = 0.9999;

    double irel_fhat = (fp_RyR - this->prm_[irel_fhat_ratio]) / (0.9586 - this->prm_[irel_fhat_ratio]);
    if (irel_fhat < 0.0) irel_fhat = 0.0;
    if (irel_fhat > 1.0) irel_fhat = 1.0;

    //// IKur
    double ikur_fhat = (this->var_[f_ikur] -  5.89379800000000009e-02) / (0.393747 -  5.89379800000000009e-02);
    if (ikur_fhat < 0.0) ikur_fhat = 0.0;
    if (ikur_fhat > 1.0) ikur_fhat = 1.0;

    this->prm_[fICaLP] = ical_f_hat;
    this->prm_[fIKsP] = iks_f_hat;
    this->prm_[fPLBP] = iup_f_pka;
    this->prm_[fTnIP] = calcium_fhat;
    this->prm_[fINaP] = ina_f_pka;
    this->prm_[fINaKP] = inak_fhat;
    this->prm_[fRyRP] = irel_fhat;
    this->prm_[fIKurP] = ikur_fhat;
    
}


void Gong2020::ComputeElectrophysiology(double v_new, double dt, double stim_current)
{
    using namespace Gng20Var;
    using namespace Gng20Prm;
    using namespace Gng20Cur;

    double INa_P = this->prm_[fINaP];
    double ICaL_P = this->prm_[fICaLP];
    double IKs_P = this->prm_[fIKsP];
    double INaK_P = this->prm_[fINaKP];
    double RyR_P = this->prm_[fRyRP];
    double SERCA_P = this->prm_[fPLBP];
    double TnI_P = this->prm_[fTnIP];
    double IKb_P = this->prm_[fIKurP];

    //CaMK constants
    this->prm_[aCaMK] = 0.05;
    this->prm_[bCaMK] = 0.00068;
    this->prm_[CaMKo] = 0.05;
    this->prm_[KmCaM] = 0.0015;

    //update CaMK
    double CaMKb = this->prm_[CaMKo]*(1.0-this->var_[CaMKt])/(1.0+this->prm_[KmCaM]/this->var_[cass]);
    double CaMKa = CaMKb+this->var_[CaMKt];

    double betaCaMKII = this->prm_[bCaMK] * (0.1 + (0.9  * this->prm_[Whole_cell_PP1] / 0.1371));
    
    // fraction of CaMK that has trapped calmodulin
    double dCaMKt = this->prm_[aCaMK]*CaMKb*(CaMKb+this->var_[CaMKt])-betaCaMKII*this->var_[CaMKt];
    if (static_cast<int>(this->prm_[PP1block]) == 1) {
        dCaMKt = this->prm_[aCaMK]*CaMKb*(CaMKb+this->var_[CaMKt])-this->prm_[bCaMK]*this->var_[CaMKt];
    }
    this->var_[CaMKt] = ALGORITHM::ForwardEuler(this->var_[CaMKt], dt, dCaMKt); 
    
    //reversal potentials
    double ENa = (this->prm_[R]*this->prm_[T]/this->prm_[F])*std::log(this->prm_[nao]/this->var_[nai]);
    double EK = (this->prm_[R]*this->prm_[T]/this->prm_[F])*std::log(this->prm_[ko]/this->var_[ki]);
    this->prm_[PKNa] = 0.01833;
    double EKs = (this->prm_[R]*this->prm_[T]/this->prm_[F])*log((this->prm_[ko]+this->prm_[PKNa]*this->prm_[nao])/(this->var_[ki]+this->prm_[PKNa]*this->var_[nai]));

    //convenient shorthand calculations
    double vffrt = v_new*this->prm_[F]*this->prm_[F]/(this->prm_[R]*this->prm_[T]);
    double vfrt = v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T]);

    //// calculate INa
    // gating NP
    double mss = 1.0/(1.0+std::exp((-(v_new+39.57))/9.871));
    double tm = 1.0/(6.765*std::exp((v_new+11.64)/34.77)+8.552*std::exp(-(v_new+77.42)/5.955));
    this->var_[m] = ALGORITHM::RushLarsen(mss, this->var_[m], dt, tm);
    
    double hss = 1.0/(1+std::exp((v_new+82.90)/6.086));
    double thf = 1.0/(1.432e-5*std::exp(-(v_new+1.196)/6.285)+6.149*std::exp((v_new+0.5096)/20.27));
    double Ahf = 0.99;
    this->var_[hf] = ALGORITHM::RushLarsen(hss, this->var_[hf], dt, thf);

    double ths = 1.0/(0.009794*std::exp(-(v_new+17.95)/28.05)+0.3343*std::exp((v_new+5.730)/56.66));
    double Ahs = 1.0-Ahf;
    this->var_[hs] = ALGORITHM::RushLarsen(hss, this->var_[hs], dt, ths);

    double h = Ahf*this->var_[hf]+Ahs*this->var_[hs];
    double jss = hss;
    double tj = 2.038+1.0/(0.02136*std::exp(-(v_new+100.6)/8.281)+0.3052*std::exp((v_new+0.9941)/38.45));
    this->var_[j] = ALGORITHM::RushLarsen(jss, this->var_[j], dt, tj);

    // gating CaMK-P
    double hssp = 1.0/(1+std::exp((v_new+89.1)/6.086));
    double thsp = 3.0*ths;
    this->var_[hsp] = ALGORITHM::RushLarsen(hssp, this->var_[hsp], dt, thsp);
    
     
    double hp = Ahf*this->var_[hf]+Ahs*this->var_[hsp];
    double tjp = 1.46*tj;
    this->var_[jp] = ALGORITHM::RushLarsen(jss, this->var_[jp], dt, tjp);

    // gating PKA-P
    // GNa increased after PKA-P
    double GNaP = this->prm_[GNa]*2.7;
    double VP_h = 5.;

    double hPss =  1.0/(1.+std::exp((v_new+82.90 + VP_h)/6.086));
    double dhPf = (hPss - this->var_[hPf])/thf;
    double dhPs = (hPss - this->var_[hPs])/ths;
    this->var_[hPf] = ALGORITHM::ForwardEuler(this->var_[hPf], dt, dhPf); 
    this->var_[hPs] = ALGORITHM::ForwardEuler(this->var_[hPs], dt, dhPs);

    double hP = Ahf*this->var_[hPf]+Ahs*this->var_[hPs];
    double jPss = hPss; 
    this->var_[jP] = ALGORITHM::RushLarsen(jPss, this->var_[jP], dt, tj);

    double hBPss = 1.0/(1.+std::exp((v_new+89.1 + VP_h)/6.086));
    double thBPf = thf; 
    double thBPs = thsp;
    this->var_[hBPf] = ALGORITHM::RushLarsen(hBPss, this->var_[hBPf], dt, thBPf);
    this->var_[hBPs] = ALGORITHM::RushLarsen(hBPss, this->var_[hBPs], dt, thBPs);

    double hBP = Ahf*this->var_[hBPf]+Ahs*this->var_[hBPs];
    double jBPss = hBPss ; 
    double tjBP = tjp ;  
    this->var_[jBP] = ALGORITHM::RushLarsen(jBPss, this->var_[jBP], dt, tjBP);

    // Putting together the channels behavior and fraction
    double fINap=(1.0/(1.0+this->prm_[KmCaMK]/CaMKa)); 
    double fINa_P = INa_P; 
    double fINa_BP = fINap*fINa_P;
    double fINa_CaMKonly = fINap-fINa_BP;
    double fINa_PKAonly = fINa_P-fINa_BP;

    double INa_NP = this->prm_[GNa]*(v_new-ENa)*std::pow(this->var_[m],3.0)*h*this->var_[j]; 
    double INa_CaMK = this->prm_[GNa]*(v_new-ENa)*std::pow(this->var_[m],3.0)*hp*this->var_[jp];
    double INa_PKA = GNaP*(v_new-ENa)*std::pow(this->var_[m],3.0)*hP*this->var_[jP];
    double INa_BP = GNaP*(v_new-ENa)*std::pow(this->var_[m],3.0)*hBP*this->var_[jBP];

    // 4 population 
    this->cur_[INa] = (1.0-this->block_coeff_[INa]) * ((1.0-fINa_CaMKonly-fINa_PKAonly-fINa_BP)*INa_NP + fINa_CaMKonly*INa_CaMK + fINa_PKAonly*INa_PKA + fINa_BP*INa_BP);

    // calculate INaL
    double mLss = 1.0/(1.0+std::exp((-(v_new+42.85))/5.264));
    double tmL = tm;
    this->var_[mL] = ALGORITHM::RushLarsen(mLss, this->var_[mL], dt, tmL);

    double hLss = 1.0/(1.0+std::exp((v_new+87.61)/7.488));
    double thL = 200.0;
    this->var_[hL] = ALGORITHM::RushLarsen(hLss, this->var_[hL], dt, thL);

    double hLssp = 1.0/(1.0+std::exp((v_new+93.81)/7.488));
    double thLp = 3.0*thL;
    this->var_[hLp] = ALGORITHM::RushLarsen(hLssp, this->var_[hLp], dt, thLp);

    double fINaLp = (1.0/(1.0+this->prm_[KmCaMK]/CaMKa));
    this->cur_[INaL] = (1.0-this->block_coeff_[INaL]) * (this->prm_[GNaL]*(v_new-ENa)*this->var_[mL]*((1.0-fINaLp)*this->var_[hL]+fINaLp*this->var_[hLp]));

    //// calculate Ito
    double ass = 1.0/(1.0+std::exp((-(v_new-14.34))/14.82));
    double ta = 1.0515/(1.0/(1.2089*(1.0+std::exp(-(v_new-18.4099)/29.3814)))+3.5/(1.0+std::exp((v_new+100.0)/29.3814)));
    this->var_[a] = ALGORITHM::RushLarsen(ass, this->var_[a], dt, ta);

    double iss = 1.0/(1.0+std::exp((v_new+43.94)/5.711));
    double delta_epi = 1.;
    if (static_cast<int>(this->prm_[celltype]) == 1)  delta_epi -= (0.95/(1.0+std::exp((v_new+70.0)/5.0)));
    
    double tiF = delta_epi * (4.562+1./(0.3933*std::exp((-(v_new+100.0))/100.0)+0.08004*std::exp((v_new+50.0)/16.59)));
    double AiF = 1.0/(1.0+std::exp((v_new-213.6)/151.2));
    this->var_[iF] = ALGORITHM::RushLarsen(iss, this->var_[iF], dt, tiF);

    double tiS = delta_epi * (23.62+1./(0.001416*std::exp((-(v_new+96.52))/59.05)+1.780e-8*std::exp((v_new+114.1)/8.079)));
    double AiS = 1.0-AiF;    
    this->var_[iS] = ALGORITHM::RushLarsen(iss, this->var_[iS], dt, tiS);

    double i = AiF*this->var_[iF]+AiS*this->var_[iS];
    double assp = 1.0/(1.0+std::exp((-(v_new-24.34))/14.82));
    this->var_[ap] = ALGORITHM::RushLarsen(assp, this->var_[ap], dt, ta);

    double dti_develop = 1.354+1.0e-4/(std::exp((v_new-167.4)/15.89)+std::exp(-(v_new-12.23)/0.2154));
    double dti_recover = 1.0-0.5/(1.0+std::exp((v_new+70.0)/20.0));
    double tiFp = dti_develop*dti_recover*tiF;
    double tiSp = dti_develop*dti_recover*tiS;
    this->var_[iFp] = ALGORITHM::RushLarsen(iss, this->var_[iFp], dt, tiFp);
    this->var_[iSp] = ALGORITHM::RushLarsen(iss, this->var_[iSp], dt, tiSp);

    double ip = AiF*this->var_[iFp]+AiS*this->var_[iSp];
    double fItop = (1.0/(1.0+this->prm_[KmCaMK]/CaMKa));
    this->cur_[Ito] = (1.0-this->block_coeff_[Ito]) * (this->prm_[Gto]*(v_new-EK)*((1.0-fItop)*this->var_[a]*i+fItop*this->var_[ap]*ip));

    //// calculate ICaL, ICaNa, ICaK
    // gating NP 
    double dss = 1.0/(1.0+std::exp((-(v_new+3.940))/4.230));
    double td = 0.6+1.0/(std::exp(-0.05*(v_new+6.0))+std::exp(0.09*(v_new+14.0)));
    this->var_[d] = ALGORITHM::RushLarsen(dss, this->var_[d], dt, td);

    double fss = 1.0/(1.0+std::exp((v_new+19.58)/3.696));
    double tff = 7.0+1.0/(0.0045*std::exp(-(v_new+20.0)/10.0)+0.0045*std::exp((v_new+20.0)/10.0));
    double Aff = 0.6;
    this->var_[ff] = ALGORITHM::RushLarsen(fss, this->var_[ff], dt, tff);

    double tfs = 1000.0+1.0/(0.000035*std::exp(-(v_new+5.0)/4.0)+0.000035*std::exp((v_new+5.0)/6.0));
    double Afs = 1.0-Aff;
    this->var_[fs] = ALGORITHM::RushLarsen(fss, this->var_[fs], dt, tfs);

    this->prm_[f] = Aff*this->var_[ff]+Afs*this->var_[fs];
    double fcass = fss; 
    double tfcaf = 7.0+1.0/(0.04*std::exp(-(v_new-4.0)/7.0)+0.04*std::exp((v_new-4.0)/7.0));
    double Afcaf = 0.3+0.6/(1.0+std::exp((v_new-10.0)/10.0));
    this->var_[fcaf] = ALGORITHM::RushLarsen(fcass, this->var_[fcaf], dt, tfcaf);

    double tfcas = 100.0+1.0/(0.00012*std::exp(-v_new/3.0)+0.00012*std::exp(v_new/7.0));
    double Afcas = 1.0-Afcaf;
    this->var_[fcas] = ALGORITHM::RushLarsen(fcass, this->var_[fcas], dt, tfcas);

    double fca = Afcaf*this->var_[fcaf]+Afcas*this->var_[fcas];
    double tjca = 75.0;
    this->var_[jca] = ALGORITHM::RushLarsen(fcass, this->var_[jca], dt, tjca); 

    // gating CaMK-P
    double tffp = 2.5*tff; 
    this->var_[ffp] = ALGORITHM::RushLarsen(fss, this->var_[ffp], dt, tffp); 

    double fp = Aff*this->var_[ffp]+Afs*this->var_[fs];
    double tfcafp = 2.5*tfcaf; 
    this->var_[fcafp] = ALGORITHM::RushLarsen(fcass, this->var_[fcafp], dt, tfcaf); 

    double fcap = Afcaf*this->var_[fcafp] + Afcas*this->var_[fcas];
    double Kmn = 0.002;
    double k2n = 1000.0;
    double km2n = this->var_[jca];
    double anca = 1.0/(k2n/km2n + std::pow((1.0+Kmn/this->var_[cass]),4.0));
    double dnca = anca*k2n-this->var_[nca]*km2n;
    this->var_[nca] = ALGORITHM::ForwardEuler(this->var_[nca], dt, dnca); 

    double PhiCaNa = vffrt*(0.75*this->var_[nass]*std::exp(1.0*vfrt)-0.75*this->prm_[nao])/(std::exp(vfrt)-1.0);
    double PhiCaK = vffrt*(0.75*this->var_[kss]*std::exp(1.0*vfrt)-0.75*this->prm_[ko])/(std::exp(vfrt)-1.0);
    double zca = 2.0;
    double PCap = 1.1*this->prm_[PCa];
    double PCaNa = 0.00125*this->prm_[PCa];
    double PCaK = 3.574e-4*this->prm_[PCa];
    double PCaNap = 0.00125*PCap;
    double PCaKp = 3.574e-4*PCap;

    // gating PKA-P
    double PCaP = this->prm_[GCal_mul]*this->prm_[PCa];
    double PCaNaP = 0.00125*PCaP;
    double PCaKP = 3.574e-4*PCaP;

    // PKA-P new steady state and variables
    double dPss = 1.0/(1.0+std::exp((-(v_new+3.940+this->prm_[VP_d]))/4.230));
    this->var_[dP] = ALGORITHM::RushLarsen(dPss, this->var_[dP], dt, td); 

    double fPss = 1.0/(1.0+exp((v_new+19.58+this->prm_[VP_f])/3.696));
    this->var_[fPf] = ALGORITHM::RushLarsen(fPss, this->var_[fPf], dt, tff); 
    this->var_[fPs] = ALGORITHM::RushLarsen(fPss, this->var_[fPs], dt, tfs); 
    

    double fcaPss = fPss;
    this->var_[fcaPf] = ALGORITHM::RushLarsen(fcaPss, this->var_[fcaPf], dt, tfcaf); 
    this->var_[fcaPs] = ALGORITHM::RushLarsen(fcaPss, this->var_[fcaPs], dt, tfcas); 
    
    double fP = Aff*this->var_[fPf] + Afs*this->var_[fPs];
    double fcaP = Afcaf*this->var_[fcaPf] + Afcas*this->var_[fcaPs];

    // Both-P population takes on dP, for this->prm_[f] gate, takes ss from PKA and tau from CaMK (only fast component was modified)
    double fBPss = fPss;
    this->var_[fBPf] = ALGORITHM::RushLarsen(fBPss, this->var_[fBPf], dt, tffp);

    double fBP = Aff*this->var_[fBPf]+Afs*this->var_[fPs];
    double fcaBPss = fcaPss;
    this->var_[fcaBPf] = ALGORITHM::RushLarsen(fcaBPss, this->var_[fcaBPf], dt, tfcafp);

    double fcaBP = Afcaf*this->var_[fcaBPf]+Afcas*this->var_[fcaPs];
    double fICaLp=(1.0/(1.0+this->prm_[KmCaMK]/CaMKa)); 
    double fICaL_P = ICaL_P; 
    double fICaL_BP = fICaLp*fICaL_P;
    double fICaL_CaMKonly = fICaLp-fICaL_BP;
    double fICaL_PKAonly = fICaL_P-fICaL_BP;

    double Cass_ceiling = 0.025;
    double PhiCaL = 4.0*vffrt*(this->var_[cass]*std::exp(2.0*vfrt)-0.341*this->prm_[cao])/(std::exp(2.0*vfrt)-1.0);
    if (this->var_[cass] > Cass_ceiling) {
       PhiCaL = 4.0*vffrt*(Cass_ceiling*std::exp(2.0*vfrt)-0.341*this->prm_[cao])/(std::exp(2.0*vfrt)-1.0);
    }
    double ICaL_NP = this->prm_[PCa]*PhiCaL*this->var_[d]*(this->prm_[f]*(1.0-this->var_[nca])+this->var_[jca]*fca*this->var_[nca]);
    double ICaL_CaMK = PCap*PhiCaL*this->var_[d]*(fp*(1.0-this->var_[nca])+this->var_[jca]*fcap*this->var_[nca]);
    double ICaL_PKA = PCaP*PhiCaL*this->var_[dP]*(fP*(1.0-this->var_[nca])+this->var_[jca]*fcaP*this->var_[nca]);
    double ICaL_BP = PCaP*PhiCaL*this->var_[dP]*(fBP*(1.0-this->var_[nca])+this->var_[jca]*fcaBP*this->var_[nca]);  

    double ICaNa_NP = PCaNa*PhiCaNa*this->var_[d]*(this->prm_[f]*(1.0-this->var_[nca])+this->var_[jca]*fca*this->var_[nca]);
    double ICaNa_CaMK = PCaNap*PhiCaNa*this->var_[d]*(fp*(1.0-this->var_[nca])+this->var_[jca]*fcap*this->var_[nca]);
    double ICaNa_PKA = PCaNaP*PhiCaNa*this->var_[dP]*(fP*(1.0-this->var_[nca])+this->var_[jca]*fcaP*this->var_[nca]);
    double ICaNa_BP = PCaNaP*PhiCaNa*this->var_[dP]*(fBP*(1.0-this->var_[nca])+this->var_[jca]*fcaBP*this->var_[nca]);

    double ICaK_NP = PCaK*PhiCaK*this->var_[d]*(this->prm_[f]*(1.0-this->var_[nca])+this->var_[jca]*fca*this->var_[nca]);
    double ICaK_CaMK = PCaKp*PhiCaK*this->var_[d]*(fp*(1.0-this->var_[nca])+this->var_[jca]*fcap*this->var_[nca]);
    double ICaK_PKA = PCaKP*PhiCaK*this->var_[dP]*(fP*(1.0-this->var_[nca])+this->var_[jca]*fcaP*this->var_[nca]);
    double ICaK_BP = PCaKP*PhiCaK*this->var_[dP]*(fBP*(1.0-this->var_[nca])+this->var_[jca]*fcaBP*this->var_[nca]);

    // 4 population
    this->cur_[ICaL]  = (1.0-this->block_coeff_[ICaL]) * ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaL_NP + fICaL_CaMKonly*ICaL_CaMK + fICaL_PKAonly*ICaL_PKA + fICaL_BP*ICaL_BP);
    this->cur_[ICaNa] = (1.0-this->block_coeff_[ICaNa]) * ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaNa_NP + fICaL_CaMKonly*ICaNa_CaMK + fICaL_PKAonly*ICaNa_PKA + fICaL_BP*ICaNa_BP);
    this->cur_[ICaK]  = (1.0-this->block_coeff_[ICaK]) * ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaK_NP + fICaL_CaMKonly*ICaK_CaMK + fICaL_PKAonly*ICaK_PKA + fICaL_BP*ICaK_BP);

    //// calculate IKr
    double xrss = 1.0/(1.0+std::exp((-(v_new+8.337))/6.789));
    double txrf = 12.98+1.0/(0.3652*std::exp((v_new-31.66)/3.869)+4.123e-5*std::exp((-(v_new-47.78))/20.38));
    double Axrf = 1.0/(1.0+std::exp((v_new+54.81)/38.21));
    this->var_[xrf] = ALGORITHM::RushLarsen(xrss, this->var_[xrf], dt, txrf);
    
    double txrs = 1.865+1.0/(0.06629*std::exp((v_new-34.70)/7.355)+1.128e-5*std::exp((-(v_new-29.74))/25.94));
    double Axrs = 1.0-Axrf;
    double dxrs = (xrss-this->var_[xrs])/txrs;
    this->var_[xrs] = ALGORITHM::ForwardEuler(this->var_[xrs], dt, dxrs); 

    double xr = Axrf*this->var_[xrf] + Axrs*this->var_[xrs];
    double rkr = 1.0/(1.0+std::exp((v_new+55.0)/75.0))*1.0/(1.0+std::exp((v_new-10.0)/30.0));
    this->cur_[IKr] = (1.0-this->block_coeff_[IKr]) * (this->prm_[GKr]*std::sqrt(this->prm_[ko]/5.4)*xr*rkr*(v_new-EK));

    //// calculate IKs
    double xs1ss = 1.0/(1.0+std::exp((-(v_new+11.60))/8.932));
    double txs1 = 817.3+1.0/(2.326e-4*std::exp((v_new+48.28)/17.80)+0.001292*std::exp((-(v_new+210.0))/230.0));
    this->var_[xs1] = ALGORITHM::RushLarsen(xs1ss, this->var_[xs1], dt, txs1);

    double xs2ss = xs1ss;
    double txs2 = 1.0/(0.01*std::exp((v_new-50.0)/20.0)+0.0193*std::exp((-(v_new+66.54))/31.0));
    this->var_[xs2] = ALGORITHM::RushLarsen(xs2ss, this->var_[xs2], dt, txs2);

    double KsCa = 1.0+0.6/(1.0+std::pow((3.8e-5/this->var_[cai]),1.4));

    // gating PKA-P
    double GKsP = this->prm_[GKS_mul]*this->prm_[GKs];
    double VP_xs = 0.;
    double raw_taux2 = 1.;
    double taux1factor = 2.75;  
    double raw_taux1 = (1.0-taux1factor) / (2.326e-4*std::exp((10.+48.28)/17.80)+0.001292*std::exp((-(10+210.0))/230.0));

    double xs1Pss = 1.0/(1.0+std::exp((-(v_new+VP_xs+11.60))/8.932));
    double txs1P = raw_taux1 + 817.3 + taux1factor / ((2.326*std::pow(10,-4)*std::exp((v_new+48.28)/17.8)) + 0.001292*std::exp(-(v_new+210.0)/230.0)) ;
    this->var_[xs1P] = ALGORITHM::RushLarsen(xs1Pss, this->var_[xs1P], dt, txs1P);

    double xs2Pss = xs1Pss; 
    double txs2P = txs2*raw_taux2 ;
    this->var_[xs2P] = ALGORITHM::RushLarsen(xs2Pss, this->var_[xs2P], dt, txs2P);

    // no CaMK-P in IKs
    double IKs_NP = this->prm_[GKs]*KsCa*this->var_[xs1]*this->var_[xs2]*(v_new-EKs);
    double IKs_PKA = GKsP*KsCa*this->var_[xs1P]*this->var_[xs2P]*(v_new-EKs);
    double fIKs_PKA = IKs_P;
    this->cur_[IKs] = (1.0-this->block_coeff_[IKs]) * ((1.0-fIKs_PKA)*IKs_NP + fIKs_PKA*IKs_PKA);

    //// calculate IK1
    double xk1ss = 1.0/(1.0+std::exp(-(v_new+2.5538*this->prm_[ko]+144.59)/(1.5692*this->prm_[ko]+3.8115)));
    double txk1 = 122.2/(std::exp((-(v_new+127.2))/20.36)+std::exp((v_new+236.8)/69.33));
    this->var_[xk1] = ALGORITHM::RushLarsen(xk1ss, this->var_[xk1], dt, txk1);

    double rk1 = 1.0/(1.0+std::exp((v_new+105.8-2.6*this->prm_[ko])/9.493));
    this->cur_[IK1] = (1.0-this->block_coeff_[IK1]) * (this->prm_[GK1]*std::sqrt(this->prm_[ko])*rk1*this->var_[xk1]*(v_new-EK));    

    //// calculate INaCa_i
    double kna1 = 15.0;
    double kna2 = 5.0;
    double kna3 = 88.12;
    double kasymm = 12.5;
    double wna = 6.0e4;
    double wca = 6.0e4;
    double wnaca = 5.0e3;
    double kcaon = 1.5e6;
    double kcaoff = 5.0e3;
    double qna = 0.5224;
    double qca = 0.1670;
    double hca = std::exp((qca*v_new*this->prm_[F])/(this->prm_[R]*this->prm_[T]));
    double hna = std::exp((qna*v_new*this->prm_[F])/(this->prm_[R]*this->prm_[T]));
    double h1 = 1+this->var_[nai]/kna3*(1+hna);
    double h2 = (this->var_[nai]*hna)/(kna3*h1);
    double h3 = 1.0/h1;
    double h4 = 1.0+this->var_[nai]/kna1*(1.+this->var_[nai]/kna2);
    double h5 = this->var_[nai]*this->var_[nai]/(h4*kna1*kna2);
    double h6 = 1.0/h4;
    double h7 = 1.0+this->prm_[nao]/kna3*(1.0+1.0/hna);
    double h8 = this->prm_[nao]/(kna3*hna*h7);
    double h9 = 1.0/h7;
    double h10 = kasymm+1.0+this->prm_[nao]/kna1*(1.0+this->prm_[nao]/kna2);
    double h11 = this->prm_[nao]*this->prm_[nao]/(h10*kna1*kna2);
    double h12 = 1.0/h10;
    double k1 = h12*this->prm_[cao]*kcaon;
    double k2 = kcaoff;
    double k3p = h9*wca;
    double k3pp = h8*wnaca;
    double k3 = k3p+k3pp;
    double k4p = h3*wca/hca;
    double k4pp = h2*wnaca;
    double k4 = k4p+k4pp;
    double k5 = kcaoff;
    double k6 = h6*this->var_[cai]*kcaon;
    double k7 = h5*h2*wna;
    double k8 = h8*h11*wna;
    double x1 = k2*k4*(k7+k6)+k5*k7*(k2+k3);
    double x2 = k1*k7*(k4+k5)+k4*k6*(k1+k8);
    double x3 = k1*k3*(k7+k6)+k8*k6*(k2+k3);
    double x4 = k2*k8*(k4+k5)+k3*k5*(k1+k8);
    double E1 = x1/(x1+x2+x3+x4);
    double E2 = x2/(x1+x2+x3+x4);
    double E3 = x3/(x1+x2+x3+x4);
    double E4 = x4/(x1+x2+x3+x4);
    double KmCaAct = 150.0e-6;
    double allo = 1.0/(1.0+std::pow((KmCaAct/this->var_[cai]),2.0));
    double zna = 1.0;
    double JncxNa = 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    double JncxCa = E2*k2-E1*k1;
    this->cur_[INaCa_i] = (1.0-this->block_coeff_[INaCa_i]) * (0.8*this->prm_[Gncx]*allo*(zna*JncxNa+zca*JncxCa));

    //// calculate INaCa_ss
    h1 = 1.+this->var_[nass]/kna3*(1+hna);
    h2 = (this->var_[nass]*hna)/(kna3*h1);
    h3 = 1.0/h1;
    h4 = 1.0+this->var_[nass]/kna1*(1.+this->var_[nass]/kna2);
    h5 = this->var_[nass]*this->var_[nass]/(h4*kna1*kna2);
    h6 = 1.0/h4;
    h7 = 1.0+this->prm_[nao]/kna3*(1.0+1.0/hna);
    h8 = this->prm_[nao]/(kna3*hna*h7);
    h9 = 1.0/h7;
    h10 = kasymm+1.0+this->prm_[nao]/kna1*(1+this->prm_[nao]/kna2);
    h11 = this->prm_[nao]*this->prm_[nao]/(h10*kna1*kna2);
    h12 = 1.0/h10;
    k1 = h12*this->prm_[cao]*kcaon;
    k2 = kcaoff;
    k3p = h9*wca;
    k3pp = h8*wnaca;
    k3 = k3p+k3pp;
    k4p = h3*wca/hca;
    k4pp = h2*wnaca;
    k4 = k4p+k4pp;
    k5 = kcaoff;
    k6 = h6*this->var_[cass]*kcaon;
    k7 = h5*h2*wna;
    k8 = h8*h11*wna;
    x1 = k2*k4*(k7+k6)+k5*k7*(k2+k3);
    x2 = k1*k7*(k4+k5)+k4*k6*(k1+k8);
    x3 = k1*k3*(k7+k6)+k8*k6*(k2+k3);
    x4 = k2*k8*(k4+k5)+k3*k5*(k1+k8);
    E1 = x1/(x1+x2+x3+x4);
    E2 = x2/(x1+x2+x3+x4);
    E3 = x3/(x1+x2+x3+x4);
    E4 = x4/(x1+x2+x3+x4);
    KmCaAct = 150.0e-6;
    allo = 1.0/(1.0+std::pow((KmCaAct/this->var_[cass]),2.0));
    JncxNa = 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa = E2*k2-E1*k1;
    this->cur_[INaCa_ss] = (1.0-this->block_coeff_[INaCa_ss]) * (0.2*this->prm_[Gncx]*allo*(zna*JncxNa+zca*JncxCa));

    
    //// calculate INaK
    double k1p = 949.5;
    double k1m = 182.4;
    double k2p = 687.2;
    double k2m = 39.4;
    k3p = 1899.0;
    double k3m = 79300.0;
    k4p = 639.0;
    double k4m = 40.0;

    double Knai0_np = 9.073;
    double Knai0_P = Knai0_np*0.7;
    double Knai0 = 9.073;
    double fINaK_PKA = INaK_P;
    if (fINaK_PKA > 0)
       Knai0 = (1.0-fINaK_PKA)*Knai0_np + fINaK_PKA*Knai0_P; 

    double Knao0 = 27.78;
    double delta = -0.1550;
    double Knai = Knai0*std::exp((delta*v_new*this->prm_[F])/(3.0*this->prm_[R]*this->prm_[T]));
    double Knao = Knao0*std::exp(((1.0-delta)*v_new*this->prm_[F])/(3.0*this->prm_[R]*this->prm_[T]));
    double Kki = 0.5;
    double Kko = 0.3582;
    double MgADP = 0.05;
    double MgATP = 9.8;
    double Kmgatp = 1.698e-7;
    double H = 1.0e-7;
    double eP = 4.2;
    double Khp = 1.698e-7;
    double Knap = 224.0;
    double Kxkur = 292.0;

    double P = eP/(1.0+H/Khp+this->var_[nai]/Knap+this->var_[ki]/Kxkur);
    double a1 = (k1p*std::pow((this->var_[nai]/Knai),3.0))/(std::pow((1.0+this->var_[nai]/Knai),3.0)+std::pow((1.0+this->var_[ki]/Kki),2.0)-1.0);
    double b1 = k1m*MgADP;
    double a2 = k2p;
    double b2 = (k2m*std::pow((this->prm_[nao]/Knao),3.0))/(std::pow((1.0+this->prm_[nao]/Knao),3.0)+std::pow((1.0+this->prm_[ko]/Kko),2.0)-1.0);
    double a3 = (k3p*std::pow((this->prm_[ko]/Kko),2.0))/(std::pow((1.0+this->prm_[nao]/Knao),3.0)+std::pow((1.0+this->prm_[ko]/Kko),2.0)-1.0);
    double b3 = (k3m*P*H)/(1.0+MgATP/Kmgatp);
    double a4 = (k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
    double b4 = (k4m*std::pow((this->var_[ki]/Kki),2.0))/(std::pow((1.0+this->var_[nai]/Knai),3.0)+std::pow((1.0+this->var_[ki]/Kki),2.0)-1.0);
    x1 = a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
    x2 = b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
    x3 = a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
    x4 = b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
    E1 = x1/(x1+x2+x3+x4);
    E2 = x2/(x1+x2+x3+x4);
    E3 = x3/(x1+x2+x3+x4);
    E4 = x4/(x1+x2+x3+x4);
    double zk = 1.0;
    double JnakNa = 3.0*(E1*a3-E2*b3);
    double JnakK = 2.0*(E4*b1-E3*a1);

    // short cut here cause INaK_P will only be given 0 or 1
    this->cur_[INaK] = (1.0-this->block_coeff_[INaK]) * (this->prm_[Pnak]*(zna*JnakNa+zk*JnakK));

    //// background current and membrane pumps
    //calculate IKb
    double xkb =1.0/(1.0+std::exp(-(v_new-14.48)/18.34));

    // added PKA effects according to O'Hara 2011
    double IKb_NP = this->prm_[GKb]*xkb*(v_new-EK);
    double GKbP = this->prm_[GKb]*2.5 ;
    double IKb_PKA = GKbP*xkb*(v_new-EK);
    double fIKb_P = IKb_P;
    this->cur_[IKb] = (1.0-this->block_coeff_[IKb]) * ((1.0-fIKb_P)*IKb_NP + fIKb_P*IKb_PKA);
    
    //calculate INab
    this->cur_[INab] = (1.0-this->block_coeff_[INab]) * (this->prm_[PNab]*vffrt*(this->var_[nai]*std::exp(vfrt)-this->prm_[nao])/(std::exp(vfrt)-1.0));

    //calculate ICab
    this->cur_[ICab] = (1.0-this->block_coeff_[ICab]) * (this->prm_[PCab]*4.0*vffrt*(this->var_[cai]*std::exp(2.0*vfrt)-0.341*this->prm_[cao])/(std::exp(2.0*vfrt)-1.0));

    //calculate IpCa
    this->cur_[IpCa] = (1.0-this->block_coeff_[IpCa]) * (this->prm_[GpCa]*this->var_[cai]/(0.0005+this->var_[cai]));

    // Compute the total Iion current.
    this->cur_[Gng20Cur::Iion] = this->cur_[INa]  + this->cur_[INaL] +  this->cur_[Ito] + this->cur_[ICaL] + this->cur_[ICaNa] + 
                               this->cur_[ICaK] + this->cur_[IKr]  +  this->cur_[IKs] + this->cur_[IK1]  + this->cur_[INaCa_i] +
                               this->cur_[INaCa_ss] + this->cur_[INaK] + this->cur_[INab] + this->cur_[IKb]  + this->cur_[IpCa] + this->cur_[ICab];

    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[Gng20Cur::Iion] - stim_current);

    //// calculate diffusion fluxes
    double JdiffNa = (this->var_[nass]-this->var_[nai])/2.0;
    double JdiffK = (this->var_[kss]-this->var_[ki])/2.0;
    double Jdiff = (this->var_[cass]-this->var_[cai])/0.2;

    //// calculate ryanodione receptor calcium induced calcium release from the jsr

    // NP channel behavior
    double bt = 4.75;
    double a_rel = 0.5*bt;
    double ICaLdivCajsr = (-this->cur_[ICaL])/(1.0+std::pow((1.5/this->var_[cajsr]),8.0));
    double Jrel_inf = a_rel * ICaLdivCajsr; 
    if (static_cast<int>(this->prm_[celltype]) == 2) Jrel_inf *= 1.7;
    double tau_rel = bt/(1.0+0.0123/this->var_[cajsr]);
    if (tau_rel < 0.001) tau_rel=0.001; 
    this->var_[Jrelnp] = ALGORITHM::RushLarsen(Jrel_inf, this->var_[Jrelnp], dt, tau_rel); 

    // CaMK-P channel behavior
    double btp = 1.25*bt;
    double a_relp = 0.5*btp;
    double Jrel_infp = a_relp * ICaLdivCajsr;
    if (static_cast<int>(this->prm_[celltype]) == 2) Jrel_infp *= 1.7;
    double tau_relp = btp/(1.0+0.0123/this->var_[cajsr]);
    if (tau_relp < 0.001) tau_relp=0.001; 
    this->var_[Jrelp] = ALGORITHM::RushLarsen(Jrel_infp, this->var_[Jrelp], dt, tau_rel); 

    // PKA-P channel behavior
    double a_relP = a_rel*2.5; 
    double Jrel_infP = a_relP * ICaLdivCajsr;
    if (static_cast<int>(this->prm_[celltype]) == 2) Jrel_infP *= 1.7;
    double tau_relP = tau_rel*0.75;
    if (tau_relP<0.001) tau_relP=0.001; 
    this->var_[JrelP] = ALGORITHM::RushLarsen(Jrel_infP, this->var_[JrelP], dt, tau_relP); 

    // Both-P
    double a_relBP = a_relp*2.5; 
    double Jrel_infBP = a_relBP * ICaLdivCajsr;
    if (static_cast<int>(this->prm_[celltype]) == 2) Jrel_infBP *= 1.7;
    double tau_relBP = tau_relp*0.75;
    if (tau_relBP < 0.001) tau_relBP = 0.001; 
    this->var_[JrelBP] = ALGORITHM::RushLarsen(Jrel_infBP, this->var_[JrelBP], dt, tau_relBP); 

    // 4 population, with both-P population takes on PKA behavior    
    double fJrelp=(1.0/(1.0+this->prm_[KmCaMK]/CaMKa));
    double fJrel_PKA = RyR_P;
    double fJrel_BP = fJrelp*fJrel_PKA;
    double fJrel_CaMKonly = fJrelp - fJrel_BP;
    double fJrel_PKAonly = fJrel_PKA - fJrel_BP;
    double Jrel = this->prm_[RyR_total]* ((1.0-fJrel_CaMKonly-fJrel_PKAonly-fJrel_BP)*this->var_[Jrelnp] + fJrel_CaMKonly*this->var_[Jrelp] + fJrel_PKAonly*this->var_[JrelP] + fJrel_BP*this->var_[JrelBP]);

    //// calculate serca pump, ca uptake flux
    // NP SERCA
    double Jupnp = 0.004375*this->var_[cai]/(this->var_[cai]+0.00092);
    
    // CaMK-P SERCA 
    double Jupp = 2.75*0.004375*this->var_[cai]/(this->var_[cai]+0.00092-0.00017);

    // PKA-P SERCA
    double mul_PKA = 0.54;
    double JupP = 0.004375*this->var_[cai]/(this->var_[cai]+0.00092*mul_PKA);
    double Jup_BP = 2.75*0.004375*this->var_[cai]/(this->var_[cai]+(0.00092-0.00017)*mul_PKA);

    if (static_cast<int>(this->prm_[celltype]) == 1) {
        Jupnp *= 1.3;
        Jupp *= 1.3;
        JupP *= 1.3;
        Jup_BP *= 1.3;
    }

    // 4 population
    double fJupp = (1.0/(1.0+this->prm_[KmCaMK]/CaMKa)); // CaMK-P channel fraction
    
    // fJupp = fCaMK ;
    double fJup_P = SERCA_P; 
    double fJup_BP = fJupp*fJup_P;
    double fJup_CaMKonly = fJupp - fJup_BP;
    double fJup_PKAonly = fJup_P - fJup_BP;


    double Jleak = this->prm_[Leak_total]*(0.0039375*this->var_[cansr]/15.0);
    double Jup = this->prm_[SERCA_total]*((1.0-fJup_CaMKonly-fJup_PKAonly-fJup_BP)*Jupnp + fJup_CaMKonly*Jupp + fJup_PKAonly*JupP + fJup_BP*Jup_BP) - Jleak;

    //calculate tranlocation flux
    double Jtr = this->prm_[Trans_total]*((this->var_[cansr]-this->var_[cajsr])/100.0);

    //// calcium buffer constants
    double cmdnmax = 0.05;
    if (static_cast<int>(this->prm_[celltype]) == 1) cmdnmax *= 1.3;
    
    double kmcmdn = 0.00238;
    double trpnmax = 0.07;
    double kmtrpn_np = 0.0005;
    double kmtrpn_P = kmtrpn_np*1.6 ;
    double kmtrpn = 0.0005;
    double fTnI_PKA = TnI_P;
    if (fTnI_PKA > 0) kmtrpn = (1-fTnI_PKA)*kmtrpn_np + fTnI_PKA*kmtrpn_P;
    
    double BSRmax = 0.047;
    double KmBSR = 0.00087;
    double BSLmax = 1.124;
    double KmBSL = 0.0087;
    double csqnmax = 10.0;
    double kmcsqn = 0.8;

    //update intracellular concentrations, using buffers for cai, cass, cajsr
    double dnai = -(this->cur_[INa]+this->cur_[INaL]+3.0*this->cur_[INaCa_i]+3.0*this->cur_[INaK]+this->cur_[INab])*this->prm_[Acap]/(this->prm_[F]*this->prm_[vmyo])+JdiffNa*this->prm_[vss]/this->prm_[vmyo];
    this->var_[nai] = ALGORITHM::ForwardEuler(this->var_[nai], dt, dnai);

    double dnass = -(this->cur_[ICaNa]+3.0*this->cur_[INaCa_ss])*this->prm_[Acap]/(this->prm_[F]*this->prm_[vss])-JdiffNa;
    this->var_[nass] = ALGORITHM::ForwardEuler(this->var_[nass], dt, dnass); 

    double dki = -(this->cur_[Ito]+this->cur_[IKr]+this->cur_[IKs]+this->cur_[IK1]+this->cur_[IKb]-stim_current-2.0*this->cur_[INaK])*this->prm_[Acap]/(this->prm_[F]*this->prm_[vmyo])+JdiffK*this->prm_[vss]/this->prm_[vmyo];
    this->var_[ki] = ALGORITHM::ForwardEuler(this->var_[ki], dt, dki); 

    double dkss = -(this->cur_[ICaK])*this->prm_[Acap]/(this->prm_[F]*this->prm_[vss])-JdiffK;
    this->var_[kss] = ALGORITHM::ForwardEuler(this->var_[kss], dt, dkss); 

    double Bcai = 1.0/(1.0+cmdnmax*kmcmdn/std::pow((kmcmdn+this->var_[cai]),2.0)+trpnmax*kmtrpn/std::pow((kmtrpn+this->var_[cai]),2.0));
    double dcai = Bcai*(-(this->cur_[IpCa]+this->cur_[ICab]-2.0*this->cur_[INaCa_i])*this->prm_[Acap]/(2.0*this->prm_[F]*this->prm_[vmyo])-Jup*this->prm_[vnsr]/this->prm_[vmyo]+Jdiff*this->prm_[vss]/this->prm_[vmyo]);
    this->var_[cai] = ALGORITHM::ForwardEuler(this->var_[cai], dt, dcai); 

    double Bcass = 1.0/(1.0+BSRmax*KmBSR/std::pow((KmBSR+this->var_[cass]),2.0)+BSLmax*KmBSL/std::pow((KmBSL+this->var_[cass]),2.0));
    double dcass = Bcass*(-(this->cur_[ICaL]-2.0*this->cur_[INaCa_ss])*this->prm_[Acap]/(2.0*this->prm_[F]*this->prm_[vss])+Jrel*this->prm_[vjsr]/this->prm_[vss]-Jdiff);
    this->var_[cass] = ALGORITHM::ForwardEuler(this->var_[cass], dt, dcass); 

    double dcansr = Jup-Jtr*this->prm_[vjsr]/this->prm_[vnsr];
    this->var_[cansr] = ALGORITHM::ForwardEuler(this->var_[cansr], dt, dcansr); 
    
    double Bcajsr = 1.0/(1.0+csqnmax*kmcsqn/std::pow((kmcsqn+this->var_[cajsr]),2.0));
    double dcajsr = Bcajsr*(Jtr-Jrel);
    this->var_[cajsr] = ALGORITHM::ForwardEuler(this->var_[cajsr], dt, dcajsr); 

}


std::string Gong2020::PrintVariables() const 
{
    using namespace Gng20Var;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "nai = " << this->var_[nai] << "\n";
    oss << "nass = " << this->var_[nass] << "\n";
    oss << "ki = " << this->var_[ki] << "\n";
    oss << "kss = " << this->var_[kss] << "\n";
    oss << "cai = " << this->var_[cai] << "\n";
    oss << "cass = " << this->var_[cass] << "\n";
    oss << "cansr = " << this->var_[cansr] << "\n";
    oss << "cajsr = " << this->var_[cajsr] << "\n";
    oss << "m = " << this->var_[m] << "\n";
    oss << "hf = " << this->var_[hf] << "\n";
    oss << "hs = " << this->var_[hs] << "\n";
    oss << "j = " << this->var_[j] << "\n";
    oss << "hsp = " << this->var_[hsp] << "\n";
    oss << "jp = " << this->var_[jp] << "\n";
    oss << "mL = " << this->var_[mL] << "\n";
    oss << "hL = " << this->var_[hL] << "\n";
    oss << "hLp = " << this->var_[hLp] << "\n";
    oss << "a = " << this->var_[a] << "\n";
    oss << "iF = " << this->var_[iF] << "\n";
    oss << "iS = " << this->var_[iS] << "\n";
    oss << "ap = " << this->var_[ap] << "\n";
    oss << "iFp = " << this->var_[iFp] << "\n";
    oss << "iSp = " << this->var_[iSp] << "\n";
    oss << "d = " << this->var_[d] << "\n";
    oss << "ff = " << this->var_[ff] << "\n";
    oss << "fs = " << this->var_[fs] << "\n";
    oss << "fcaf = " << this->var_[fcaf] << "\n";
    oss << "fcas = " << this->var_[fcas] << "\n";
    oss << "jca = " << this->var_[jca] << "\n";
    oss << "nca = " << this->var_[nca] << "\n";
    oss << "ffp = " << this->var_[ffp] << "\n";
    oss << "fcafp = " << this->var_[fcafp] << "\n";
    oss << "xrf = " << this->var_[xrf] << "\n";
    oss << "xrs = " << this->var_[xrs] << "\n";
    oss << "xs1 = " << this->var_[xs1] << "\n";
    oss << "xs2 = " << this->var_[xs2] << "\n";
    oss << "xk1 = " << this->var_[xk1] << "\n";
    oss << "Jrelnp = " << this->var_[Jrelnp] << "\n";
    oss << "Jrelp = " << this->var_[Jrelp] << "\n";
    oss << "CaMKt = " << this->var_[CaMKt] << "\n";
    oss << "hPf = " << this->var_[hPf] << "\n";
    oss << "hPs = " << this->var_[hPs] << "\n";
    oss << "jP = " << this->var_[jP] << "\n";
    oss << "hBPf = " << this->var_[hBPf] << "\n";
    oss << "hBPs = " << this->var_[hBPs] << "\n";
    oss << "jBP = " << this->var_[jBP] << "\n";
    oss << "dP = " << this->var_[dP] << "\n";
    oss << "fPf = " << this->var_[fPf] << "\n";
    oss << "fPs = " << this->var_[fPs] << "\n";
    oss << "fcaPf = " << this->var_[fcaPf] << "\n";
    oss << "fcaPs = " << this->var_[fcaPs] << "\n";
    oss << "fBPf = " << this->var_[fBPf] << "\n";
    oss << "fcaBPf = " << this->var_[fcaBPf] << "\n";
    oss << "xs1P = " << this->var_[xs1P] << "\n";
    oss << "xs2P = " << this->var_[xs2P] << "\n";
    oss << "JrelP = " << this->var_[JrelP] << "\n";
    oss << "JrelBP = " << this->var_[JrelBP] << "\n";
    oss << "beta_cav_Gs_aGTP = " << this->var_[beta_cav_Gs_aGTP] << "\n";
    oss << "beta_eca_Gs_aGTP = " << this->var_[beta_eca_Gs_aGTP] << "\n";
    oss << "beta_cyt_Gs_aGTP = " << this->var_[beta_cyt_Gs_aGTP] << "\n";
    oss << "beta_cav_Gs_bg = " << this->var_[beta_cav_Gs_bg] << "\n";
    oss << "beta_eca_Gs_bg = " << this->var_[beta_eca_Gs_bg] << "\n";
    oss << "beta_cyt_Gs_bg = " << this->var_[beta_cyt_Gs_bg] << "\n";
    oss << "beta_cav_Gs_aGDP = " << this->var_[beta_cav_Gs_aGDP] << "\n";
    oss << "beta_eca_Gs_aGDP = " << this->var_[beta_eca_Gs_aGDP] << "\n";
    oss << "beta_cyt_Gs_aGDP = " << this->var_[beta_cyt_Gs_aGDP] << "\n";
    oss << "cAMP_cav = " << this->var_[cAMP_cav] << "\n";
    oss << "cAMP_eca = " << this->var_[cAMP_eca] << "\n";
    oss << "cAMP_cyt = " << this->var_[cAMP_cyt] << "\n";
    oss << "beta_cav_Rb1_pka_tot = " << this->var_[beta_cav_Rb1_pka_tot] << "\n";
    oss << "beta_eca_Rb1_pka_tot = " << this->var_[beta_eca_Rb1_pka_tot] << "\n";
    oss << "beta_cyt_Rb1_pka_tot = " << this->var_[beta_cyt_Rb1_pka_tot] << "\n";
    oss << "beta_cav_Rb1_grk_tot = " << this->var_[beta_cav_Rb1_grk_tot] << "\n";
    oss << "beta_eca_Rb1_grk_tot = " << this->var_[beta_eca_Rb1_grk_tot] << "\n";
    oss << "beta_cyt_Rb1_grk_tot = " << this->var_[beta_cyt_Rb1_grk_tot] << "\n";
    oss << "pka_cav_ARC = " << this->var_[pka_cav_ARC] << "\n";
    oss << "pka_cav_A2RC = " << this->var_[pka_cav_A2RC] << "\n";
    oss << "pka_cav_A2R = " << this->var_[pka_cav_A2R] << "\n";
    oss << "pka_cav_C = " << this->var_[pka_cav_C] << "\n";
    oss << "pka_cav_PKIC = " << this->var_[pka_cav_PKIC] << "\n";
    oss << "pka_eca_ARC = " << this->var_[pka_eca_ARC] << "\n";
    oss << "pka_eca_A2RC = " << this->var_[pka_eca_A2RC] << "\n";
    oss << "pka_eca_A2R = " << this->var_[pka_eca_A2R] << "\n";
    oss << "pka_eca_C = " << this->var_[pka_eca_C] << "\n";
    oss << "pka_eca_PKIC = " << this->var_[pka_eca_PKIC] << "\n";
    oss << "pka_cyt_ARC = " << this->var_[pka_cyt_ARC] << "\n";
    oss << "pka_cyt_A2RC = " << this->var_[pka_cyt_A2RC] << "\n";
    oss << "pka_cyt_A2R = " << this->var_[pka_cyt_A2R] << "\n";
    oss << "pka_cyt_C = " << this->var_[pka_cyt_C] << "\n";
    oss << "pka_cyt_PKIC = " << this->var_[pka_cyt_PKIC] << "\n";
    oss << "PDE3_P_cav = " << this->var_[PDE3_P_cav] << "\n";
    oss << "PDE3_P_cyt = " << this->var_[PDE3_P_cyt] << "\n";
    oss << "PDE4_P_cav = " << this->var_[PDE4_P_cav] << "\n";
    oss << "PDE4_P_eca = " << this->var_[PDE4_P_eca] << "\n";
    oss << "PDE4_P_cyt = " << this->var_[PDE4_P_cyt] << "\n";
    oss << "inhib1_p = " << this->var_[inhib1_p] << "\n";
    oss << "ICaLp = " << this->var_[ICaLp] << "\n";
    oss << "IKsp = " << this->var_[IKsp] << "\n";
    oss << "iup_f_plb = " << this->var_[iup_f_plb] << "\n";
    oss << "f_tni = " << this->var_[f_tni] << "\n";
    oss << "ina_f_ina = " << this->var_[ina_f_ina] << "\n";
    oss << "f_inak = " << this->var_[f_inak] << "\n";
    oss << "RyRp = " << this->var_[RyRp] << "\n";
    oss << "f_ikur = " << this->var_[f_ikur] << "\n";
    oss << "beta_cav_Rb2_pka_tot = " << this->var_[beta_cav_Rb2_pka_tot] << "\n";
    oss << "beta_cav_Rb2_grk_tot = " << this->var_[beta_cav_Rb2_grk_tot] << "\n";
    oss << "beta_cav_Gi_aGTP = " << this->var_[beta_cav_Gi_aGTP] << "\n";
    oss << "beta_cav_Gi_bg = " << this->var_[beta_cav_Gi_bg] << "\n";
    oss << "beta_cav_Gi_aGDP = " << this->var_[beta_cav_Gi_aGDP] << "\n";
    oss << "beta_eca_Rb2_pka_tot = " << this->var_[beta_eca_Rb2_pka_tot] << "\n";
    oss << "beta_eca_Rb2_grk_tot = " << this->var_[beta_eca_Rb2_grk_tot] << "\n";
    oss << "beta_eca_Gi_aGTP = " << this->var_[beta_eca_Gi_aGTP] << "\n";
    oss << "beta_eca_Gi_bg = " << this->var_[beta_eca_Gi_bg] << "\n";
    oss << "beta_eca_Gi_aGDP = " << this->var_[beta_eca_Gi_aGDP];
    return oss.str();

}


std::string Gong2020::PrintParameters() const 
{
    using namespace Gng20Prm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "signaling = " << this->prm_[signaling] << "\n";
    oss << "electrophys = " << this->prm_[electrophys] << "\n";
    oss << "cipa = " << this->prm_[cipa] << "\n";
    oss << "GKS_mul = " << this->prm_[GKS_mul] << "\n";
    oss << "VP_d = " << this->prm_[VP_d] << "\n";
    oss << "VP_f = " << this->prm_[VP_f] << "\n";
    oss << "GCal_mul = " << this->prm_[GCal_mul] << "\n";
    oss << "CiPA_scale_IKr = " << this->prm_[CiPA_scale_IKr] << "\n";
    oss << "CiPA_scale_IKs = " << this->prm_[CiPA_scale_IKs] << "\n";
    oss << "CiPA_scale_IK1 = " << this->prm_[CiPA_scale_IK1] << "\n";
    oss << "CiPA_scale_ICaL = " << this->prm_[CiPA_scale_ICaL] << "\n";
    oss << "CiPA_scale_INaL = " << this->prm_[CiPA_scale_INaL] << "\n";
    oss << "nao = " << this->prm_[nao] << "\n";
    oss << "cao = " << this->prm_[cao] << "\n";
    oss << "ko = " << this->prm_[ko] << "\n";
    oss << "R = " << this->prm_[R] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "F = " << this->prm_[F] << "\n";
    oss << "L = " << this->prm_[L] << "\n";
    oss << "rad = " << this->prm_[rad] << "\n";
    oss << "vcell = " << this->prm_[vcell] << "\n";
    oss << "Ageo = " << this->prm_[Ageo] << "\n";
    oss << "Acap = " << this->prm_[Acap] << "\n";
    oss << "vmyo = " << this->prm_[vmyo] << "\n";
    oss << "vnsr = " << this->prm_[vnsr] << "\n";
    oss << "vjsr = " << this->prm_[vjsr] << "\n";
    oss << "vss = " << this->prm_[vss] << "\n";
    oss << "aCaMK = " << this->prm_[aCaMK] << "\n";
    oss << "bCaMK = " << this->prm_[bCaMK] << "\n";
    oss << "CaMKo = " << this->prm_[CaMKo] << "\n";
    oss << "PKNa = " << this->prm_[PKNa] << "\n";
    oss << "GNa = " << this->prm_[GNa] << "\n";
    oss << "GNaL = " << this->prm_[GNaL] << "\n";
    oss << "Gto = " << this->prm_[Gto] << "\n";
    oss << "GKr = " << this->prm_[GKr] << "\n";
    oss << "GKs = " << this->prm_[GKs] << "\n";
    oss << "GK1 = " << this->prm_[GK1] << "\n";
    oss << "Gncx = " << this->prm_[Gncx] << "\n";
    oss << "GKb = " << this->prm_[GKb] << "\n";
    oss << "GpCa = " << this->prm_[GpCa] << "\n";
    oss << "PCa = " << this->prm_[PCa] << "\n";
    oss << "Pnak = " << this->prm_[Pnak] << "\n";
    oss << "PNab = " << this->prm_[PNab] << "\n";
    oss << "PCab = " << this->prm_[PCab] << "\n";
    oss << "KmCaMK = " << this->prm_[KmCaMK] << "\n";
    oss << "KmCaM = " << this->prm_[KmCaM] << "\n";
    oss << "SERCA_total = " << this->prm_[SERCA_total] << "\n";
    oss << "RyR_total = " << this->prm_[RyR_total] << "\n";
    oss << "Leak_total = " << this->prm_[Leak_total] << "\n";
    oss << "Trans_total = " << this->prm_[Trans_total] << "\n";
    oss << "iso_L = " << this->prm_[iso_L] << "\n";
    oss << "ibmx = " << this->prm_[ibmx] << "\n";
    oss << "phys_T = " << this->prm_[phys_T] << "\n";
    oss << "phys_R = " << this->prm_[phys_R] << "\n";
    oss << "length = " << this->prm_[length] << "\n";
    oss << "cell_pi = " << this->prm_[cell_pi] << "\n";
    oss << "radius = " << this->prm_[radius] << "\n";
    oss << "geoArea = " << this->prm_[geoArea] << "\n";
    oss << "volume = " << this->prm_[volume] << "\n";
    oss << "v_cav = " << this->prm_[v_cav] << "\n";
    oss << "v_eca = " << this->prm_[v_eca] << "\n";
    oss << "v_cyt = " << this->prm_[v_cyt] << "\n";
    oss << "vr_cav = " << this->prm_[vr_cav] << "\n";
    oss << "vr_eca = " << this->prm_[vr_eca] << "\n";
    oss << "vr_cyt = " << this->prm_[vr_cyt] << "\n";
    oss << "j_cav_eca = " << this->prm_[j_cav_eca] << "\n";
    oss << "j_cav_cyt = " << this->prm_[j_cav_cyt] << "\n";
    oss << "j_eca_cyt = " << this->prm_[j_eca_cyt] << "\n";
    oss << "PKA_tot = " << this->prm_[PKA_tot] << "\n";
    oss << "f_cav = " << this->prm_[f_cav] << "\n";
    oss << "f_eca = " << this->prm_[f_eca] << "\n";
    oss << "f_cyt = " << this->prm_[f_cyt] << "\n";
    oss << "PKA_eca = " << this->prm_[PKA_eca] << "\n";
    oss << "PKA_cav = " << this->prm_[PKA_cav] << "\n";
    oss << "PKA_cyt = " << this->prm_[PKA_cyt] << "\n";
    oss << "PKI_tot = " << this->prm_[PKI_tot] << "\n";
    oss << "f_pki_cav = " << this->prm_[f_pki_cav] << "\n";
    oss << "f_pki_eca = " << this->prm_[f_pki_eca] << "\n";
    oss << "f_pki_cyt = " << this->prm_[f_pki_cyt] << "\n";
    oss << "PKI_cav = " << this->prm_[PKI_cav] << "\n";
    oss << "PKI_eca = " << this->prm_[PKI_eca] << "\n";
    oss << "PKI_cyt = " << this->prm_[PKI_cyt] << "\n";
    oss << "f_pki = " << this->prm_[f_pki] << "\n";
    oss << "K_pki = " << this->prm_[K_pki] << "\n";
    oss << "b_pki = " << this->prm_[b_pki] << "\n";
    oss << "pka_cav_f1 = " << this->prm_[pka_cav_f1] << "\n";
    oss << "pka_cav_f2 = " << this->prm_[pka_cav_f2] << "\n";
    oss << "pka_cav_f3 = " << this->prm_[pka_cav_f3] << "\n";
    oss << "pka_cav_K1 = " << this->prm_[pka_cav_K1] << "\n";
    oss << "pka_cav_K2 = " << this->prm_[pka_cav_K2] << "\n";
    oss << "pka_cav_K3 = " << this->prm_[pka_cav_K3] << "\n";
    oss << "pka_cav_b1 = " << this->prm_[pka_cav_b1] << "\n";
    oss << "pka_cav_b2 = " << this->prm_[pka_cav_b2] << "\n";
    oss << "pka_cav_b3 = " << this->prm_[pka_cav_b3] << "\n";
    oss << "pka_eca_f1 = " << this->prm_[pka_eca_f1] << "\n";
    oss << "pka_eca_f2 = " << this->prm_[pka_eca_f2] << "\n";
    oss << "pka_eca_f3 = " << this->prm_[pka_eca_f3] << "\n";
    oss << "pka_eca_K1 = " << this->prm_[pka_eca_K1] << "\n";
    oss << "pka_eca_K2 = " << this->prm_[pka_eca_K2] << "\n";
    oss << "pka_eca_K3 = " << this->prm_[pka_eca_K3] << "\n";
    oss << "pka_eca_b1 = " << this->prm_[pka_eca_b1] << "\n";
    oss << "pka_eca_b2 = " << this->prm_[pka_eca_b2] << "\n";
    oss << "pka_eca_b3 = " << this->prm_[pka_eca_b3] << "\n";
    oss << "pka_cyt_f1 = " << this->prm_[pka_cyt_f1] << "\n";
    oss << "pka_cyt_f2 = " << this->prm_[pka_cyt_f2] << "\n";
    oss << "pka_cyt_f3 = " << this->prm_[pka_cyt_f3] << "\n";
    oss << "pka_cyt_K1 = " << this->prm_[pka_cyt_K1] << "\n";
    oss << "pka_cyt_K2 = " << this->prm_[pka_cyt_K2] << "\n";
    oss << "pka_cyt_K3 = " << this->prm_[pka_cyt_K3] << "\n";
    oss << "pka_cyt_b1 = " << this->prm_[pka_cyt_b1] << "\n";
    oss << "pka_cyt_b2 = " << this->prm_[pka_cyt_b2] << "\n";
    oss << "pka_cyt_b3 = " << this->prm_[pka_cyt_b3] << "\n";
    oss << "f = " << this->prm_[f] << "\n";
    oss << "kp = " << this->prm_[kp] << "\n";
    oss << "kdp = " << this->prm_[kdp] << "\n";
    oss << "Kp = " << this->prm_[Kp] << "\n";
    oss << "Kdp = " << this->prm_[Kdp] << "\n";
    oss << "PP1_eca = " << this->prm_[PP1_eca] << "\n";
    oss << "PP1_cav = " << this->prm_[PP1_cav] << "\n";
    oss << "PP1_cyt = " << this->prm_[PP1_cyt] << "\n";
    oss << "pp1_K = " << this->prm_[pp1_K] << "\n";
    oss << "inhib1_tot = " << this->prm_[inhib1_tot] << "\n";
    oss << "PP2A = " << this->prm_[PP2A] << "\n";
    oss << "PDE2_tot = " << this->prm_[PDE2_tot] << "\n";
    oss << "f_pde2_cav = " << this->prm_[f_pde2_cav] << "\n";
    oss << "f_pde2_eca = " << this->prm_[f_pde2_eca] << "\n";
    oss << "f_pde2_cyt = " << this->prm_[f_pde2_cyt] << "\n";
    oss << "f_pde2_part = " << this->prm_[f_pde2_part] << "\n";
    oss << "f_pde_part = " << this->prm_[f_pde_part] << "\n";
    oss << "f_pde4_part = " << this->prm_[f_pde4_part] << "\n";
    oss << "f_pde4_cav = " << this->prm_[f_pde4_cav] << "\n";
    oss << "f_pde4_eca = " << this->prm_[f_pde4_eca] << "\n";
    oss << "f_pde4_cyt = " << this->prm_[f_pde4_cyt] << "\n";
    oss << "kPDE2 = " << this->prm_[kPDE2] << "\n";
    oss << "kPDE3 = " << this->prm_[kPDE3] << "\n";
    oss << "kPDE4 = " << this->prm_[kPDE4] << "\n";
    oss << "KmPDE2 = " << this->prm_[KmPDE2] << "\n";
    oss << "KmPDE3 = " << this->prm_[KmPDE3] << "\n";
    oss << "KmPDE4 = " << this->prm_[KmPDE4] << "\n";
    oss << "KmIbmxPde2 = " << this->prm_[KmIbmxPde2] << "\n";
    oss << "h_ibmx_pde2 = " << this->prm_[h_ibmx_pde2] << "\n";
    oss << "h_ibmx_pde3 = " << this->prm_[h_ibmx_pde3] << "\n";
    oss << "h_ibmx_pde4 = " << this->prm_[h_ibmx_pde4] << "\n";
    oss << "KmIbmxPde3 = " << this->prm_[KmIbmxPde3] << "\n";
    oss << "KmIbmxPde4 = " << this->prm_[KmIbmxPde4] << "\n";
    oss << "KPDEp = " << this->prm_[KPDEp] << "\n";
    oss << "delta_k_pde34 = " << this->prm_[delta_k_pde34] << "\n";
    oss << "ff_pde3_cyt = " << this->prm_[ff_pde3_cyt] << "\n";
    oss << "kfPDEp = " << this->prm_[kfPDEp] << "\n";
    oss << "r_pde34_frac = " << this->prm_[r_pde34_frac] << "\n";
    oss << "r_pde3_cyt = " << this->prm_[r_pde3_cyt] << "\n";
    oss << "kbPDEp = " << this->prm_[kbPDEp] << "\n";
    oss << "pde_PDE3_tot_alpha = " << this->prm_[pde_PDE3_tot_alpha] << "\n";
    oss << "pde_PDE3_tot_beta = " << this->prm_[pde_PDE3_tot_beta] << "\n";
    oss << "PDE3_tot = " << this->prm_[PDE3_tot] << "\n";
    oss << "PDE4_tot = " << this->prm_[PDE4_tot] << "\n";
    oss << "ibmx_h2 = " << this->prm_[ibmx_h2] << "\n";
    oss << "ibmx_h3 = " << this->prm_[ibmx_h3] << "\n";
    oss << "ibmx_h4 = " << this->prm_[ibmx_h4] << "\n";
    oss << "ibmx2 = " << this->prm_[ibmx2] << "\n";
    oss << "ibmx3 = " << this->prm_[ibmx3] << "\n";
    oss << "ibmx4 = " << this->prm_[ibmx4] << "\n";
    oss << "f_pde3_cav = " << this->prm_[f_pde3_cav] << "\n";
    oss << "f_pde3_cyt = " << this->prm_[f_pde3_cyt] << "\n";
    oss << "PDE2_cav = " << this->prm_[PDE2_cav] << "\n";
    oss << "PDE2_eca = " << this->prm_[PDE2_eca] << "\n";
    oss << "PDE2_cyt = " << this->prm_[PDE2_cyt] << "\n";
    oss << "PDE3_cav = " << this->prm_[PDE3_cav] << "\n";
    oss << "PDE3_cyt = " << this->prm_[PDE3_cyt] << "\n";
    oss << "PDE4_cav = " << this->prm_[PDE4_cav] << "\n";
    oss << "PDE4_eca = " << this->prm_[PDE4_eca] << "\n";
    oss << "PDE4_cyt = " << this->prm_[PDE4_cyt] << "\n";
    oss << "beta_R_b1_tot = " << this->prm_[beta_R_b1_tot] << "\n";
    oss << "beta_R_b2_tot = " << this->prm_[beta_R_b2_tot] << "\n";
    oss << "Gs_tot = " << this->prm_[Gs_tot] << "\n";
    oss << "Gi_tot = " << this->prm_[Gi_tot] << "\n";
    oss << "f_Gs_eca = " << this->prm_[f_Gs_eca] << "\n";
    oss << "f_Gs_cav = " << this->prm_[f_Gs_cav] << "\n";
    oss << "f_Gs_cyt = " << this->prm_[f_Gs_cyt] << "\n";
    oss << "f_Gi_cav = " << this->prm_[f_Gi_cav] << "\n";
    oss << "f_Gi_eca = " << this->prm_[f_Gi_eca] << "\n";
    oss << "f_Rb1_cav = " << this->prm_[f_Rb1_cav] << "\n";
    oss << "f_Rb1_eca = " << this->prm_[f_Rb1_eca] << "\n";
    oss << "f_Rb1_cyt = " << this->prm_[f_Rb1_cyt] << "\n";
    oss << "f_Rb2_cav = " << this->prm_[f_Rb2_cav] << "\n";
    oss << "f_Rb2_eca = " << this->prm_[f_Rb2_eca] << "\n";
    oss << "k_b1_l = " << this->prm_[k_b1_l] << "\n";
    oss << "k_b1_c = " << this->prm_[k_b1_c] << "\n";
    oss << "k_b1_h = " << this->prm_[k_b1_h] << "\n";
    oss << "k_b2_n = " << this->prm_[k_b2_n] << "\n";
    oss << "k_b2_h = " << this->prm_[k_b2_h] << "\n";
    oss << "k_b2_a = " << this->prm_[k_b2_a] << "\n";
    oss << "k_b2_c = " << this->prm_[k_b2_c] << "\n";
    oss << "k_b2_l = " << this->prm_[k_b2_l] << "\n";
    oss << "k_b2_f = " << this->prm_[k_b2_f] << "\n";
    oss << "k_act1_Gs = " << this->prm_[k_act1_Gs] << "\n";
    oss << "k_act2_Gs = " << this->prm_[k_act2_Gs] << "\n";
    oss << "k_act1_Gi = " << this->prm_[k_act1_Gi] << "\n";
    oss << "k_act2_Gi = " << this->prm_[k_act2_Gi] << "\n";
    oss << "k_hydr_Gs = " << this->prm_[k_hydr_Gs] << "\n";
    oss << "k_hydr_Gi = " << this->prm_[k_hydr_Gi] << "\n";
    oss << "k_reas_Gs = " << this->prm_[k_reas_Gs] << "\n";
    oss << "k_reas_Gi = " << this->prm_[k_reas_Gi] << "\n";
    oss << "rate_bds = " << this->prm_[rate_bds] << "\n";
    oss << "k_grk_dp = " << this->prm_[k_grk_dp] << "\n";
    oss << "k_grk_p = " << this->prm_[k_grk_p] << "\n";
    oss << "k_pka_p = " << this->prm_[k_pka_p] << "\n";
    oss << "k_pka_dp = " << this->prm_[k_pka_dp] << "\n";
    oss << "beta_cav_GRK = " << this->prm_[beta_cav_GRK] << "\n";
    oss << "beta_cav_R_b1_tot = " << this->prm_[beta_cav_R_b1_tot] << "\n";
    oss << "beta_cav_k_GsAct_b2 = " << this->prm_[beta_cav_k_GsAct_b2] << "\n";
    oss << "beta_cav_R_b2_tot = " << this->prm_[beta_cav_R_b2_tot] << "\n";
    oss << "beta_cav_Rb2_pka_f_a = " << this->prm_[beta_cav_Rb2_pka_f_a] << "\n";
    oss << "beta_cav_Gs_f_c22 = " << this->prm_[beta_cav_Gs_f_c22] << "\n";
    oss << "beta_cav_Gs_f_a = " << this->prm_[beta_cav_Gs_f_a] << "\n";
    oss << "beta_cav_Gs_f_c33 = " << this->prm_[beta_cav_Gs_f_c33] << "\n";
    oss << "beta_cav_Gs_f_c11 = " << this->prm_[beta_cav_Gs_f_c11] << "\n";
    oss << "beta_eca_GRK = " << this->prm_[beta_eca_GRK] << "\n";
    oss << "beta_eca_k_GsAct_b2 = " << this->prm_[beta_eca_k_GsAct_b2] << "\n";
    oss << "beta_eca_R_b2_tot = " << this->prm_[beta_eca_R_b2_tot] << "\n";
    oss << "beta_eca_R_b1_tot = " << this->prm_[beta_eca_R_b1_tot] << "\n";
    oss << "beta_eca_Rb2_pka_f_a = " << this->prm_[beta_eca_Rb2_pka_f_a] << "\n";
    oss << "beta_eca_Gs_f_c11 = " << this->prm_[beta_eca_Gs_f_c11] << "\n";
    oss << "beta_eca_Gs_f_c33 = " << this->prm_[beta_eca_Gs_f_c33] << "\n";
    oss << "beta_eca_Gs_f_a = " << this->prm_[beta_eca_Gs_f_a] << "\n";
    oss << "beta_eca_Gs_f_c22 = " << this->prm_[beta_eca_Gs_f_c22] << "\n";
    oss << "beta_cyt_R_b1_tot = " << this->prm_[beta_cyt_R_b1_tot] << "\n";
    oss << "beta_cyt_GRK = " << this->prm_[beta_cyt_GRK] << "\n";
    oss << "beta_cyt_Rb1_np_f_a = " << this->prm_[beta_cyt_Rb1_np_f_a] << "\n";
    oss << "ATP = " << this->prm_[ATP] << "\n";
    oss << "KmATP = " << this->prm_[KmATP] << "\n";
    oss << "AC_tot = " << this->prm_[AC_tot] << "\n";
    oss << "f_AC47_eca = " << this->prm_[f_AC47_eca] << "\n";
    oss << "f_AC56_cav = " << this->prm_[f_AC56_cav] << "\n";
    oss << "f_AC56_AC47 = " << this->prm_[f_AC56_AC47] << "\n";
    oss << "hGsAC47 = " << this->prm_[hGsAC47] << "\n";
    oss << "hGsAC56 = " << this->prm_[hGsAC56] << "\n";
    oss << "hGsGiAC56 = " << this->prm_[hGsGiAC56] << "\n";
    oss << "KmGsAC47 = " << this->prm_[KmGsAC47] << "\n";
    oss << "KmGsAC56 = " << this->prm_[KmGsAC56] << "\n";
    oss << "KmGiAC56 = " << this->prm_[KmGiAC56] << "\n";
    oss << "KmGsGiAC56 = " << this->prm_[KmGsGiAC56] << "\n";
    oss << "basalAC47 = " << this->prm_[basalAC47] << "\n";
    oss << "basalAC56 = " << this->prm_[basalAC56] << "\n";
    oss << "afAC47 = " << this->prm_[afAC47] << "\n";
    oss << "afAC56 = " << this->prm_[afAC56] << "\n";
    oss << "vGsGiAC56 = " << this->prm_[vGsGiAC56] << "\n";
    oss << "AC47_cyt = " << this->prm_[AC47_cyt] << "\n";
    oss << "AC56_cav = " << this->prm_[AC56_cav] << "\n";
    oss << "fATP = " << this->prm_[fATP] << "\n";
    oss << "AC47_eca = " << this->prm_[AC47_eca] << "\n";
    oss << "AC56_cyt = " << this->prm_[AC56_cyt] << "\n";
    oss << "ka_inak = " << this->prm_[ka_inak] << "\n";
    oss << "kp_inak = " << this->prm_[kp_inak] << "\n";
    oss << "Ka_inak = " << this->prm_[Ka_inak] << "\n";
    oss << "Kp_inak = " << this->prm_[Kp_inak] << "\n";
    oss << "Ka_ina = " << this->prm_[Ka_ina] << "\n";
    oss << "Kp_ina = " << this->prm_[Kp_ina] << "\n";
    oss << "ka_ina = " << this->prm_[ka_ina] << "\n";
    oss << "kp_ina = " << this->prm_[kp_ina] << "\n";
    oss << "ka_plb = " << this->prm_[ka_plb] << "\n";
    oss << "kp_plb = " << this->prm_[kp_plb] << "\n";
    oss << "Ka_plb = " << this->prm_[Ka_plb] << "\n";
    oss << "Kp_plb = " << this->prm_[Kp_plb] << "\n";
    oss << "ka_ikur = " << this->prm_[ka_ikur] << "\n";
    oss << "kp_ikur = " << this->prm_[kp_ikur] << "\n";
    oss << "Ka_ikur = " << this->prm_[Ka_ikur] << "\n";
    oss << "Kp_ikur = " << this->prm_[Kp_ikur] << "\n";
    oss << "ka_iks = " << this->prm_[ka_iks] << "\n";
    oss << "kp_iks = " << this->prm_[kp_iks] << "\n";
    oss << "Ka_iks = " << this->prm_[Ka_iks] << "\n";
    oss << "Kp_iks = " << this->prm_[Kp_iks] << "\n";
    oss << "M = " << this->prm_[M] << "\n";
    oss << "iks_sig_L = " << this->prm_[iks_sig_L] << "\n";
    oss << "iks_sig_K = " << this->prm_[iks_sig_K] << "\n";
    oss << "Yotiao = " << this->prm_[Yotiao] << "\n";
    oss << "IKs_tot = " << this->prm_[IKs_tot] << "\n";
    oss << "iks_sig_PKAf_sum = " << this->prm_[iks_sig_PKAf_sum] << "\n";
    oss << "iks_sig_PKAf = " << this->prm_[iks_sig_PKAf] << "\n";
    oss << "iks_sig_PP1f_eca_sum = " << this->prm_[iks_sig_PP1f_eca_sum] << "\n";
    oss << "PP1f_eca = " << this->prm_[PP1f_eca] << "\n";
    oss << "iks_sig_IKsf_sum = " << this->prm_[iks_sig_IKsf_sum] << "\n";
    oss << "IKsf = " << this->prm_[IKsf] << "\n";
    oss << "Yotiaof = " << this->prm_[Yotiaof] << "\n";
    oss << "IKs_arn = " << this->prm_[IKs_arn] << "\n";
    oss << "IKs_arp = " << this->prm_[IKs_arp] << "\n";
    oss << "ka_tni = " << this->prm_[ka_tni] << "\n";
    oss << "kp_tni = " << this->prm_[kp_tni] << "\n";
    oss << "Ka_tni = " << this->prm_[Ka_tni] << "\n";
    oss << "Kp_tni = " << this->prm_[Kp_tni] << "\n";
    oss << "RyR_tot = " << this->prm_[RyR_tot] << "\n";
    oss << "RyR_akap = " << this->prm_[RyR_akap] << "\n";
    oss << "ka_ryr = " << this->prm_[ka_ryr] << "\n";
    oss << "kp_ryr = " << this->prm_[kp_ryr] << "\n";
    oss << "Ka_ryr = " << this->prm_[Ka_ryr] << "\n";
    oss << "Kp_ryr = " << this->prm_[Kp_ryr] << "\n";
    oss << "Mr = " << this->prm_[Mr] << "\n";
    oss << "Lr = " << this->prm_[Lr] << "\n";
    oss << "Kr = " << this->prm_[Kr] << "\n";
    oss << "akap_sig_RyRf_sum = " << this->prm_[akap_sig_RyRf_sum] << "\n";
    oss << "RyRf = " << this->prm_[RyRf] << "\n";
    oss << "ICaL_tot = " << this->prm_[ICaL_tot] << "\n";
    oss << "ICaL_akap = " << this->prm_[ICaL_akap] << "\n";
    oss << "ka_ical = " << this->prm_[ka_ical] << "\n";
    oss << "kp_ical = " << this->prm_[kp_ical] << "\n";
    oss << "Ka_ical = " << this->prm_[Ka_ical] << "\n";
    oss << "Kp_ical = " << this->prm_[Kp_ical] << "\n";
    oss << "Mi = " << this->prm_[Mi] << "\n";
    oss << "Li = " << this->prm_[Li] << "\n";
    oss << "Ki = " << this->prm_[Ki] << "\n";
    oss << "akap_sig_ICaLf_sum = " << this->prm_[akap_sig_ICaLf_sum] << "\n";
    oss << "ICaLf = " << this->prm_[ICaLf] << "\n";
    oss << "akap_sig_PP1f_cav_b = " << this->prm_[akap_sig_PP1f_cav_b] << "\n";
    oss << "akap_sig_PP1f_cav_c = " << this->prm_[akap_sig_PP1f_cav_c] << "\n";
    oss << "akap_sig_PP1f_cav_d = " << this->prm_[akap_sig_PP1f_cav_d] << "\n";
    oss << "akap_sig_PP1f_cav_rr = " << this->prm_[akap_sig_PP1f_cav_rr] << "\n";
    oss << "akap_sig_PP1f_cav_yi = " << this->prm_[akap_sig_PP1f_cav_yi] << "\n";
    oss << "akap_sig_PP1f_cav_yr = " << this->prm_[akap_sig_PP1f_cav_yr] << "\n";
    oss << "akap_sig_PP1f_cav_mag = " << this->prm_[akap_sig_PP1f_cav_mag] << "\n";
    oss << "akap_sig_PP1f_cav_arg = " << this->prm_[akap_sig_PP1f_cav_arg] << "\n";
    oss << "akap_sig_PP1f_cav_x = " << this->prm_[akap_sig_PP1f_cav_x] << "\n";
    oss << "PP1f_cav = " << this->prm_[PP1f_cav] << "\n";
    oss << "akap_sig_PKAf_d = " << this->prm_[akap_sig_PKAf_d] << "\n";
    oss << "akap_sig_PKAf_b = " << this->prm_[akap_sig_PKAf_b] << "\n";
    oss << "akap_sig_PKAf_c = " << this->prm_[akap_sig_PKAf_c] << "\n";
    oss << "akap_sig_PKAf_rr = " << this->prm_[akap_sig_PKAf_rr] << "\n";
    oss << "akap_sig_PKAf_yr = " << this->prm_[akap_sig_PKAf_yr] << "\n";
    oss << "akap_sig_PKAf_yi = " << this->prm_[akap_sig_PKAf_yi] << "\n";
    oss << "akap_sig_PKAf_mag = " << this->prm_[akap_sig_PKAf_mag] << "\n";
    oss << "akap_sig_PKAf_arg = " << this->prm_[akap_sig_PKAf_arg] << "\n";
    oss << "akap_sig_PKAf_x = " << this->prm_[akap_sig_PKAf_x] << "\n";
    oss << "akap_sig_PKAf = " << this->prm_[akap_sig_PKAf] << "\n";
    oss << "RyR_akapf = " << this->prm_[RyR_akapf] << "\n";
    oss << "RyR_arn = " << this->prm_[RyR_arn] << "\n";
    oss << "RyR_arp = " << this->prm_[RyR_arp] << "\n";
    oss << "ICaL_akapf = " << this->prm_[ICaL_akapf] << "\n";
    oss << "ICaL_arn = " << this->prm_[ICaL_arn] << "\n";
    oss << "ICaL_arp = " << this->prm_[ICaL_arp] << "\n";
    oss << "beta_0 = " << this->prm_[beta_0] << "\n";
    oss << "irel_fhat_ratio = " << this->prm_[irel_fhat_ratio] << "\n";
    oss << "ical_f_hat_ratio = " << this->prm_[ical_f_hat_ratio] << "\n";
    oss << "iks_f_hat_ratio = " << this->prm_[iks_f_hat_ratio] << "\n";
    oss << "Whole_cell_PP1 = " << this->prm_[Whole_cell_PP1] << "\n";
    oss << "PP1block = " << this->prm_[PP1block] << "\n";
    oss << "fICaLP = " << this->prm_[fICaLP] << "\n";
    oss << "fIKsP = " << this->prm_[fIKsP] << "\n";
    oss << "fPLBP = " << this->prm_[fPLBP] << "\n";
    oss << "fTnIP = " << this->prm_[fTnIP] << "\n";
    oss << "fINaP = " << this->prm_[fINaP] << "\n";
    oss << "fINaKP = " << this->prm_[fINaKP] << "\n";
    oss << "fRyRP = " << this->prm_[fRyRP] << "\n";
    oss << "fIKurP = " << this->prm_[fIKurP];
    return oss.str();

}


std::string Gong2020::PrintCurrents() const
{
    using namespace Gng20Cur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "INaL = " << this->cur_[INaL] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "ICaNa = " << this->cur_[ICaNa] << "\n";
    oss << "ICaK = " << this->cur_[ICaK] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "INaCa_i = " << this->cur_[INaCa_i] << "\n";
    oss << "INaCa_ss = " << this->cur_[INaCa_ss] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "IKb = " << this->cur_[IKb] << "\n";
    oss << "INab = " << this->cur_[INab] << "\n";
    oss << "ICab = " << this->cur_[ICab] << "\n";
    oss << "IpCa = " << this->cur_[IpCa] << "\n";
    oss << "Iion = " << this->cur_[Gng20Cur::Iion];
    return oss.str();

}


std::string Gong2020::PrintBlockCoeffs() const
{
    using namespace Gng20Cur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->block_coeff_[INa] << "\n";
    oss << "INaL = " << this->block_coeff_[INaL] << "\n";
    oss << "Ito = " << this->block_coeff_[Ito] << "\n";
    oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
    oss << "ICaNa = " << this->block_coeff_[ICaNa] << "\n";
    oss << "ICaK = " << this->block_coeff_[ICaK] << "\n";
    oss << "IKr = " << this->block_coeff_[IKr] << "\n";
    oss << "IKs = " << this->block_coeff_[IKs] << "\n";
    oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
    oss << "INaCa_ss = " << this->block_coeff_[INaCa_i] << "\n";
    oss << "INaCa_ss = " << this->block_coeff_[INaCa_ss] << "\n";
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "IKb = " << this->block_coeff_[IKb] << "\n";
    oss << "INab = " << this->block_coeff_[INab] << "\n";
    oss << "ICab = " << this->block_coeff_[ICab] << "\n";
    oss << "IpCa = " << this->block_coeff_[IpCa];
    return oss.str();
}

} // End of namespace ELECTRA
