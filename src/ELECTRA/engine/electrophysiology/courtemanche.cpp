/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/electrophysiology/courtemanche.hpp"


namespace ELECTRA {


void Courtemanche::SetDataMapping() 
{
    using namespace CourteVar;
    using namespace CourtePrm;
    using namespace CourteCur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["g_u"] = static_cast<std::size_t>(g_u);
    this->mapped_data_["g_v"] = static_cast<std::size_t>(g_v);
    this->mapped_data_["g_w"] = static_cast<std::size_t>(g_w);
    this->mapped_data_["d"] = static_cast<std::size_t>(d);
    this->mapped_data_["f_Ca"] = static_cast<std::size_t>(f_Ca);
    this->mapped_data_["f"] = static_cast<std::size_t>(f);
    this->mapped_data_["h"] = static_cast<std::size_t>(h);
    this->mapped_data_["j"] = static_cast<std::size_t>(j);
    this->mapped_data_["m"] = static_cast<std::size_t>(m);
    this->mapped_data_["Ca_i"] = static_cast<std::size_t>(Ca_i);
    this->mapped_data_["Ca_rel"] = static_cast<std::size_t>(Ca_rel);
    this->mapped_data_["Ca_up"] = static_cast<std::size_t>(Ca_up);
    this->mapped_data_["K_i"] = static_cast<std::size_t>(K_i);
    this->mapped_data_["Na_i"] = static_cast<std::size_t>(Na_i);
    this->mapped_data_["xr"] = static_cast<std::size_t>(xr);
    this->mapped_data_["xs"] = static_cast<std::size_t>(xs);
    this->mapped_data_["oa"] = static_cast<std::size_t>(oa);
    this->mapped_data_["oi"] = static_cast<std::size_t>(oi);
    this->mapped_data_["ua"] = static_cast<std::size_t>(ua);
    this->mapped_data_["ui"] = static_cast<std::size_t>(ui);
    this->mapped_data_["y"] = static_cast<std::size_t>(y);

    // Set parameters mapping.
    this->mapped_data_["rkr"] = static_cast<std::size_t>(rkr);
    this->mapped_data_["TYPE1"] = static_cast<std::size_t>(TYPE1);
    this->mapped_data_["TYPE3"] = static_cast<std::size_t>(TYPE3);
    this->mapped_data_["ISOdlf"] = static_cast<std::size_t>(ISOdlf);
    this->mapped_data_["cAF"] = static_cast<std::size_t>(cAF);
    this->mapped_data_["cAF2"] = static_cast<std::size_t>(cAF2);
    this->mapped_data_["SR"] = static_cast<std::size_t>(SR);
    this->mapped_data_["Ach"] = static_cast<std::size_t>(Ach);
    this->mapped_data_["drugvector"] = static_cast<std::size_t>(drugvector);
    this->mapped_data_["Engel"] = static_cast<std::size_t>(Engel);
    this->mapped_data_["Penaranda"] = static_cast<std::size_t>(Penaranda);
    this->mapped_data_["Cm"] = static_cast<std::size_t>(Cm);
    this->mapped_data_["CMDN_max"] = static_cast<std::size_t>(CMDN_max);
    this->mapped_data_["CSQN_max"] = static_cast<std::size_t>(CSQN_max);
    this->mapped_data_["Km_CMDN"] = static_cast<std::size_t>(Km_CMDN);
    this->mapped_data_["Km_CSQN"] = static_cast<std::size_t>(Km_CSQN);
    this->mapped_data_["Km_TRPN"] = static_cast<std::size_t>(Km_TRPN);
    this->mapped_data_["TRPN_max"] = static_cast<std::size_t>(TRPN_max);
    this->mapped_data_["Ca_up_max"] = static_cast<std::size_t>(Ca_up_max);
    this->mapped_data_["K_rel"] = static_cast<std::size_t>(K_rel);
    this->mapped_data_["I_up_max"] = static_cast<std::size_t>(I_up_max);
    this->mapped_data_["K_up"] = static_cast<std::size_t>(K_up);
    this->mapped_data_["I_NaCa_max"] = static_cast<std::size_t>(I_NaCa_max);
    this->mapped_data_["K_mCa"] = static_cast<std::size_t>(K_mCa);
    this->mapped_data_["K_mNa"] = static_cast<std::size_t>(K_mNa);
    this->mapped_data_["K_sat"] = static_cast<std::size_t>(K_sat);
    this->mapped_data_["gamma"] = static_cast<std::size_t>(gamma);
    this->mapped_data_["g_B_Ca"] = static_cast<std::size_t>(g_B_Ca);
    this->mapped_data_["g_B_K"] = static_cast<std::size_t>(g_B_K);
    this->mapped_data_["g_B_Na"] = static_cast<std::size_t>(g_B_Na);
    this->mapped_data_["g_Na"] = static_cast<std::size_t>(g_Na);
    this->mapped_data_["g_Kca_eb"] = static_cast<std::size_t>(g_Kca_eb);
    this->mapped_data_["g_Kca_pb"] = static_cast<std::size_t>(g_Kca_pb);
    this->mapped_data_["V_cell"] = static_cast<std::size_t>(V_cell);
    this->mapped_data_["V_i"] = static_cast<std::size_t>(V_i);
    this->mapped_data_["V_up"] = static_cast<std::size_t>(V_up);
    this->mapped_data_["V_rel"] = static_cast<std::size_t>(V_rel);
    this->mapped_data_["F"] = static_cast<std::size_t>(F);
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["i_CaP_max"] = static_cast<std::size_t>(i_CaP_max);
    this->mapped_data_["Km_K_o"] = static_cast<std::size_t>(Km_K_o);
    this->mapped_data_["Km_Na_i"] = static_cast<std::size_t>(Km_Na_i);
    this->mapped_data_["i_NaK_max"] = static_cast<std::size_t>(i_NaK_max);
    this->mapped_data_["Ca_o"] = static_cast<std::size_t>(Ca_o);
    this->mapped_data_["K_o"] = static_cast<std::size_t>(K_o);
    this->mapped_data_["Na_o"] = static_cast<std::size_t>(Na_o);
    this->mapped_data_["tau_tr"] = static_cast<std::size_t>(tau_tr);
    this->mapped_data_["K_Q10"] = static_cast<std::size_t>(K_Q10);
    this->mapped_data_["EAD"] = static_cast<std::size_t>(EAD);

    // Set currents mapping.
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);
    this->mapped_data_["Ito"] = static_cast<std::size_t>(Ito);
    this->mapped_data_["IKur"] = static_cast<std::size_t>(IKur);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["ICaL"] = static_cast<std::size_t>(ICaL);
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["IKach"] = static_cast<std::size_t>(IKach);
    this->mapped_data_["IKca_E"] = static_cast<std::size_t>(IKca_E);
    this->mapped_data_["IKca_P"] = static_cast<std::size_t>(IKca_P);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);
    this->mapped_data_["ICap"] = static_cast<std::size_t>(ICap);
    this->mapped_data_["IRel"] = static_cast<std::size_t>(IRel);
    this->mapped_data_["Itr"] = static_cast<std::size_t>(Itr);
    this->mapped_data_["IupLeak"] = static_cast<std::size_t>(IupLeak);
    this->mapped_data_["Iup"] = static_cast<std::size_t>(Iup);

}


Courtemanche::Courtemanche()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Courtemanche;
    this->dt_stable_ = 0.02;
    this->var_.resize(23, 0.);
    this->prm_.resize(50, 0.);
    this->cur_.resize(20, 0.);
    this->block_coeff_.resize(19, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


Courtemanche::~Courtemanche()
{}


void Courtemanche::Initialize(CellType cell_type)
{
    using namespace CourteVar;
    using namespace CourtePrm;

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(23, 0.);
    this->prm_.clear();           this->prm_.resize(50, 0.);
    this->cur_.clear();           this->cur_.resize(20, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(19, 0.);

    // Set the model variables.
    this->var_[v] = -86.2149239355137;
    this->var_[dvdt] = 0.;
    this->var_[g_u] = -6.3166405898509e-21;
    this->var_[g_v] = 1.;
    this->var_[g_w] = 0.999403848287359;
    this->var_[d] = 7.28606563749481e-05;
    this->var_[f_Ca] = 0.805807026738844;
    this->var_[f] = 0.973413382551336;
    this->var_[h] = 0.988486881470332;
    this->var_[j] = 0.992618996204717;
    this->var_[m] = 0.00126181348739113;
    this->var_[Ca_i] = 8.43005909327977e-05;
    this->var_[Ca_rel] = 0.663346896943366;
    this->var_[Ca_up] = 1.05131464949419;
    this->var_[K_i] = 139.439310636634;
    this->var_[Na_i] = 12.158477246837;
    this->var_[xr] = 0.000240649802414847;
    this->var_[xs] = 0.0154077609021343;
    this->var_[oa] = 0.0230157860030608;
    this->var_[oi] = 0.999706921551885;
    this->var_[ua] = 0.00286417924596937;
    this->var_[ui] = 0.99195261769904;
    this->var_[y] = 0.00113765125606322;

    // Set model parameters according to the cell type.
    switch (cell_type)
    {
        case CellType::left_atrial :
            this->prm_[rkr] = 1.6;
            break;

        case CellType::right_atrial :
            this->prm_[rkr] = 0.;
            break;
    
        default:
            std::string error_str = "Could not initialize the Courtemanche atrial ap model. Expected cell type: ELECTRA::CellType::left_atrial or ELECTRA::CellType::right_atrial";
            throw std::invalid_argument(Logger::Error(error_str));
            break;
    }

    // Set the univarsal model parameters.
    this->prm_[TYPE1] = 0.;
    this->prm_[TYPE3] = 0.;
    this->prm_[ISOdlf] = 0.;
    this->prm_[cAF] = 0.;
    this->prm_[cAF2] = 0.;
    this->prm_[SR] = 1.;
    this->prm_[Ach] = 0.1;
    this->prm_[drugvector] = 0.;
    this->prm_[Engel] = 1.;
    this->prm_[Penaranda] = 0.;
    this->prm_[Cm] = 100.;
    this->prm_[CMDN_max] = 0.05;
    this->prm_[CSQN_max] = 10.;
    this->prm_[Km_CMDN] = 0.00238;
    this->prm_[Km_CSQN] = 0.8;
    this->prm_[Km_TRPN] = 0.0005;
    this->prm_[TRPN_max] = 0.07;
    this->prm_[Ca_up_max] = 15.;
    this->prm_[K_rel] = 30.;
    this->prm_[I_up_max] = 0.005;
    this->prm_[K_up] = 0.00092;
    this->prm_[I_NaCa_max] = 1600.;
    this->prm_[K_mCa] = 1.38;
    this->prm_[K_mNa] = 87.5;
    this->prm_[K_sat] = 0.1;
    this->prm_[gamma] = 0.35;
    this->prm_[g_B_Ca] = 0.001131;
    this->prm_[g_B_K] = 0.;
    this->prm_[g_B_Na] = 0.0006744375;
    this->prm_[g_Na] = 7.8;
    this->prm_[g_Kca_eb] = 2.8;
    this->prm_[g_Kca_pb] = 0.01;
    this->prm_[V_cell] = 20100.;
    this->prm_[V_i] = 13668.;
    this->prm_[V_up] = 1109.52;
    this->prm_[V_rel] = 96.48;
    this->prm_[F] = 96.4867;
    this->prm_[R] = 8.3143;
    this->prm_[T] = 310.;
    this->prm_[i_CaP_max] = 0.275;
    this->prm_[Km_K_o] = 1.5;
    this->prm_[Km_Na_i] = 10.;
    this->prm_[i_NaK_max] = 0.59933874;
    this->prm_[Ca_o] = 1.8;
    this->prm_[K_o] = 5.4;
    this->prm_[Na_o] = 140.;
    this->prm_[tau_tr] = 180.;
    this->prm_[K_Q10] = 3.;
    this->prm_[EAD] = 0.;

}


void Courtemanche::Compute(double v_new, double dt, double stim_current)
{
    using namespace CourteVar;
    using namespace CourtePrm;
    using namespace CourteCur;

    // ISO fitting curves --> given the parameter prm.ISOdlf they give the
    // factor by which increase the ICaL, Ito and IKur currents
    double fit_ICaL_SR_001 = this->prm_[SR]*(58.*std::atan(0.12*(0.01-13.))+70.);
    double fit_Ito_SR_001 = this->prm_[SR]*(40.533*std::atan(0.676*(0.01+0.653))-15.075);
    double fit_IKs_SR_001 = this->prm_[SR]*(22341.478*std::atan(497.52*(0.01+0.679))-35028.19);
    double fit_ICaL_cAF_001 = this->prm_[cAF]*(102.23*std::atan(0.733*(0.01-4.614))+141.65);
    double fit_Ito_cAF_001 = this->prm_[cAF]*(45.922*std::atan(1.327*(0.01+0.196))-10.667);
    double fit_IKs_cAF_001 = this->prm_[cAF]*(68.55*std::atan(5.168*(0.01+0.075))-28.46);

    double fit_ICaL_SR = this->prm_[SR]*(58.*std::atan(0.12*(this->prm_[ISOdlf]-13.))+70.);
    double fit_Ito_SR = this->prm_[SR]*(40.533*std::atan(0.676*(this->prm_[ISOdlf]+0.653))-15.075);
    double fit_IKs_SR = this->prm_[SR]*(22341.478*std::atan(497.52*(this->prm_[ISOdlf]+0.679))-35028.19);
    double fit_ICaL_cAF = this->prm_[cAF]*(102.23*std::atan(0.733*(this->prm_[ISOdlf]-4.614))+141.65);
    double fit_Ito_cAF = this->prm_[cAF]*(45.922*std::atan(1.327*(this->prm_[ISOdlf]+0.196))-10.667);
    double fit_IKs_cAF = this->prm_[cAF]*(68.55*std::atan(5.168*(this->prm_[ISOdlf]+0.075))-28.46);
    if (this->prm_[ISOdlf] < 0.01) {
        fit_ICaL_SR = this->prm_[SR]*(this->prm_[ISOdlf]*fit_ICaL_SR_001*100.);
        fit_Ito_SR = this->prm_[SR]*(this->prm_[ISOdlf]*fit_Ito_SR_001*100.);
        fit_IKs_SR = this->prm_[SR]*(this->prm_[ISOdlf]*fit_IKs_SR_001*100.);
        fit_ICaL_cAF = this->prm_[cAF]*(this->prm_[ISOdlf]*fit_ICaL_cAF_001*100.);
        fit_Ito_cAF = this->prm_[cAF]*(this->prm_[ISOdlf]*fit_Ito_cAF_001*100.);
        fit_IKs_cAF = this->prm_[cAF]*(this->prm_[ISOdlf]*fit_IKs_cAF_001*100.);
    }

    double g_Kr = this->prm_[rkr]*0.029411765;
    double g_Ks = 0.12941176*(1.+fit_IKs_cAF/100.)*(1.+fit_IKs_SR/100.);
    double g_to = 0.1652*(1.0-0.7*this->prm_[TYPE1])*(1.0-fit_Ito_cAF/100.)*(1.0-fit_Ito_SR/100.)*(1.0-0.50*this->prm_[cAF]-0.15*this->prm_[cAF2]);
    double g_Ca_L = 0.12375*(1.0-0.7*this->prm_[TYPE3])*(1.0+fit_ICaL_cAF/100.)*(1.0+fit_ICaL_SR/100.)*(1.0-0.70*this->prm_[cAF]+0.05*this->prm_[cAF2]);
    double g_K1 = 0.09*(1.0+1.1*this->prm_[cAF2]);

    // Equilibrium potentials
    double E_Na = this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log(this->prm_[Na_o]/this->var_[Na_i]);
    double E_Ca = this->prm_[R]*this->prm_[T]/(2.*this->prm_[F])*std::log(this->prm_[Ca_o]/this->var_[Ca_i]);
    double E_K = this->prm_[R]*this->prm_[T]/this->prm_[F]*std::log(this->prm_[K_o]/this->var_[K_i]);

    // Fast Na+ Current
    this->cur_[INa] = (1.0-this->block_coeff_[INa]) * (this->prm_[Cm] * this->prm_[g_Na] * this->var_[m]*this->var_[m]*this->var_[m] * this->var_[h] * this->var_[j] * (v_new-E_Na));
    
    // m
    double alpha_m = 0.;
    if (v_new == -47.13) { alpha_m = 3.2; }
    else { alpha_m = 0.32*(v_new+47.13)/(1.-std::exp(-0.1*(v_new+47.13))); }
    double beta_m = 0.08*std::exp(-v_new/11.);
    double m_inf = alpha_m/(alpha_m+beta_m);
    double tau_m = std::pow(alpha_m+beta_m, -1);
    this->var_[m] = ALGORITHM::RushLarsen(m_inf, this->var_[m], dt, tau_m);

    // h
    double alpha_h = 0.;
    double beta_h = 0.;
    if (v_new < -40.) { 
        alpha_h = 0.135*std::exp((v_new+80.)/(-6.8)); 
        beta_h = 3.56*std::exp(0.079*v_new)+3.1e5*std::exp(0.35*v_new); 
    }
    else { 
        beta_h = 1./(0.13*(1.+std::exp((v_new+10.66)/-11.1))); 
    }
    double h_inf = alpha_h/(alpha_h+beta_h);
    double tau_h = 1./(alpha_h+beta_h);
    this->var_[h] = ALGORITHM::RushLarsen(h_inf, this->var_[h], dt, tau_h);
    
    // j 
    double alpha_j = 0.;
    double beta_j = 0.;
    if (v_new < -40.) { 
        alpha_j = (-1.2714e5*std::exp(0.2444*v_new)-3.474e-5*std::exp(-0.04391*v_new))*(v_new+37.78)/(1.+std::exp(0.311*(v_new+79.23)));
        beta_j = 0.1212*std::exp(-0.01052*v_new)/(1.+std::exp(-0.1378*(v_new+40.14)));
    }
    else { 
        beta_j = 0.3*std::exp(-2.535e-7*v_new)/(1.+std::exp(-0.1*(v_new+32.)));
    }
    double j_inf = alpha_j/(alpha_j+beta_j);
    double tau_j = 1. / (alpha_j+beta_j);
    this->var_[j] = ALGORITHM::RushLarsen(j_inf, this->var_[j], dt, tau_j);

    // Time-Independent K+ Current
    this->cur_[IK1] =  (1.0-this->block_coeff_[IK1]) * (this->prm_[Cm]*this->prm_[g_K1]*(v_new-E_K)/(1.+std::exp(0.07*(v_new+80.))));

    // Transient Outward K+ Current
    this->cur_[Ito] =  (1.0-this->block_coeff_[Ito]) * (this->prm_[Cm] * g_to * this->var_[oa]*this->var_[oa]*this->var_[oa] * this->var_[oi] * (v_new-E_K));

    // oa
    double alpha_oa = 0.65*std::pow((std::exp((v_new+10.)/-8.5)+std::exp((v_new-30.)/-59.)), -1.);
    double beta_oa = 0.65*std::pow((2.5+std::exp((v_new+82.)/17.)), -1.);
    double tau_oa = std::pow((alpha_oa+beta_oa), -1.)/this->prm_[K_Q10];
    double oa_inf = std::pow((1.+std::exp((v_new+20.47)/-17.54)), -1.);
    this->var_[oa] = ALGORITHM::RushLarsen(oa_inf, this->var_[oa], dt, tau_oa);
    
    // oi
    double alpha_oi = std::pow((18.53+1.0*std::exp((v_new+113.7)/10.95)), -1.);
    double beta_oi = std::pow((35.56+1.0*std::exp((v_new+1.26)/-7.44)), -1.);
    double tau_oi = std::pow((alpha_oi+beta_oi), -1.)/this->prm_[K_Q10];
    double oi_inf = std::pow((1.+std::exp((v_new+43.1)/5.3)), -1.);
    this->var_[oi] = ALGORITHM::RushLarsen(oi_inf, this->var_[oi], dt, tau_oi);

    // Ultra rapid Delayed Rectifier K+ Current
    double g_Kur = (0.005+0.05/(1.+std::exp((v_new-15.)/-13.)))*(1.0 - 0.5*this->prm_[cAF]+0.01*this->prm_[cAF2]);
    this->cur_[IKur] =  (1.0-this->block_coeff_[IKur]) * (this->prm_[Cm]* g_Kur * this->var_[ua]*this->var_[ua]*this->var_[ua] * this->var_[ui] * (v_new-E_K));
    
    // ua
    double alpha_ua = 0.65*std::pow((std::exp((v_new+10.)/-8.5)+std::exp((v_new-30.)/-59.)), -1.);
    double beta_ua = 0.65*std::pow((2.5+std::exp((v_new+82.)/17.)), -1.);
    double tau_ua = std::pow((alpha_ua+beta_ua), -1.0)/this->prm_[K_Q10];
    double ua_inf = std::pow((1.+std::exp((v_new+30.03)/-9.6)), -1.);
    this->var_[ua] = ALGORITHM::RushLarsen(ua_inf, this->var_[ua], dt, tau_ua);

    // ui
    double alpha_ui = std::pow((21.+1.*std::exp((v_new-185.)/-28.0)), -1.);
    double beta_ui = 1.0*std::exp((v_new-158.)/(16.));
    double tau_ui = std::pow((alpha_ui+beta_ui),-1.)/this->prm_[K_Q10];
    double ui_inf = std::pow((1.+std::exp((v_new-99.45)/27.48)), -1.);
    this->var_[ui] = ALGORITHM::RushLarsen(ui_inf, this->var_[ui], dt, tau_ui);

    // Rapid Delayed Outward Rectifier K+ Current
    this->cur_[IKr] =  (1.0-this->block_coeff_[IKr]) * (this->prm_[Cm] * g_Kr * this->var_[xr]*(v_new-E_K) / (1.+std::exp((v_new+15.)/22.4)));
    
    // xr
    double alpha_xr = 0.0003*((v_new+14.1)/(1-std::exp((v_new+14.1)/(-5.))));
    double beta_xr = 7.3898e-5*((v_new-3.3328)/(std::exp((v_new-3.3328)/5.1237)-1.));
    double tau_xr = std::pow((alpha_xr+beta_xr), -1.);
    double xr_inf = std::pow((1.+std::exp((v_new+14.1)/-6.5)), -1.);
    this->var_[xr] = ALGORITHM::RushLarsen(xr_inf, this->var_[xr], dt, tau_xr);
            
    // Slow Delayed Outward Rectifier K+ Current
    this->cur_[IKs] =  (1.0-this->block_coeff_[IKs]) * (this->prm_[Cm] * g_Ks * this->var_[xs]*this->var_[xs] * (v_new-E_K));

    double alpha_xs = 4e-5*((v_new-19.9)/(1.-std::exp((v_new-19.9)/(-17))));
    double beta_xs = 3.5e-5*((v_new-19.9)/(std::exp((v_new-19.9)/9)-1));
    double tau_xs = 0.5*std::pow((alpha_xs+beta_xs), -1.);
    double xs_inf = std::pow((1.+std::exp((v_new-19.9)/-12.7)), -0.5);
    this->var_[xs] = ALGORITHM::RushLarsen(xs_inf, this->var_[xs], dt, tau_xs);

    // L-Type Ca2+ Current
    this->cur_[ICaL] =  (1.0-this->block_coeff_[ICaL]) * (this->prm_[Cm] * g_Ca_L * this->var_[d] * this->var_[f] * this->var_[f_Ca] * (v_new-65.));
    
    // d
    double tau_d = (1.-std::exp((v_new+10.)/(-6.24)))/(0.035*(v_new+10.)*(1.+std::exp((v_new+10.)/(-6.24))));
    double d_inf = std::pow((1.+std::exp((v_new+10.)/(-8.))), -1.);
    this->var_[d] = ALGORITHM::RushLarsen(d_inf, this->var_[d], dt, tau_d);
    
    // f
    double tau_f = 9.*std::pow((0.0197*std::exp(-std::pow(0.0337, 2.)*std::pow((v_new+10.), 2.))+0.02), -1.);
    double f_inf = std::pow((1+exp((v_new+28)/6.9)), -1);
    this->var_[f] = ALGORITHM::RushLarsen(f_inf, this->var_[f], dt, tau_f);
    
    // fCa
    double tau_f_Ca = 2.;
    double f_Ca_inf = std::pow((1.0+this->var_[Ca_i]/0.00035), -1.);
    this->var_[f_Ca] = ALGORITHM::RushLarsen(f_Ca_inf, this->var_[f_Ca], dt, tau_f_Ca);

    // Na+-K- Pump Current
    double sigma = 1./7.*(std::exp(this->prm_[Na_o]/67.3)-1.);
    double f_NaK = std::pow((1.+0.1245*std::exp(-0.1*this->prm_[F]*v_new/(this->prm_[R]*this->prm_[T]))+0.0365*sigma*std::exp(-this->prm_[F]*v_new/(this->prm_[R]*this->prm_[T]))), -1.);
    this->cur_[INaK] =  (1.0-this->block_coeff_[INaK]) * (this->prm_[Cm]*this->prm_[i_NaK_max]*f_NaK*(1./(1.+std::pow((this->prm_[Km_Na_i]/this->var_[Na_i]), 1.5)))*(this->prm_[K_o]/(this->prm_[K_o]+this->prm_[Km_K_o])));

    // IKACh
    this->cur_[IKach] =  (1.0-this->block_coeff_[IKach]) * (this->prm_[Cm] * 10./(1.+(9.13652/std::pow((this->prm_[Ach]*0.1), (0.477811)))) * (v_new - E_K) * (0.0517 + 5.0 / (1. + std::exp((v_new + 85.)/5.))));

    // IKCa given by two formulations, Engel or Penaranda.
    double tau_y = (3.*this->prm_[Engel]) + (5.*this->prm_[Penaranda]);
    double KKCa = 0.0007;
    double gKCa_E = this->prm_[g_Kca_eb]*this->prm_[drugvector]*(1.+0.1304*this->prm_[TYPE1])*(1.+4.2174*this->prm_[TYPE3]);
    double gKCa_P = this->prm_[g_Kca_pb]*this->prm_[drugvector];
    double numzss = std::pow(this->var_[Ca_i]/0.0025, 2.); 
    double denomzss = 1. + numzss;
    double yss = (numzss/denomzss)*this->prm_[Engel] + (this->var_[Ca_i]*this->var_[Ca_i]/(KKCa*KKCa+this->var_[Ca_i]*this->var_[Ca_i]))*this->prm_[Penaranda];
    this->var_[y] = ALGORITHM::RushLarsen(yss, this->var_[y], dt, tau_y);


    // IKCa, as in Engel:
    this->cur_[IKca_E] = (1.0-this->block_coeff_[IKca_E]) * ((this->prm_[Cm] * gKCa_E * (this->var_[y]*this->var_[y])*(v_new - E_K))*this->prm_[Engel]);
    // IKCa, as in Penaranda:
    this->cur_[IKca_P] = (1.0-this->block_coeff_[IKca_P]) * ((this->prm_[Cm] * gKCa_P * this->var_[y]*(v_new - E_K))*this->prm_[Penaranda]);


    // Na+/Ca2+ Exchanger Current
    this->cur_[INaCa] =  (1.0-this->block_coeff_[INaCa]) * (this->prm_[Cm]*this->prm_[I_NaCa_max]*(std::exp(this->prm_[gamma]*this->prm_[F]*v_new/(this->prm_[R]*this->prm_[T]))*std::pow(this->var_[Na_i], 3.)*this->prm_[Ca_o]-std::exp((this->prm_[gamma]-1.)*this->prm_[F]*v_new/(this->prm_[R]*this->prm_[T]))*std::pow(this->prm_[Na_o], 3.)*this->var_[Ca_i])/((std::pow(this->prm_[K_mNa], 3.)+std::pow(this->prm_[Na_o], 3.))*(this->prm_[K_mCa]+this->prm_[Ca_o])*(1.+this->prm_[K_sat]*std::exp((this->prm_[gamma]-1.)*v_new*this->prm_[F]/(this->prm_[R]*this->prm_[T])))));

    // Background Currents
    this->cur_[INab] =  (1.0-this->block_coeff_[INab]) * (this->prm_[Cm]*this->prm_[g_B_Na]*(v_new-E_Na));
    this->cur_[ICab] =  (1.0-this->block_coeff_[ICab]) * (this->prm_[Cm]*this->prm_[g_B_Ca]*(v_new-E_Ca));
    double i_B_K = this->prm_[Cm]*this->prm_[g_B_K]*(v_new-E_K);

    // Ca2+ Pump Current
    this->cur_[ICap] =  (1.0-this->block_coeff_[ICap]) * (this->prm_[Cm]*this->prm_[i_CaP_max]*this->var_[Ca_i]/(0.0005+this->var_[Ca_i]));

    // Ca2+ Release Current from JSR
    this->cur_[IRel] =  (1.0-this->block_coeff_[IRel]) * (this->prm_[K_rel]*this->var_[g_u]*this->var_[g_u]*this->var_[g_v]*this->var_[g_w]*(this->var_[Ca_rel]-this->var_[Ca_i]));
    double Fn = 1e-12*this->prm_[V_rel]*this->cur_[IRel] - ((5e-13)/this->prm_[F])*(0.5*this->cur_[ICaL]-0.2*this->cur_[INaCa]);

    // u
    double tau_u = 8.;
    double u_inf = std::pow((1.+std::exp(-(Fn-3.4175e-13)/13.67e-16)), -1.);
    this->var_[g_u] = ALGORITHM::RushLarsen(u_inf, this->var_[g_u], dt, tau_u);
    
    // v
    double tau_v = 1.91+2.09*std::pow((1.+std::exp(-(Fn-3.4175e-13)/13.67e-16)), -1.0);
    double v_inf = 1.0-std::pow((1.+std::exp(-(Fn-6.835e-14)/13.67e-16)), -1.0);
    this->var_[g_v] = ALGORITHM::RushLarsen(v_inf, this->var_[g_v], dt, tau_v);

    // w
    double tau_w = 6.*(1.0-std::exp(-(v_new-7.9)/5.))/((1.+0.3*std::exp(-(v_new-7.9)/5.))*1.0*(v_new-7.9));
    double w_inf = 1.0-std::pow((1.+std::exp(-(v_new-40.)/17.)), -1.);
    this->var_[g_w] = ALGORITHM::RushLarsen(w_inf, this->var_[g_w], dt, tau_w);

    // Transfer Current from NSR to JSR
    this->cur_[Itr] =  (1.0-this->block_coeff_[Itr]) * ((this->var_[Ca_up]-this->var_[Ca_rel])/this->prm_[tau_tr]);

    // Ca2+ Leak Current by the NSR
    this->cur_[IupLeak] =  (1.0-this->block_coeff_[IupLeak]) * (this->prm_[I_up_max]*this->var_[Ca_up]/this->prm_[Ca_up_max]);

    // Ca2+ Uptake Current by the NSR
    this->cur_[Iup] =  (1.0-this->block_coeff_[Iup]) * (this->prm_[I_up_max]/(1.+this->prm_[K_up]/this->var_[Ca_i]));

    // Total ion current
    this->cur_[CourteCur::Iion] = (this->cur_[INa]+this->cur_[IK1]+this->cur_[Ito]+this->cur_[IKur]+this->cur_[IKr]+this->cur_[IKs]+this->cur_[IKach]+this->cur_[IKca_E]+this->cur_[IKca_P]+this->cur_[INab]+this->cur_[ICab]+this->cur_[INaK]+this->cur_[ICap]+this->cur_[INaCa]+this->cur_[ICaL])/this->prm_[Cm];

    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[CourteCur::Iion] - stim_current);

    double dNa_i = (-3.0*this->cur_[INaK]-(3.*this->cur_[INaCa]+this->cur_[INab]+this->cur_[INa])) / (this->prm_[V_i]*this->prm_[F]);
    this->var_[Na_i] = ALGORITHM::ForwardEuler(this->var_[Na_i], dt, dNa_i);

    double dK_i = (2.0*this->cur_[INaK]-(this->cur_[IK1]+this->cur_[Ito]+this->cur_[IKur]+this->cur_[IKr]+this->cur_[IKs]+i_B_K)) / (this->prm_[V_i]*this->prm_[F]);
    this->var_[K_i] = ALGORITHM::ForwardEuler(this->var_[K_i], dt, dK_i);

    double B1 = (2.*this->cur_[INaCa]-(this->cur_[ICap]+this->cur_[ICaL]+this->cur_[ICab]))/(2.*this->prm_[V_i]*this->prm_[F])+(this->prm_[V_up]*(this->cur_[IupLeak]-this->cur_[Iup])+this->cur_[IRel]*this->prm_[V_rel])/this->prm_[V_i];
    double B2 = 1.0+this->prm_[TRPN_max]*this->prm_[Km_TRPN]/std::pow((this->var_[Ca_i]+this->prm_[Km_TRPN]), 2.) + this->prm_[CMDN_max]*this->prm_[Km_CMDN]/std::pow((this->var_[Ca_i]+this->prm_[Km_CMDN]), 2.);
    this->var_[Ca_i] = ALGORITHM::ForwardEuler(this->var_[Ca_i], dt, B1/B2);

    double dCa_up = this->cur_[Iup]-(this->cur_[IupLeak]+this->cur_[Itr]*this->prm_[V_rel]/this->prm_[V_up]);
    this->var_[Ca_up] = ALGORITHM::ForwardEuler(this->var_[Ca_up], dt, dCa_up);

    double dCa_rel = (this->cur_[Itr]-this->cur_[IRel]) * std::pow((1.+this->prm_[CSQN_max]*this->prm_[Km_CSQN]/std::pow((this->var_[Ca_rel]+this->prm_[Km_CSQN]), 2.)), -1.);
    this->var_[Ca_rel] = ALGORITHM::ForwardEuler(this->var_[Ca_rel], dt, dCa_rel);

}



std::string Courtemanche::PrintVariables() const 
{
    using namespace CourteVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "g_u = " << this->var_[g_u] << "\n";
    oss << "g_v = " << this->var_[g_v] << "\n";
    oss << "g_w = " << this->var_[g_w] << "\n";
    oss << "d = " << this->var_[d] << "\n";
    oss << "f_Ca = " << this->var_[f_Ca] << "\n";
    oss << "f = " << this->var_[f] << "\n";
    oss << "h = " << this->var_[h] << "\n";
    oss << "j = " << this->var_[j] << "\n";
    oss << "m = " << this->var_[m] << "\n";
    oss << "Ca_i = " << this->var_[Ca_i] << "\n";
    oss << "Ca_rel = " << this->var_[Ca_rel] << "\n";
    oss << "Ca_up = " << this->var_[Ca_up] << "\n";
    oss << "K_i = " << this->var_[K_i] << "\n";
    oss << "Na_i = " << this->var_[Na_i] << "\n";
    oss << "xr = " << this->var_[xr] << "\n";
    oss << "xs = " << this->var_[xs] << "\n";
    oss << "oa = " << this->var_[oa] << "\n";
    oss << "oi = " << this->var_[oi] << "\n";
    oss << "ua = " << this->var_[ua] << "\n";
    oss << "ui = " << this->var_[ui] << "\n";
    oss << "y = " << this->var_[y];
    return oss.str();

}


std::string Courtemanche::PrintParameters() const
{
    using namespace CourtePrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "rkr = " << this->prm_[rkr] << "\n";
    oss << "TYPE1 = " << this->prm_[TYPE1] << "\n";
    oss << "TYPE3 = " << this->prm_[TYPE3] << "\n";
    oss << "ISOdlf = " << this->prm_[ISOdlf] << "\n";
    oss << "cAF = " << this->prm_[cAF] << "\n";
    oss << "cAF2 = " << this->prm_[cAF2] << "\n";
    oss << "SR = " << this->prm_[SR] << "\n";
    oss << "Ach = " << this->prm_[Ach] << "\n";
    oss << "drugvector = " << this->prm_[drugvector] << "\n";
    oss << "Engel = " << this->prm_[Engel] << "\n";
    oss << "Penaranda = " << this->prm_[Penaranda] << "\n";
    oss << "Cm = " << this->prm_[Cm] << "\n";
    oss << "CMDN_max = " << this->prm_[CMDN_max] << "\n";
    oss << "CSQN_max = " << this->prm_[CSQN_max] << "\n";
    oss << "Km_CMDN = " << this->prm_[Km_CMDN] << "\n";
    oss << "Km_CSQN = " << this->prm_[Km_CSQN] << "\n";
    oss << "Km_TRPN = " << this->prm_[Km_TRPN] << "\n";
    oss << "TRPN_max = " << this->prm_[TRPN_max] << "\n";
    oss << "Ca_up_max = " << this->prm_[Ca_up_max] << "\n";
    oss << "K_rel = " << this->prm_[K_rel] << "\n";
    oss << "I_up_max = " << this->prm_[I_up_max] << "\n";
    oss << "K_up = " << this->prm_[K_up] << "\n";
    oss << "I_NaCa_max = " << this->prm_[I_NaCa_max] << "\n";
    oss << "K_mCa = " << this->prm_[K_mCa] << "\n";
    oss << "K_mNa = " << this->prm_[K_mNa] << "\n";
    oss << "K_sat = " << this->prm_[K_sat] << "\n";
    oss << "gamma = " << this->prm_[gamma] << "\n";
    oss << "g_B_Ca = " << this->prm_[g_B_Ca] << "\n";
    oss << "g_B_K = " << this->prm_[g_B_K] << "\n";
    oss << "g_B_Na = " << this->prm_[g_B_Na] << "\n";
    oss << "g_Na = " << this->prm_[g_Na] << "\n";
    oss << "g_Kca_eb = " << this->prm_[g_Kca_eb] << "\n";
    oss << "g_Kca_pb = " << this->prm_[g_Kca_pb] << "\n";
    oss << "V_cell = " << this->prm_[V_cell] << "\n";
    oss << "V_i = " << this->prm_[V_i] << "\n";
    oss << "V_up = " << this->prm_[V_up] << "\n";
    oss << "V_rel = " << this->prm_[V_rel] << "\n";
    oss << "F = " << this->prm_[F] << "\n";
    oss << "R = " << this->prm_[R] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "i_CaP_max = " << this->prm_[i_CaP_max] << "\n";
    oss << "Km_K_o = " << this->prm_[Km_K_o] << "\n";
    oss << "Km_Na_i = " << this->prm_[Km_Na_i] << "\n";
    oss << "i_NaK_max = " << this->prm_[i_NaK_max] << "\n";
    oss << "Ca_o = " << this->prm_[Ca_o] << "\n";
    oss << "K_o = " << this->prm_[K_o] << "\n";
    oss << "Na_o = " << this->prm_[Na_o] << "\n";
    oss << "tau_tr = " << this->prm_[tau_tr] << "\n";
    oss << "K_Q10 = " << this->prm_[K_Q10] << "\n";
    oss << "EAD = " << this->prm_[EAD];
    return oss.str();

}


std::string Courtemanche::PrintCurrents() const
{
    using namespace CourteCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "IKur = " << this->cur_[IKur] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "IKach = " << this->cur_[IKach] << "\n";
    oss << "IKca_E = " << this->cur_[IKca_E] << "\n";
    oss << "IKca_P = " << this->cur_[IKca_P] << "\n";
    oss << "INaCa = " << this->cur_[INaCa] << "\n";
    oss << "INab = " << this->cur_[INab] << "\n";
    oss << "ICab = " << this->cur_[ICab] << "\n";
    oss << "ICap = " << this->cur_[ICap] << "\n";
    oss << "IRel = " << this->cur_[IRel] << "\n";
    oss << "Itr = " << this->cur_[Itr] << "\n";
    oss << "IupLeak = " << this->cur_[IupLeak] << "\n";
    oss << "Iup = " << this->cur_[Iup] << "\n";
    oss << "Iion = " << this->cur_[CourteCur::Iion]; 
    return oss.str();

}


std::string Courtemanche::PrintBlockCoeffs() const
{
    using namespace CourteCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->block_coeff_[INa] << "\n";
    oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
    oss << "Ito = " << this->block_coeff_[Ito] << "\n";
    oss << "IKur = " << this->block_coeff_[IKur] << "\n";
    oss << "IKr = " << this->block_coeff_[IKr] << "\n";
    oss << "IKs = " << this->block_coeff_[IKs] << "\n";
    oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "IKach = " << this->block_coeff_[IKach] << "\n";
    oss << "IKca_E = " << this->block_coeff_[IKca_E] << "\n";
    oss << "IKca_P = " << this->block_coeff_[IKca_P] << "\n";
    oss << "INaCa = " << this->block_coeff_[INaCa] << "\n";
    oss << "INab = " << this->block_coeff_[INab] << "\n";
    oss << "ICab = " << this->block_coeff_[ICab] << "\n";
    oss << "ICap = " << this->block_coeff_[ICap] << "\n";
    oss << "IRel = " << this->block_coeff_[IRel] << "\n";
    oss << "Itr = " << this->block_coeff_[Itr] << "\n";
    oss << "IupLeak = " << this->block_coeff_[IupLeak] << "\n";
    oss << "Iup = " << this->block_coeff_[Iup];
    return oss.str();

}



} // End of namespace ELECTRA