/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/electrophysiology/paci_ventri.hpp"


namespace ELECTRA {

void PaciVentri::SetDataMapping()
{
    using namespace PciVtrVar;
    using namespace PciVtrPrm;
    using namespace PciVtrCur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["Nai"] = static_cast<std::size_t>(Nai);
    this->mapped_data_["Cai"] = static_cast<std::size_t>(Cai);
    this->mapped_data_["m"] = static_cast<std::size_t>(m);
    this->mapped_data_["h"] = static_cast<std::size_t>(h);
    this->mapped_data_["j"] = static_cast<std::size_t>(j);
    this->mapped_data_["d"] = static_cast<std::size_t>(d);
    this->mapped_data_["f1"] = static_cast<std::size_t>(f1);
    this->mapped_data_["f2"] = static_cast<std::size_t>(f2);
    this->mapped_data_["fCa"] = static_cast<std::size_t>(fCa);
    this->mapped_data_["Xr1"] = static_cast<std::size_t>(Xr1);
    this->mapped_data_["Xr2"] = static_cast<std::size_t>(Xr2);
    this->mapped_data_["Xs"] = static_cast<std::size_t>(Xs);
    this->mapped_data_["Xf"] = static_cast<std::size_t>(Xf);
    this->mapped_data_["q"] = static_cast<std::size_t>(q);
    this->mapped_data_["r"] = static_cast<std::size_t>(r);
    this->mapped_data_["Ca_SR"] = static_cast<std::size_t>(Ca_SR);
    this->mapped_data_["g"] = static_cast<std::size_t>(g);

    // Set parameters mapping.
    this->mapped_data_["Cm"] = static_cast<std::size_t>(Cm);
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["F"] = static_cast<std::size_t>(F);
    this->mapped_data_["Nao"] = static_cast<std::size_t>(Nao);
    this->mapped_data_["Cao"] = static_cast<std::size_t>(Cao);
    this->mapped_data_["Ki"] = static_cast<std::size_t>(Ki);
    this->mapped_data_["Ko"] = static_cast<std::size_t>(Ko);
    this->mapped_data_["PkNa"] = static_cast<std::size_t>(PkNa);
    this->mapped_data_["g_Na"] = static_cast<std::size_t>(g_Na);
    this->mapped_data_["g_CaL"] = static_cast<std::size_t>(g_CaL);
    this->mapped_data_["tau_fCa"] = static_cast<std::size_t>(tau_fCa);
    this->mapped_data_["g_Kr"] = static_cast<std::size_t>(g_Kr);
    this->mapped_data_["L0"] = static_cast<std::size_t>(L0);
    this->mapped_data_["Q"] = static_cast<std::size_t>(Q);
    this->mapped_data_["g_Ks"] = static_cast<std::size_t>(g_Ks);
    this->mapped_data_["g_K1"] = static_cast<std::size_t>(g_K1);
    this->mapped_data_["g_f"] = static_cast<std::size_t>(g_f);
    this->mapped_data_["E_f"] = static_cast<std::size_t>(E_f);
    this->mapped_data_["g_b_Na"] = static_cast<std::size_t>(g_b_Na);
    this->mapped_data_["g_b_Ca"] = static_cast<std::size_t>(g_b_Ca);
    this->mapped_data_["Km_K"] = static_cast<std::size_t>(Km_K);
    this->mapped_data_["Km_Na"] = static_cast<std::size_t>(Km_Na);
    this->mapped_data_["PNaK"] = static_cast<std::size_t>(PNaK);
    this->mapped_data_["kNaCa"] = static_cast<std::size_t>(kNaCa);
    this->mapped_data_["alpha"] = static_cast<std::size_t>(alpha);
    this->mapped_data_["gamma"] = static_cast<std::size_t>(gamma);
    this->mapped_data_["Ksat"] = static_cast<std::size_t>(Ksat);
    this->mapped_data_["KmCa"] = static_cast<std::size_t>(KmCa);
    this->mapped_data_["KmNai"] = static_cast<std::size_t>(KmNai);
    this->mapped_data_["g_PCa"] = static_cast<std::size_t>(g_PCa);
    this->mapped_data_["KPCa"] = static_cast<std::size_t>(KPCa);
    this->mapped_data_["g_to"] = static_cast<std::size_t>(g_to);
    this->mapped_data_["Vc"] = static_cast<std::size_t>(Vc);
    this->mapped_data_["V_SR"] = static_cast<std::size_t>(V_SR);
    this->mapped_data_["a_rel"] = static_cast<std::size_t>(a_rel);
    this->mapped_data_["b_rel"] = static_cast<std::size_t>(b_rel);
    this->mapped_data_["c_rel"] = static_cast<std::size_t>(c_rel);
    this->mapped_data_["tau_g"] = static_cast<std::size_t>(tau_g);
    this->mapped_data_["Kup"] = static_cast<std::size_t>(Kup);
    this->mapped_data_["Buf_C"] = static_cast<std::size_t>(Buf_C);
    this->mapped_data_["Buf_SR"] = static_cast<std::size_t>(Buf_SR);
    this->mapped_data_["Kbuf_C"] = static_cast<std::size_t>(Kbuf_C);
    this->mapped_data_["Kbuf_SR"] = static_cast<std::size_t>(Kbuf_SR);
    this->mapped_data_["VmaxUp"] = static_cast<std::size_t>(VmaxUp);
    this->mapped_data_["V_leak"] = static_cast<std::size_t>(V_leak);
    this->mapped_data_["E_K"] = static_cast<std::size_t>(E_K);
    this->mapped_data_["constf2"] = static_cast<std::size_t>(constf2);
    this->mapped_data_["V_half"] = static_cast<std::size_t>(V_half);
    this->mapped_data_["m_inf"] = static_cast<std::size_t>(m_inf);
    this->mapped_data_["h_inf"] = static_cast<std::size_t>(h_inf);
    this->mapped_data_["j_inf"] = static_cast<std::size_t>(j_inf);
    this->mapped_data_["d_inf"] = static_cast<std::size_t>(d_inf);
    this->mapped_data_["f1_inf"] = static_cast<std::size_t>(f1_inf);
    this->mapped_data_["f2_inf"] = static_cast<std::size_t>(f2_inf);
    this->mapped_data_["alpha_fCa"] = static_cast<std::size_t>(alpha_fCa);
    this->mapped_data_["Xr1_inf"] = static_cast<std::size_t>(Xr1_inf);
    this->mapped_data_["Xr2_inf"] = static_cast<std::size_t>(Xr2_inf);
    this->mapped_data_["Xs_inf"] = static_cast<std::size_t>(Xs_inf);
    this->mapped_data_["Xf_inf"] = static_cast<std::size_t>(Xf_inf);
    this->mapped_data_["q_inf"] = static_cast<std::size_t>(q_inf);
    this->mapped_data_["r_inf"] = static_cast<std::size_t>(r_inf);
    this->mapped_data_["g_inf"] = static_cast<std::size_t>(g_inf);
    this->mapped_data_["alpha_m"] = static_cast<std::size_t>(alpha_m);
    this->mapped_data_["alpha_h"] = static_cast<std::size_t>(alpha_h);
    this->mapped_data_["alpha_j"] = static_cast<std::size_t>(alpha_j);
    this->mapped_data_["alpha_d"] = static_cast<std::size_t>(alpha_d);
    this->mapped_data_["constf1"] = static_cast<std::size_t>(constf1);
    this->mapped_data_["tau_f2"] = static_cast<std::size_t>(tau_f2);
    this->mapped_data_["beta_fCa"] = static_cast<std::size_t>(beta_fCa);
    this->mapped_data_["alpha_Xr1"] = static_cast<std::size_t>(alpha_Xr1);
    this->mapped_data_["alpha_Xr2"] = static_cast<std::size_t>(alpha_Xr2);
    this->mapped_data_["alpha_Xs"] = static_cast<std::size_t>(alpha_Xs);
    this->mapped_data_["E_Na"] = static_cast<std::size_t>(E_Na);
    this->mapped_data_["tau_Xf"] = static_cast<std::size_t>(tau_Xf);
    this->mapped_data_["tau_q"] = static_cast<std::size_t>(tau_q);
    this->mapped_data_["tau_r"] = static_cast<std::size_t>(tau_r);
    this->mapped_data_["const2"] = static_cast<std::size_t>(const2);
    this->mapped_data_["beta_m"] = static_cast<std::size_t>(beta_m);
    this->mapped_data_["beta_h"] = static_cast<std::size_t>(beta_h);
    this->mapped_data_["beta_j"] = static_cast<std::size_t>(beta_j);
    this->mapped_data_["beta_d"] = static_cast<std::size_t>(beta_d);
    this->mapped_data_["tau_f1"] = static_cast<std::size_t>(tau_f1);
    this->mapped_data_["gamma_fCa"] = static_cast<std::size_t>(gamma_fCa);
    this->mapped_data_["beta_Xr1"] = static_cast<std::size_t>(beta_Xr1);
    this->mapped_data_["beta_Xr2"] = static_cast<std::size_t>(beta_Xr2);
    this->mapped_data_["beta_Xs"] = static_cast<std::size_t>(beta_Xs);
    this->mapped_data_["E_Ks"] = static_cast<std::size_t>(E_Ks);
    this->mapped_data_["tau_m"] = static_cast<std::size_t>(tau_m);
    this->mapped_data_["tau_h"] = static_cast<std::size_t>(tau_h);
    this->mapped_data_["tau_j"] = static_cast<std::size_t>(tau_j);
    this->mapped_data_["gamma_d"] = static_cast<std::size_t>(gamma_d);
    this->mapped_data_["fCa_inf"] = static_cast<std::size_t>(fCa_inf);
    this->mapped_data_["tau_Xr1"] = static_cast<std::size_t>(tau_Xr1);
    this->mapped_data_["tau_Xr2"] = static_cast<std::size_t>(tau_Xr2);
    this->mapped_data_["tau_Xs"] = static_cast<std::size_t>(tau_Xs);
    this->mapped_data_["E_Ca"] = static_cast<std::size_t>(E_Ca);
    this->mapped_data_["tau_d"] = static_cast<std::size_t>(tau_d);
    this->mapped_data_["constfCa"] = static_cast<std::size_t>(constfCa);
    this->mapped_data_["alpha_K1"] = static_cast<std::size_t>(alpha_K1);
    this->mapped_data_["beta_K1"] = static_cast<std::size_t>(beta_K1);
    this->mapped_data_["XK1_inf"] = static_cast<std::size_t>(XK1_inf);
    this->mapped_data_["Cai_bufc"] = static_cast<std::size_t>(Cai_bufc);
    this->mapped_data_["Ca_SR_bufSR"] = static_cast<std::size_t>(Ca_SR_bufSR);
    this->mapped_data_["Xf_infinity"] = static_cast<std::size_t>(Xf_infinity);

    // Set currents mapping.
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["ICaL"] = static_cast<std::size_t>(ICaL);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);
    this->mapped_data_["If"] = static_cast<std::size_t>(If);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["Ito"] = static_cast<std::size_t>(Ito);
    this->mapped_data_["IPCa"] = static_cast<std::size_t>(IPCa);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);
    this->mapped_data_["Irel"] = static_cast<std::size_t>(Irel);
    this->mapped_data_["Iup"] = static_cast<std::size_t>(Iup);
    this->mapped_data_["Ileak"] = static_cast<std::size_t>(Ileak);

}


PaciVentri::PaciVentri()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::PaciVentri;
    this->dt_stable_ = 0.02;
    this->var_.resize(19, 0.);
    this->prm_.resize(105, 0.);
    this->cur_.resize(16, 0.);
    this->block_coeff_.resize(15, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


PaciVentri::~PaciVentri()
{}


void PaciVentri::Initialize(CellType cell_type)
{
    using namespace PciVtrVar;
    using namespace PciVtrPrm;

    if (cell_type != CellType::ventricular) {
        std::string error_str = "Could not initialize PaciVentri ap model. Expected: CellType::ventricular";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(19, 0.);
    this->prm_.clear();           this->prm_.resize(105, 0.);
    this->cur_.clear();           this->cur_.resize(16, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(15, 0.);

    // Set the model variables.
    this->var_[v] = -74.3340057623841;
    this->var_[dvdt] = 0.;
    this->var_[Nai] = 10.9248496211574;     
    this->var_[Cai] = 1.80773974140477e-5;  
    this->var_[m] = 0.102953468725004;      
    this->var_[h] = 0.786926637881461;      
    this->var_[j] = 0.253943221774722;      
    this->var_[d] = 8.96088425225182e-5;    
    this->var_[f1] = 0.970411811263976;     
    this->var_[f2] = 0.999965815466749;     
    this->var_[fCa] = 0.998925296531804;    
    this->var_[Xr1] = 0.00778547011240132;  
    this->var_[Xr2] = 0.432162576531617;    
    this->var_[Xs] = 0.0322944866983666;    
    this->var_[Xf] = 0.100615100568753;     
    this->var_[q] = 0.839295925773219;      
    this->var_[r] = 0.00573289893326379;    
    this->var_[Ca_SR] = 0.2734234751931;    
    this->var_[g] = 0.999999981028517;        

    // Set the model parameters.            
    this->prm_[Cm] = 9.87109e-11;               
    this->prm_[R] = 8.314472;                   
    this->prm_[T] = 310;                        
    this->prm_[F] = 96485.3415;                 
    this->prm_[Nao] = 151;                      
    this->prm_[Cao] = 1.8;                      
    this->prm_[Ki] = 150;                       
    this->prm_[Ko] = 5.4;                       
    this->prm_[PkNa] = 0.03;                    
    this->prm_[g_Na] = 3671.2302;               
    this->prm_[g_CaL] = 8.635702e-5;            
    this->prm_[tau_fCa] = 0.002;                
    this->prm_[g_Kr] = 29.8667;                 
    this->prm_[L0] = 0.025;                     
    this->prm_[Q] = 2.3;                        
    this->prm_[g_Ks] = 2.041;                   
    this->prm_[g_K1] = 28.1492;                 
    this->prm_[g_f] = 30.10312;                 
    this->prm_[E_f] = -0.017;                   
    this->prm_[g_b_Na] = 0.9;                   
    this->prm_[g_b_Ca] = 0.69264;               
    this->prm_[Km_K] = 1;                       
    this->prm_[Km_Na] = 40;                     
    this->prm_[PNaK] = 1.841424;                
    this->prm_[kNaCa] = 4900;                   
    this->prm_[alpha] = 2.8571432;              
    this->prm_[gamma] = 0.35;                   
    this->prm_[Ksat] = 0.1;                     
    this->prm_[KmCa] = 1.38;                    
    this->prm_[KmNai] = 87.5;                   
    this->prm_[g_PCa] = 0.4125;                 
    this->prm_[KPCa] = 0.0005;                  
    this->prm_[g_to] = 29.9038;                 
    this->prm_[Vc] = 8800;                      
    this->prm_[V_SR] = 583.73;                  
    this->prm_[a_rel] = 16.464;                 
    this->prm_[b_rel] = 0.25;                   
    this->prm_[c_rel] = 8.232;                  
    this->prm_[tau_g] = 0.002;                  
    this->prm_[Kup] = 0.00025;                  
    this->prm_[Buf_C] = 0.25;                   
    this->prm_[Buf_SR] = 10;                    
    this->prm_[Kbuf_C] = 0.001;                 
    this->prm_[Kbuf_SR] = 0.3;                  
    this->prm_[VmaxUp] = 0.56064;               
    this->prm_[V_leak] = 0.00044444;            
    this->prm_[E_K] = ((this->prm_[R]*this->prm_[T])/this->prm_[F])*std::log(this->prm_[Ko]/this->prm_[Ki]);
    this->prm_[constf2] = 1;                    
    this->prm_[V_half] = 1000*(((-this->prm_[R]*this->prm_[T])/(this->prm_[F]*this->prm_[Q]))*std::log(std::pow(1+this->prm_[Cao]/2.6, 4)/(this->prm_[L0]*std::pow(1+this->prm_[Cao]/0.58, 4))) - 0.019);
    this->prm_[m_inf] = 0;              
    this->prm_[h_inf] = 0;              
    this->prm_[j_inf] = 0;              
    this->prm_[d_inf] = 0;              
    this->prm_[f1_inf] = 0;             
    this->prm_[f2_inf] = 0;             
    this->prm_[alpha_fCa] = 0;          
    this->prm_[Xr1_inf] = 0;            
    this->prm_[Xr2_inf] = 0;            
    this->prm_[Xs_inf] = 0;             
    this->prm_[Xf_inf] = 0;             
    this->prm_[q_inf] = 0;              
    this->prm_[r_inf] = 0;              
    this->prm_[g_inf] = 0;              
    this->prm_[alpha_m] = 0;            
    this->prm_[alpha_h] = 0;            
    this->prm_[alpha_j] = 0;            
    this->prm_[alpha_d] = 0;            
    this->prm_[constf1] = 0;            
    this->prm_[tau_f2] = 0;             
    this->prm_[beta_fCa] = 0;           
    this->prm_[alpha_Xr1] = 0;          
    this->prm_[alpha_Xr2] = 0;          
    this->prm_[alpha_Xs] = 0;           
    this->prm_[E_Na] = 0;               
    this->prm_[tau_Xf] = 0;             
    this->prm_[tau_q] = 0;              
    this->prm_[tau_r] = 0;              
    this->prm_[const2] = 0;             
    this->prm_[beta_m] = 0;             
    this->prm_[beta_h] = 0;                 
    this->prm_[beta_j] = 0;             
    this->prm_[beta_d] = 0;             
    this->prm_[tau_f1] = 0;             
    this->prm_[gamma_fCa] = 0;          
    this->prm_[beta_Xr1] = 0;           
    this->prm_[beta_Xr2] = 0;           
    this->prm_[beta_Xs] = 0;              
    this->prm_[E_Ks] = 0;               
    this->prm_[tau_m] = 0;                  
    this->prm_[tau_h] = 0;              
    this->prm_[tau_j] = 0;              
    this->prm_[gamma_d] = 0;            
    this->prm_[fCa_inf] = 0;            
    this->prm_[tau_Xr1] = 0;            
    this->prm_[tau_Xr2] = 0;            
    this->prm_[tau_Xs] = 0;             
    this->prm_[E_Ca] = 0;               
    this->prm_[tau_d] = 0;              
    this->prm_[constfCa] = 0;           
    this->prm_[alpha_K1] = 0;           
    this->prm_[beta_K1] = 0;            
    this->prm_[XK1_inf] = 0;            
    this->prm_[Cai_bufc] = 0;           
    this->prm_[Ca_SR_bufSR] = 0;        
    this->prm_[Xf_infinity] = 0;

}


void PaciVentri::Compute(double v_new, double dt, double stim_current)
{
    using namespace PciVtrVar;
    using namespace PciVtrPrm;
    using namespace PciVtrCur;

    // Convert to SI.
    dt *= 1.e-3;
    v_new *= 1.e-3;
    stim_current *= 1.e-9;

    this->prm_[f2_inf] = 0.33 + 0.67/(1. + std::exp((1000.*v_new+35.) / 4.));
    this->prm_[tau_f2] = ((600.*std::exp(-std::pow(1000.*v_new+25., 2.)/170.) + (31./(1. + std::exp((25.-1000.*v_new)/10.)) + 16./(1. + std::exp((30.+1000.*v_new)/10.)))) * this->prm_[constf2]) / 1000.;
    this->var_[f2] = this->prm_[f2_inf] - (this->prm_[f2_inf] - this->var_[f2])*std::exp(-dt/this->prm_[tau_f2]);
    
    this->prm_[Xf_infinity] = 1. / (1. + std::exp((1000.*v_new+77.85)/5.));
    this->prm_[tau_Xf] = (1900. / (1. + std::exp((1000.*v_new+15.)/10.))) / 1000.;
    this->var_[Xf] = this->var_[Xf] + dt*(this->prm_[Xf_infinity] - this->var_[Xf]) / this->prm_[tau_Xf];
    
    this->prm_[q_inf] = 1. / (1. + std::exp((1000.*v_new+53.) / 13.));
    this->prm_[tau_q] = (6.06 + 39.102/(0.57*exp(-0.08*(1000.*v_new+44.)) + 0.065*std::exp(0.1*(1000.*v_new+45.93)))) / 1000.;
    this->var_[q] = this->prm_[q_inf] - (this->prm_[q_inf] - this->var_[q]) * std::exp(-dt/this->prm_[tau_q]);
    
    this->prm_[r_inf] = 1. / (1. + std::exp(-(1000.*v_new - 22.3) / 18.75));
    this->prm_[tau_r] = (2.75352 + 14.4052/(1.037*std::exp(0.09*(1000.*v_new+30.61)) + 0.369*std::exp(-0.12*(1000.*v_new+23.84))))/1000.;
    this->var_[r] = this->prm_[r_inf] - (this->prm_[r_inf]-this->var_[r])*std::exp(-dt/this->prm_[tau_r]);
    
    if (this->var_[Cai] <= 0.00035) {
        this->prm_[g_inf] = 1. / (1. + std::pow(this->var_[Cai]/0.00035, 6.));
    } else {
        this->prm_[g_inf] = 1. / (1. + std::pow(this->var_[Cai]/0.00035, 16.));
    }
        
    if (this->prm_[g_inf] > this->var_[g] && v_new > -0.06) {
        this->prm_[const2] = 0.;
    } else {
        this->prm_[const2] = 1.;
    }
    this->var_[g] = this->var_[g] + dt*(this->prm_[const2]*(this->prm_[g_inf] - this->var_[g])) / this->prm_[tau_g];
    
    this->prm_[f1_inf] = 1. / (1. + std::exp((1000.*v_new+26.)/3.));
    
    if (this->prm_[f1_inf] - this->var_[f1] > 0.) {
        this->prm_[constf1] = 1. + 1433.*(this->var_[Cai] - 50.*1e-6);
    } else {
        this->prm_[constf1] = 1.;
    }
    this->prm_[tau_f1] = ((20. + (1102.5*std::exp(-std::pow(std::pow(1000.*v_new+27., 2.)/15., 2.)) + (200./(1.+std::exp((13.-1000.*v_new)/10.)) + 180./(1. + std::exp((30.+1000.*v_new)/10.)))))*this->prm_[constf1])/1000.;
    this->var_[f1] = this->prm_[f1_inf] - (this->prm_[f1_inf]-this->var_[f1])*std::exp(-dt/this->prm_[tau_f1]);
    
    this->prm_[m_inf] = 1. / std::pow( 1 +std::exp((-1000.*v_new - 34.1)/5.9), 1./3.);
    this->prm_[alpha_m] = 1. / (1. + std::exp((-1000.*v_new - 60.)/5.));
    this->prm_[beta_m] = 0.1/(1. + std::exp((1000.*v_new + 35.)/5.)) + 0.1/(1+std::exp((1000.*v_new - 50.)/200.));
    this->prm_[tau_m] = (this->prm_[alpha_m]*this->prm_[beta_m]) / 1000.;
    this->var_[m] = this->prm_[m_inf] - (this->prm_[m_inf]-this->var_[m])*std::exp(-dt/this->prm_[tau_m]);
    
    this->prm_[h_inf] = 1. / std::pow((1.+std::exp((1000.*v_new+72.1)/5.7)), 0.5);
    if (v_new < -0.04) {
        this->prm_[alpha_h] = 0.057*std::exp(-(1000.*v_new+80.)/6.8);
        this->prm_[beta_h] = 2.7*std::exp(0.079*(1000.*v_new)) + 3.1*(std::pow(10., 5.)*std::exp(0.3485*(1000.*v_new)));
        this->prm_[tau_h] = 1.5 / ((this->prm_[alpha_h]+this->prm_[beta_h])*1000.);
    } else {
        this->prm_[alpha_h] = 0.;
        this->prm_[beta_h] = 0.77 / (0.13*(1. + std::exp((1000.*v_new+10.66)/ - 11.1)));
        this->prm_[tau_h] = 2.542/1000.;
    }
    this->var_[h] = this->prm_[h_inf] - (this->prm_[h_inf]-this->var_[h])*std::exp(-dt/this->prm_[tau_h]);
    
    this->prm_[j_inf] = 1. / std::pow((1. + std::exp((1000.*v_new+72.1)/5.7)), 0.5);
    if (v_new < -0.04) {
        this->prm_[alpha_j] = ((-25428.*std::exp(0.2444*(1000.*v_new)) - 6.948*(std::pow(10., -6.)*std::exp(-0.04391*(1000.*v_new)))) * (1000.*v_new+37.78)) / (1.+std::exp(0.311*(1000.*v_new+79.23)));
        this->prm_[beta_j] = (0.02424*std::exp(-0.01052*(1000.*v_new))) / (1.+std::exp(-0.1378*(1000.*v_new+40.14)));
    } else {
        this->prm_[alpha_j] = 0.;
        this->prm_[beta_j] = (0.6*std::exp(0.057*(1000.*v_new))) / (1+std::exp(-0.1*(1000.*v_new+32.)));
    }
    this->prm_[tau_j] = 7. / (1000.*(this->prm_[alpha_j]+this->prm_[beta_j]));
    this->var_[j] = this->prm_[j_inf] - (this->prm_[j_inf]-this->var_[j])*std::exp(-dt/this->prm_[tau_j]);
    
    this->prm_[Xr1_inf] = 1. / (1. + std::exp((this->prm_[V_half] - 1000.*v_new)/4.9));
    this->prm_[alpha_Xr1] = 450. / (1. + std::exp((-45. - 1000.*v_new)/10.));
    this->prm_[beta_Xr1] = 6. / (1. + std::exp((30. + 1000.*v_new)/11.));
    this->prm_[tau_Xr1] = (this->prm_[alpha_Xr1]*this->prm_[beta_Xr1])/1000.;
    this->var_[Xr1] = this->prm_[Xr1_inf] - (this->prm_[Xr1_inf]-this->var_[Xr1]) * std::exp(-dt/this->prm_[tau_Xr1]);
    
    this->prm_[Xr2_inf] = 1. / (1. + std::exp((1000.*v_new+88.)/50.));
    this->prm_[alpha_Xr2] = 3. / (1. + std::exp((-60. - 1000.*v_new)/20.));
    this->prm_[beta_Xr2] = 1.12 / (1. + std::exp((-60. + 1000.*v_new)/20.));
    this->prm_[tau_Xr2] = (this->prm_[alpha_Xr2]*this->prm_[beta_Xr2])/1000.;
    this->var_[Xr2] = this->prm_[Xr2_inf] - (this->prm_[Xr2_inf]-this->var_[Xr2])*std::exp(-dt/this->prm_[tau_Xr2]);
    
    this->prm_[Xs_inf] = 1. / (1. + std::exp((-1000.*v_new - 20.)/16.));
    this->prm_[alpha_Xs] = 1100. / std::pow((1. + std::exp((-10. - 1000.*v_new)/6.)), 0.5);
    this->prm_[beta_Xs] = 1. / (1. + std::exp((-60. + 1000.*v_new)/20.));
    this->prm_[tau_Xs] = (this->prm_[alpha_Xs]*this->prm_[beta_Xs])/1000.;
    this->var_[Xs] = this->prm_[Xs_inf] - (this->prm_[Xs_inf]-this->var_[Xs])*std::exp(-dt/this->prm_[tau_Xs]);
    
    this->prm_[d_inf] = 1. / (1. + std::exp(-(1000.*v_new + 9.1) / 7.));
    this->prm_[alpha_d] = 0.25 + 1.4/(1+std::exp((-1000.*v_new - 35.) / 13.));
    this->prm_[beta_d] = 1.4 / (1+std::exp((1000.*v_new + 5.) / 5.));
    this->prm_[gamma_d] = 1. / (1+std::exp((-1000*v_new + 50.) / 20.));
    this->prm_[tau_d] = (this->prm_[alpha_d]*this->prm_[beta_d]+this->prm_[gamma_d]) / 1000.;
    this->var_[d] = this->prm_[d_inf] - (this->prm_[d_inf]-this->var_[d])*std::exp(-dt/this->prm_[tau_d]);
    
    this->prm_[alpha_fCa] = 1. / (1. + std::pow(this->var_[Cai]/0.0006, 8.));
    this->prm_[beta_fCa] = 0.1 / (1. + std::exp((this->var_[Cai] - 0.0009) / 0.0001));
    this->prm_[gamma_fCa] = 0.3 / (1. + std::exp((this->var_[Cai] - 0.00075) / 0.0008));
    this->prm_[fCa_inf] = (this->prm_[alpha_fCa] + (this->prm_[beta_fCa]+this->prm_[gamma_fCa])) / 1.3156;
    if (v_new > - 0.06 && this->prm_[fCa_inf] > this->var_[fCa]) {
        this->prm_[constfCa] = 0.;
    } else {
        this->prm_[constfCa] = 1.;
    }
    this->var_[fCa] = this->prm_[fCa_inf] - (this->prm_[fCa_inf]-this->var_[fCa])*std::exp(-dt/this->prm_[tau_fCa]);
    
    this->prm_[E_Na] = ((this->prm_[R]*this->prm_[T])/this->prm_[F]) * std::log(this->prm_[Nao]/this->var_[Nai]);
    
    this->cur_[INa] = (1-this->block_coeff_[INa]) * (this->prm_[g_Na]*(std::pow(this->var_[m], 3) * (this->var_[h]*(this->var_[j]*(v_new-this->prm_[E_Na])))));
    this->cur_[INaK] = (1-this->block_coeff_[INaK]) * (((((this->prm_[PNaK]*this->prm_[Ko])/(this->prm_[Ko]+this->prm_[Km_K])) * this->var_[Nai]) / (this->var_[Nai] + this->prm_[Km_Na])) / (1.+(0.1245*std::exp((-0.1*(v_new*this->prm_[F]))/(this->prm_[R]*this->prm_[T])) + 0.0353*std::exp((-v_new*this->prm_[F])/(this->prm_[R]*this->prm_[T])))));
    this->cur_[INaCa] = (1-this->block_coeff_[INaCa]) * ((this->prm_[kNaCa]*(std::exp((this->prm_[gamma]*(v_new*this->prm_[F])) / (this->prm_[R]*this->prm_[T]))*(std::pow(this->var_[Nai], 3.)*this->prm_[Cao]) - std::exp(((this->prm_[gamma] - 1.)*(v_new*this->prm_[F]))/(this->prm_[R]*this->prm_[T]))*(std::pow(this->prm_[Nao], 3.) * (this->var_[Cai]*this->prm_[alpha])))) / ((std::pow(this->prm_[KmNai], 3.) + std::pow(this->prm_[Nao], 3.))*((this->prm_[KmCa]+this->prm_[Cao])*(1. + this->prm_[Ksat]*std::exp(((this->prm_[gamma] - 1.)*(v_new*this->prm_[F])) / (this->prm_[R]*this->prm_[T]))))));
    this->cur_[INab] = (1-this->block_coeff_[INab]) * (this->prm_[g_b_Na]*(v_new - this->prm_[E_Na]));
    this->var_[Nai] = this->var_[Nai] + dt*(-this->prm_[Cm]*(this->cur_[INa]+ this->cur_[INab] + 3.*this->cur_[INaK] + 3.*this->cur_[INaCa])) / (this->prm_[F]*(this->prm_[Vc]*1.e-18));
    
    this->cur_[ICaL] = (1.-this->block_coeff_[ICaL]) * (((((this->prm_[g_CaL] * (4.*(v_new*std::pow(this->prm_[F], 2.)))) / (this->prm_[R]*this->prm_[T])) * (this->var_[Cai]*std::exp((2.*(v_new*this->prm_[F])) / (this->prm_[R]*this->prm_[T])) - 0.341*this->prm_[Cao])) / (std::exp((2.*(v_new*this->prm_[F]))/(this->prm_[R]*this->prm_[T])) - 1.)) * (this->var_[d]*this->var_[f1]*this->var_[f2]*this->var_[fCa]));
    
    this->prm_[alpha_K1] = 3.91 / (1. + std::exp(0.5942*((1000.*v_new - 1000.*this->prm_[E_K]) - 200.)));
    this->prm_[beta_K1] = (-1.509*std::exp(0.0002*((1000.*v_new - 1000.*this->prm_[E_K])+100.)) + std::exp(0.5886*((1000.*v_new - 1000.*this->prm_[E_K]) - 10.))) / (1.+std::exp(0.4547*(1000.*v_new - 1000.*this->prm_[E_K])));
    this->prm_[XK1_inf] = this->prm_[alpha_K1] / (this->prm_[alpha_K1]+this->prm_[beta_K1]);
    this->cur_[IK1] =  (1. - this->block_coeff_[IK1]) * (this->prm_[g_K1]*(this->prm_[XK1_inf] * ((v_new - this->prm_[E_K])*std::pow((this->prm_[Ko]/5.4), 0.5))));
    this->cur_[If] =  (1. - this->block_coeff_[If]) * (this->prm_[g_f]*(this->var_[Xf] * (v_new - this->prm_[E_f])));
    this->cur_[IKr] = (1. - this->block_coeff_[IKr]) * (this->prm_[g_Kr]*((v_new - this->prm_[E_K]) * (this->var_[Xr1]*(this->var_[Xr2]*std::pow((this->prm_[Ko]/5.4), 0.5)))));
    this->prm_[E_Ks] = ((this->prm_[R]*this->prm_[T])/this->prm_[F]) * std::log((this->prm_[Ko] + this->prm_[PkNa]*this->prm_[Nao]) / (this->prm_[Ki] + this->prm_[PkNa]*this->var_[Nai]));
    
    this->prm_[E_Ca] =  ((0.5*(this->prm_[R]*this->prm_[T]))/this->prm_[F]) * std::log(this->prm_[Cao]/this->var_[Cai]);
    this->cur_[IKs] = (1. - this->block_coeff_[IKs]) * (this->prm_[g_Ks]*((v_new - this->prm_[E_Ks])*( std::pow(this->var_[Xs], 2.) * (1. + 0.6/(1+std::pow((3.8*1.e-05)/this->var_[Cai], 1.4))))));
    this->cur_[Ito] = (1. - this->block_coeff_[Ito]) * (this->prm_[g_to]*((v_new - this->prm_[E_K])*(this->var_[q]*this->var_[r])));
    this->cur_[IPCa] = (1. - this->block_coeff_[IPCa]) * ((this->prm_[g_PCa]*this->var_[Cai]) / (this->var_[Cai]+this->prm_[KPCa]));
    this->cur_[ICab] = (1. - this->block_coeff_[ICab]) * (this->prm_[g_b_Ca] * (v_new - this->prm_[E_Ca]));
        
    this->cur_[PciVtrCur::Iion] = this->cur_[IK1] + this->cur_[Ito] + this->cur_[IKr] + this->cur_[IKs] + this->cur_[ICaL] + 
                                  this->cur_[INaK] + this->cur_[INa] + this->cur_[INaCa] + this->cur_[IPCa] + this->cur_[If] + 
                                  this->cur_[INab] + this->cur_[ICab];

    this->var_[dvdt] =  - (this->cur_[PciVtrCur::Iion] - stim_current/this->prm_[Cm]);

    this->prm_[Cai_bufc] = 1. / (1. + (this->prm_[Buf_C]*this->prm_[Kbuf_C]) / std::pow(this->var_[Cai]+this->prm_[Kbuf_C], 2.));
    this->cur_[Irel] = (1. - this->block_coeff_[Irel]) * ((this->prm_[c_rel] + (this->prm_[a_rel]*std::pow(this->var_[Ca_SR], 2.)) / (std::pow(this->prm_[b_rel], 2.) + std::pow(this->var_[Ca_SR], 2.))) * (this->var_[d]*(this->var_[g]*0.0411)));
    this->cur_[Iup] = (1. - this->block_coeff_[Iup]) * (this->prm_[VmaxUp] / (1. + std::pow(this->prm_[Kup], 2.) / std::pow(this->var_[Cai], 2.)));
    this->cur_[Ileak] = (1. - this->block_coeff_[Ileak]) * ((this->var_[Ca_SR] - this->var_[Cai])*this->prm_[V_leak]);
    
    this->var_[Cai] =  this->var_[Cai] + dt*this->prm_[Cai_bufc]*((this->cur_[Ileak] - this->cur_[Iup] + this->cur_[Irel]) - (((this->cur_[ICaL] + this->cur_[ICab] + this->cur_[IPCa]) -  2.*this->cur_[INaCa])*this->prm_[Cm])/(2.*(this->prm_[Vc]*(this->prm_[F]*1.e-18))));
    this->prm_[Ca_SR_bufSR] = 1. / (1. + (this->prm_[Buf_SR]*this->prm_[Kbuf_SR]) / std::pow(this->var_[Ca_SR]+this->prm_[Kbuf_SR], 2.));
    this->var_[Ca_SR] = this->var_[Ca_SR] + dt*((this->prm_[Ca_SR_bufSR]*this->prm_[Vc]) / this->prm_[V_SR]) * (this->cur_[Iup] - (this->cur_[Irel]+this->cur_[Ileak]));


}


std::string PaciVentri::PrintVariables() const 
{
    using namespace PciVtrVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v: " << this->var_[v] << "\n";
    oss << "dvdt: " << this->var_[dvdt] << "\n";
    oss << "Nai: " << this->var_[Nai] << "\n";
    oss << "Cai: " << this->var_[Cai] << "\n";
    oss << "m: " << this->var_[m] << "\n";
    oss << "h: " << this->var_[h] << "\n";
    oss << "j: " << this->var_[j] << "\n";
    oss << "d: " << this->var_[d] << "\n";
    oss << "f1: " << this->var_[f1] << "\n";
    oss << "f2: " << this->var_[f2] << "\n";
    oss << "fCa: " << this->var_[fCa] << "\n";
    oss << "Xr1: " << this->var_[Xr1] << "\n";
    oss << "Xr2: " << this->var_[Xr2] << "\n";
    oss << "Xs: " << this->var_[Xs] << "\n";
    oss << "Xf: " << this->var_[Xf] << "\n";
    oss << "q: " << this->var_[q] << "\n";
    oss << "r: " << this->var_[r] << "\n";
    oss << "Ca_SR: " << this->var_[Ca_SR] << "\n";
    oss << "g: " << this->var_[g];
    return oss.str();

}


std::string PaciVentri::PrintParameters() const 
{
    using namespace PciVtrPrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "Cm: " << this->prm_[Cm] << "\n";
    oss << "R: " << this->prm_[R] << "\n";
    oss << "T: " << this->prm_[T] << "\n";
    oss << "F: " << this->prm_[F] << "\n";
    oss << "Nao: " << this->prm_[Nao] << "\n";
    oss << "Cao: " << this->prm_[Cao] << "\n";
    oss << "Ki: " << this->prm_[Ki] << "\n";
    oss << "Ko: " << this->prm_[Ko] << "\n";
    oss << "PkNa: " << this->prm_[PkNa] << "\n";
    oss << "g_Na: " << this->prm_[g_Na] << "\n";
    oss << "g_CaL: " << this->prm_[g_CaL] << "\n";
    oss << "tau_fCa: " << this->prm_[tau_fCa] << "\n";
    oss << "g_Kr: " << this->prm_[g_Kr] << "\n";
    oss << "L0: " << this->prm_[L0] << "\n";
    oss << "Q: " << this->prm_[Q] << "\n";
    oss << "g_Ks: " << this->prm_[g_Ks] << "\n";
    oss << "g_K1: " << this->prm_[g_K1] << "\n";
    oss << "g_f: " << this->prm_[g_f] << "\n";
    oss << "E_f: " << this->prm_[E_f] << "\n";
    oss << "g_b_Na: " << this->prm_[g_b_Na] << "\n";
    oss << "g_b_Ca: " << this->prm_[g_b_Ca] << "\n";
    oss << "Km_K: " << this->prm_[Km_K] << "\n";
    oss << "Km_Na: " << this->prm_[Km_Na] << "\n";
    oss << "PNaK: " << this->prm_[PNaK] << "\n";
    oss << "kNaCa: " << this->prm_[kNaCa] << "\n";
    oss << "alpha: " << this->prm_[alpha] << "\n";
    oss << "gamma: " << this->prm_[gamma] << "\n";
    oss << "Ksat: " << this->prm_[Ksat] << "\n";
    oss << "KmCa: " << this->prm_[KmCa] << "\n";
    oss << "KmNai: " << this->prm_[KmNai] << "\n";
    oss << "g_PCa: " << this->prm_[g_PCa] << "\n";
    oss << "KPCa: " << this->prm_[KPCa] << "\n";
    oss << "g_to: " << this->prm_[g_to] << "\n";
    oss << "Vc: " << this->prm_[Vc] << "\n";
    oss << "V_SR: " << this->prm_[V_SR] << "\n";
    oss << "a_rel: " << this->prm_[a_rel] << "\n";
    oss << "b_rel: " << this->prm_[b_rel] << "\n";
    oss << "c_rel: " << this->prm_[c_rel] << "\n";
    oss << "tau_g: " << this->prm_[tau_g] << "\n";
    oss << "Kup: " << this->prm_[Kup] << "\n";
    oss << "Buf_C: " << this->prm_[Buf_C] << "\n";
    oss << "Buf_SR: " << this->prm_[Buf_SR] << "\n";
    oss << "Kbuf_C: " << this->prm_[Kbuf_C] << "\n";
    oss << "Kbuf_SR: " << this->prm_[Kbuf_SR] << "\n";
    oss << "VmaxUp: " << this->prm_[VmaxUp] << "\n";
    oss << "V_leak: " << this->prm_[V_leak] << "\n";
    oss << "E_K: " << this->prm_[E_K] << "\n";
    oss << "constf2: " << this->prm_[constf2] << "\n";
    oss << "V_half: " << this->prm_[V_half] << "\n";
    oss << "m_inf: " << this->prm_[m_inf] << "\n";
    oss << "h_inf: " << this->prm_[h_inf] << "\n";
    oss << "j_inf: " << this->prm_[j_inf] << "\n";
    oss << "d_inf: " << this->prm_[d_inf] << "\n";
    oss << "f1_inf: " << this->prm_[f1_inf] << "\n";
    oss << "f2_inf: " << this->prm_[f2_inf] << "\n";
    oss << "alpha_fCa: " << this->prm_[alpha_fCa] << "\n";
    oss << "Xr1_inf: " << this->prm_[Xr1_inf] << "\n";
    oss << "Xr2_inf: " << this->prm_[Xr2_inf] << "\n";
    oss << "Xs_inf: " << this->prm_[Xs_inf] << "\n";
    oss << "Xf_inf: " << this->prm_[Xf_inf] << "\n";
    oss << "q_inf: " << this->prm_[q_inf] << "\n";
    oss << "r_inf: " << this->prm_[r_inf] << "\n";
    oss << "g_inf: " << this->prm_[g_inf] << "\n";
    oss << "alpha_m: " << this->prm_[alpha_m] << "\n";
    oss << "alpha_h: " << this->prm_[alpha_h] << "\n";
    oss << "alpha_j: " << this->prm_[alpha_j] << "\n";
    oss << "alpha_d: " << this->prm_[alpha_d] << "\n";
    oss << "constf1: " << this->prm_[constf1] << "\n";
    oss << "tau_f2: " << this->prm_[tau_f2] << "\n";
    oss << "beta_fCa: " << this->prm_[beta_fCa] << "\n";
    oss << "alpha_Xr1: " << this->prm_[alpha_Xr1] << "\n";
    oss << "alpha_Xr2: " << this->prm_[alpha_Xr2] << "\n";
    oss << "alpha_Xs: " << this->prm_[alpha_Xs] << "\n";
    oss << "E_Na: " << this->prm_[E_Na] << "\n";
    oss << "tau_Xf: " << this->prm_[tau_Xf] << "\n";
    oss << "tau_q: " << this->prm_[tau_q] << "\n";
    oss << "tau_r: " << this->prm_[tau_r] << "\n";
    oss << "const2: " << this->prm_[const2] << "\n";
    oss << "beta_m: " << this->prm_[beta_m] << "\n";
    oss << "beta_h: " << this->prm_[beta_h] << "\n";
    oss << "beta_j: " << this->prm_[beta_j] << "\n";
    oss << "beta_d: " << this->prm_[beta_d] << "\n";
    oss << "tau_f1: " << this->prm_[tau_f1] << "\n";
    oss << "gamma_fCa: " << this->prm_[gamma_fCa] << "\n";
    oss << "beta_Xr1: " << this->prm_[beta_Xr1] << "\n";
    oss << "beta_Xr2: " << this->prm_[beta_Xr2] << "\n";
    oss << "beta_Xs: " << this->prm_[beta_Xs] << "\n";
    oss << "E_Ks: " << this->prm_[E_Ks] << "\n";
    oss << "tau_m: " << this->prm_[tau_m] << "\n";
    oss << "tau_h: " << this->prm_[tau_h] << "\n";
    oss << "tau_j: " << this->prm_[tau_j] << "\n";
    oss << "gamma_d: " << this->prm_[gamma_d] << "\n";
    oss << "fCa_inf: " << this->prm_[fCa_inf] << "\n";
    oss << "tau_Xr1: " << this->prm_[tau_Xr1] << "\n";
    oss << "tau_Xr2: " << this->prm_[tau_Xr2] << "\n";
    oss << "tau_Xs: " << this->prm_[tau_Xs] << "\n";
    oss << "E_Ca: " << this->prm_[E_Ca] << "\n";
    oss << "tau_d: " << this->prm_[tau_d] << "\n";
    oss << "constfCa: " << this->prm_[constfCa] << "\n";
    oss << "alpha_K1: " << this->prm_[alpha_K1] << "\n";
    oss << "beta_K1: " << this->prm_[beta_K1] << "\n";
    oss << "XK1_inf: " << this->prm_[XK1_inf] << "\n";
    oss << "Cai_bufc: " << this->prm_[Cai_bufc] << "\n";
    oss << "Ca_SR_bufSR: " << this->prm_[Ca_SR_bufSR];
    return oss.str();

}


std::string PaciVentri::PrintCurrents() const
{
    using namespace PciVtrCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa: " << this->cur_[INa] << "\n";
    oss << "INaK: " << this->cur_[INaK] << "\n";
    oss << "INaCa: " << this->cur_[INaCa] << "\n";
    oss << "INab: " << this->cur_[INab] << "\n";
    oss << "ICaL: " << this->cur_[ICaL] << "\n";
    oss << "IK1: " << this->cur_[IK1] << "\n";
    oss << "If: " << this->cur_[If] << "\n";
    oss << "IKr: " << this->cur_[IKr] << "\n";
    oss << "IKs: " << this->cur_[IKs] << "\n";
    oss << "Ito: " << this->cur_[Ito] << "\n";
    oss << "IPCa: " << this->cur_[IPCa] << "\n";
    oss << "ICab: " << this->cur_[ICab] << "\n";
    oss << "Irel: " << this->cur_[Irel] << "\n";
    oss << "Iup: " << this->cur_[Iup] << "\n";
    oss << "Ileak: " << this->cur_[Ileak] << "\n";
    oss << "Iion: " << this->cur_[PciVtrCur::Iion];
    return oss.str();

}


std::string PaciVentri::PrintBlockCoeffs() const
{
    using namespace PciVtrCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa: " << this->block_coeff_[INa] << "\n";
    oss << "INaK: " << this->block_coeff_[INaK] << "\n";
    oss << "INaCa: " << this->block_coeff_[INaCa] << "\n";
    oss << "INab: " << this->block_coeff_[INab] << "\n";
    oss << "ICaL: " << this->block_coeff_[ICaL] << "\n";
    oss << "IK1: " << this->block_coeff_[IK1] << "\n";
    oss << "If: " << this->block_coeff_[If] << "\n";
    oss << "IKr: " << this->block_coeff_[IKr] << "\n";
    oss << "IKs: " << this->block_coeff_[IKs] << "\n";
    oss << "Ito: " << this->block_coeff_[Ito] << "\n";
    oss << "IPCa: " << this->block_coeff_[IPCa] << "\n";
    oss << "ICab: " << this->block_coeff_[ICab] << "\n";
    oss << "Irel: " << this->block_coeff_[Irel] << "\n";
    oss << "Iup: " << this->block_coeff_[Iup] << "\n";
    oss << "Ileak: " << this->block_coeff_[Ileak];
    return oss.str();

}

} // End of namespace ELECTRA
