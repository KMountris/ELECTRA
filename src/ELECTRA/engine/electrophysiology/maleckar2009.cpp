/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/electrophysiology/maleckar2009.hpp"


namespace ELECTRA {


void Maleckar2009::SetDataMapping() 
{
    using namespace MlcrVar;
    using namespace MlcrPrm;
    using namespace MlcrCur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["Na_c"] = static_cast<std::size_t>(Na_c);
    this->mapped_data_["Na_i"] = static_cast<std::size_t>(Na_i);
    this->mapped_data_["m"] = static_cast<std::size_t>(m);
    this->mapped_data_["h1"] = static_cast<std::size_t>(h1);
    this->mapped_data_["h2"] = static_cast<std::size_t>(h2);
    this->mapped_data_["Ca_d"] = static_cast<std::size_t>(Ca_d);
    this->mapped_data_["d_L"] = static_cast<std::size_t>(d_L);
    this->mapped_data_["f_L1"] = static_cast<std::size_t>(f_L1);
    this->mapped_data_["f_L2"] = static_cast<std::size_t>(f_L2);
    this->mapped_data_["K_c"] = static_cast<std::size_t>(K_c);
    this->mapped_data_["K_i"] = static_cast<std::size_t>(K_i);
    this->mapped_data_["r"] = static_cast<std::size_t>(r);
    this->mapped_data_["s"] = static_cast<std::size_t>(s);
    this->mapped_data_["a_ur"] = static_cast<std::size_t>(a_ur);
    this->mapped_data_["i_ur"] = static_cast<std::size_t>(i_ur);
    this->mapped_data_["n"] = static_cast<std::size_t>(n);
    this->mapped_data_["pa"] = static_cast<std::size_t>(pa);
    this->mapped_data_["Ca_c"] = static_cast<std::size_t>(Ca_c);
    this->mapped_data_["Ca_i"] = static_cast<std::size_t>(Ca_i);
    this->mapped_data_["O_C"] = static_cast<std::size_t>(O_C);
    this->mapped_data_["O_TC"] = static_cast<std::size_t>(O_TC);
    this->mapped_data_["O_TMgC"] = static_cast<std::size_t>(O_TMgC);
    this->mapped_data_["O_TMgMg"] = static_cast<std::size_t>(O_TMgMg);
    this->mapped_data_["O"] = static_cast<std::size_t>(O);
    this->mapped_data_["Ca_rel"] = static_cast<std::size_t>(Ca_rel);
    this->mapped_data_["Ca_up"] = static_cast<std::size_t>(Ca_up);
    this->mapped_data_["O_Calse"] = static_cast<std::size_t>(O_Calse);
    this->mapped_data_["F1"] = static_cast<std::size_t>(F1);
    this->mapped_data_["F2"] = static_cast<std::size_t>(F2);
   

    // Set parameters mapping.
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["F"] = static_cast<std::size_t>(F);
    this->mapped_data_["Cm"] = static_cast<std::size_t>(Cm);
    this->mapped_data_["P_Na"] = static_cast<std::size_t>(P_Na);
    this->mapped_data_["g_Ca_L"] = static_cast<std::size_t>(g_Ca_L);
    this->mapped_data_["E_Ca_app"] = static_cast<std::size_t>(E_Ca_app);
    this->mapped_data_["k_Ca"] = static_cast<std::size_t>(k_Ca);
    this->mapped_data_["g_t"] = static_cast<std::size_t>(g_t);
    this->mapped_data_["g_kur"] = static_cast<std::size_t>(g_kur);
    this->mapped_data_["g_K1"] = static_cast<std::size_t>(g_K1);
    this->mapped_data_["g_Ks"] = static_cast<std::size_t>(g_Ks);
    this->mapped_data_["g_Kr"] = static_cast<std::size_t>(g_Kr);
    this->mapped_data_["g_B_Na"] = static_cast<std::size_t>(g_B_Na);
    this->mapped_data_["g_B_Ca"] = static_cast<std::size_t>(g_B_Ca);
    this->mapped_data_["K_NaK_K"] = static_cast<std::size_t>(K_NaK_K);
    this->mapped_data_["i_NaK_max"] = static_cast<std::size_t>(i_NaK_max);
    this->mapped_data_["pow_K_NaK_Na_15"] = static_cast<std::size_t>(pow_K_NaK_Na_15);
    this->mapped_data_["i_CaP_max"] = static_cast<std::size_t>(i_CaP_max);
    this->mapped_data_["k_CaP"] = static_cast<std::size_t>(k_CaP);
    this->mapped_data_["K_NaCa"] = static_cast<std::size_t>(K_NaCa);
    this->mapped_data_["d_NaCa"] = static_cast<std::size_t>(d_NaCa);
    this->mapped_data_["gamma_Na"] = static_cast<std::size_t>(gamma_Na);
    this->mapped_data_["ACh"] = static_cast<std::size_t>(ACh);
    this->mapped_data_["phi_Na_en"] = static_cast<std::size_t>(phi_Na_en);
    this->mapped_data_["Vol_i"] = static_cast<std::size_t>(Vol_i);
    this->mapped_data_["Vol_d"] = static_cast<std::size_t>(Vol_d);
    this->mapped_data_["tau_di"] = static_cast<std::size_t>(tau_di);
    this->mapped_data_["Mg_i"] = static_cast<std::size_t>(Mg_i);
    this->mapped_data_["Vol_c"] = static_cast<std::size_t>(Vol_c);
    this->mapped_data_["tau_Na"] = static_cast<std::size_t>(tau_Na);
    this->mapped_data_["tau_K"] = static_cast<std::size_t>(tau_K);
    this->mapped_data_["tau_Ca"] = static_cast<std::size_t>(tau_Ca);
    this->mapped_data_["Na_b"] = static_cast<std::size_t>(Na_b);
    this->mapped_data_["Ca_b"] = static_cast<std::size_t>(Ca_b);
    this->mapped_data_["K_b"] = static_cast<std::size_t>(K_b);
    this->mapped_data_["I_up_max"] = static_cast<std::size_t>(I_up_max);
    this->mapped_data_["k_cyca"] = static_cast<std::size_t>(k_cyca);
    this->mapped_data_["k_srca"] = static_cast<std::size_t>(k_srca);
    this->mapped_data_["k_xcs"] = static_cast<std::size_t>(k_xcs);
    this->mapped_data_["alpha_rel"] = static_cast<std::size_t>(alpha_rel);
    this->mapped_data_["Vol_up"] = static_cast<std::size_t>(Vol_up);
    this->mapped_data_["Vol_rel"] = static_cast<std::size_t>(Vol_rel);
    this->mapped_data_["r_recov"] = static_cast<std::size_t>(r_recov);
    this->mapped_data_["tau_tr"] = static_cast<std::size_t>(tau_tr);
    this->mapped_data_["k_rel_i"] = static_cast<std::size_t>(k_rel_i);
    this->mapped_data_["k_rel_d"] = static_cast<std::size_t>(k_rel_d);
  
    // Set currents mapping.
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["ICaL"] = static_cast<std::size_t>(ICaL);
    this->mapped_data_["Ito"] = static_cast<std::size_t>(Ito);
    this->mapped_data_["IKur"] = static_cast<std::size_t>(IKur);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["ICaP"] = static_cast<std::size_t>(ICaP);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["IKACh"] = static_cast<std::size_t>(IKACh);

}


Maleckar2009::Maleckar2009()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Maleckar2009;
    this->dt_stable_ = 0.005;
    this->var_.resize(31, 0.);
    this->prm_.resize(47, 0.);
    this->cur_.resize(14, 0.);
    this->block_coeff_.resize(13, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


Maleckar2009::~Maleckar2009()
{}


void Maleckar2009::Initialize(CellType cell_type)
{
    using namespace MlcrVar;
    using namespace MlcrPrm;

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(31, 0.);
    this->prm_.clear();           this->prm_.resize(47, 0.);
    this->cur_.clear();           this->cur_.resize(14, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(13, 0.);

    // Check cell type.
    if (cell_type != CellType::atrial) {
            std::string error_str = "Could not initialize the Maleckar2009 atrial ap model. Expected cell type: ELECTRA::CellType::atrial";
            throw std::invalid_argument(Logger::Error(error_str));
    }

    // Set the model variables.
    this->var_[v] = -74.031982;
    this->var_[dvdt] = 0.;
    this->var_[Na_c] = 130.022096;
    this->var_[Na_i] = 8.516766;
    this->var_[m] = 0.003289;
    this->var_[h1] = 0.877202;
    this->var_[h2] = 0.873881;
    this->var_[Ca_d] = 7.1e-5;
    this->var_[d_L] = 0.000014;
    this->var_[f_L1] = 0.998597;
    this->var_[f_L2] = 0.998586;
    this->var_[K_c] = 5.560224;
    this->var_[K_i] = 129.485991;
    this->var_[r] = 0.001089;
    this->var_[s] = 0.948597;
    this->var_[a_ur] = 0.000367;
    this->var_[i_ur] = 0.96729;
    this->var_[n] = 0.004374;
    this->var_[pa] = 0.000053;
    this->var_[Ca_c] = 1.815768;
    this->var_[Ca_i] = 6.5e-5;
    this->var_[O_C] = 0.026766;
    this->var_[O_TC] = 0.012922;
    this->var_[O_TMgC] = 0.190369;
    this->var_[O_TMgMg] = 0.714463;
    this->var_[O] = 1.38222;
    this->var_[Ca_rel] = 0.632613;
    this->var_[Ca_up] = 0.649195;
    this->var_[O_Calse] = 0.431547;
    this->var_[F1] = 0.470055;
    this->var_[F2] = 0.002814;
    
    // Set the univarsal model parameters.
    this->prm_[R] = 8314.;
    this->prm_[T] = 306.15;
    this->prm_[F] = 96487;
    this->prm_[Cm] = 50.;
    this->prm_[P_Na] = 0.0018;
    this->prm_[g_Ca_L] = 6.75;
    this->prm_[E_Ca_app] = 60.;
    this->prm_[k_Ca] = 0.025;
    this->prm_[g_t] = 8.25;
    this->prm_[g_kur] = 2.25;
    this->prm_[g_K1] = 3.1;
    this->prm_[g_Ks] = 1.;
    this->prm_[g_Kr] = 0.5;
    this->prm_[g_B_Na] = 0.060599;
    this->prm_[g_B_Ca] = 0.078681;
    this->prm_[K_NaK_K] = 1.;
    this->prm_[i_NaK_max] = 68.55;
    this->prm_[pow_K_NaK_Na_15] = 36.4829;
    this->prm_[i_CaP_max] = 4.;
    this->prm_[k_CaP] = 0.0002;
    this->prm_[K_NaCa] = 0.0374842;
    this->prm_[d_NaCa] = 0.0003;
    this->prm_[gamma_Na] = 0.45;
    this->prm_[ACh] = 1e-24;
    this->prm_[phi_Na_en] = 0.;
    this->prm_[Vol_i] = 0.005884;
    this->prm_[Vol_d] = 0.00011768;
    this->prm_[tau_di] = 0.01;
    this->prm_[Mg_i] = 2.5;
    this->prm_[Vol_c] = 0.000800224;
    this->prm_[tau_Na] = 14.3;
    this->prm_[tau_K] = 10.;
    this->prm_[tau_Ca] = 24.7;
    this->prm_[Na_b] = 130.;
    this->prm_[Ca_b] = 1.8;
    this->prm_[K_b] = 5.4;
    this->prm_[I_up_max] = 2800.;
    this->prm_[k_cyca] = 0.0003;
    this->prm_[k_srca] = 0.5;
    this->prm_[k_xcs] = 0.4;
    this->prm_[alpha_rel] = 200000.;
    this->prm_[Vol_up] = 0.0003969;
    this->prm_[Vol_rel] = 0.0000441;
    this->prm_[r_recov] = 0.815;
    this->prm_[tau_tr] = 0.01;
    this->prm_[k_rel_i] = 0.0003;
    this->prm_[k_rel_d] = 0.003;

}


void Maleckar2009::Compute(double v_new, double dt, double stim_current)
{
    using namespace MlcrVar;
    using namespace MlcrPrm;
    using namespace MlcrCur;   

    dt /= 1000.;
    
    double dO_TMgMg = 2000.00*this->prm_[Mg_i]*((1.00000 - this->var_[O_TMgC]) - this->var_[O_TMgMg]) -  666.000*this->var_[O_TMgMg];
    this->var_[O_TMgMg] = this->var_[O_TMgMg] + dt*dO_TMgMg;
    
    double tau_r =  0.00350000*std::exp(((  - v_new*v_new)/30.0000)/30.0000)+0.00150000;
    double r_infinity = 1.00000/(1.00000+std::exp((v_new - 1.00000)/ - 11.0000));
    this->var_[r] = r_infinity - (r_infinity - this->var_[r]) * std::exp(-dt/tau_r);
    
    double a_ur_infinity = 1.00000/(1.00000+std::exp( - (v_new+6.00000)/8.60000));
    double tau_a_ur = 0.00900000/(1.00000+std::exp((v_new+5.00000)/12.0000))+0.000500000;
    this->var_[a_ur] = a_ur_infinity - (a_ur_infinity - this->var_[a_ur]) * std::exp(-dt/tau_a_ur);
    
    double i_ur_infinity = 1.00000/(1.00000+std::exp((v_new+7.50000)/10.0000));
    double tau_i_ur = 0.590000/(1.00000+std::exp((v_new+60.0000)/10.0000))+3.05000;
    this->var_[i_ur] = i_ur_infinity - (i_ur_infinity - this->var_[i_ur]) * std::exp(-dt/tau_i_ur);
    
    double m_infinity = 1.00000/(1.00000+std::exp((v_new+27.1200)/ - 8.21000));
    double m_factor = (v_new+25.5700)/28.8000;
    double tau_m =  4.20000e-05*std::exp(  - m_factor*m_factor)+2.40000e-05;
    this->var_[m] = m_infinity - (m_infinity - this->var_[m]) * std::exp(-dt/tau_m);
    
    double h_infinity = 1.00000/(1.00000+std::exp((v_new+63.6000)/5.30000));
    double h_factor = 1.00000/(1.00000+std::exp((v_new+35.1000)/3.20000));
    double tau_h1 =  0.0300000*h_factor+0.000300000;
    this->var_[h1] = h_infinity - (h_infinity - this->var_[h1]) * std::exp(-dt/tau_h1);    
    
    double tau_h2 =  0.120000*h_factor+0.00300000;
    this->var_[h2] = h_infinity - (h_infinity - this->var_[h2]) * std::exp(-dt/tau_h2);
    
    double d_L_infinity = 1.00000/(1.00000+std::exp((v_new+9.00000)/ - 5.80000));
    double d_L_factor = (v_new+35.0000)/30.0000;
    double tau_d_L =  0.00270000*std::exp(  - d_L_factor*d_L_factor)+0.00200000;
    this->var_[d_L] = d_L_infinity - (d_L_infinity - this->var_[d_L]) * std::exp(-dt/tau_d_L);
    
    double f_L_infinity = 1.00000/(1.00000+std::exp((v_new+27.4000)/7.10000));
    double f_L_factor = v_new+40.0000;
    double tau_f_L1 =  0.161000*std::exp(((  - f_L_factor*f_L_factor)/14.4000)/14.4000)+0.0100000;
    this->var_[f_L1] = f_L_infinity - (f_L_infinity - this->var_[f_L1]) * std::exp(-dt/tau_f_L1);
    
    double tau_f_L2 =  1.33230*std::exp(((  - f_L_factor*f_L_factor)/14.2000)/14.2000)+0.0626000;
    this->var_[f_L2] = f_L_infinity - (f_L_infinity - this->var_[f_L2]) * std::exp(-dt/tau_f_L2);

    double s_factor = (v_new+52.4500)/15.8827;
    double tau_s =  0.0256350*std::exp(  - s_factor*s_factor)+0.0141400;
    double s_infinity = 1.00000/(1.00000+std::exp((v_new+40.5000)/11.5000));
    this->var_[s] = s_infinity - (s_infinity - this->var_[s]) * std::exp(-dt/tau_s);

    double n_factor = (v_new - 20.0000)/20.0000;
    double tau_n = 0.700000+ 0.400000*std::exp(  - n_factor*n_factor);
    double n_infinity = 1.00000/(1.00000+std::exp((v_new - 19.9000)/ - 12.7000));
    this->var_[n] = n_infinity - (n_infinity - this->var_[n]) * std::exp(-dt/tau_n);
    
    double pa_factor = (v_new+20.1376)/22.1996;
    double tau_pa = 0.0311800+ 0.217180*std::exp(  - pa_factor*pa_factor);
    double p_a_infinity = 1.00000/(1.00000+std::exp((v_new+15.0000)/ - 6.00000));
    this->var_[pa] = p_a_infinity - (p_a_infinity - this->var_[pa]) * std::exp(-dt/tau_pa);
    
    double r_Ca_d_term = this->var_[Ca_d]/(this->var_[Ca_d]+this->prm_[k_rel_d]);
    double r_Ca_d_factor =  r_Ca_d_term*r_Ca_d_term*r_Ca_d_term*r_Ca_d_term;
    double r_Ca_i_term = this->var_[Ca_i]/(this->var_[Ca_i]+this->prm_[k_rel_i]);
    double r_Ca_i_factor =  r_Ca_i_term*r_Ca_i_term*r_Ca_i_term*r_Ca_i_term;
    double r_act =  203.800*(r_Ca_i_factor+r_Ca_d_factor);
    double dF1 = this->prm_[r_recov]*((1.00000 - this->var_[F1]) - this->var_[F2]) -  r_act*this->var_[F1];
    this->var_[F1] = this->var_[F1] + dt*dF1;
    
    double r_inact = 33.9600+ 339.600*r_Ca_i_factor;
    double dF2 = r_act*this->var_[F1] -  r_inact*this->var_[F2];
    this->var_[F2] = this->var_[F2] + dt*dF2;
    
    double E_K = (( this->prm_[R]*this->prm_[T])/this->prm_[F])*std::log(this->var_[K_c]/this->var_[K_i]);
    this->cur_[Ito] = (1.0-this->block_coeff_[Ito]) * (this->prm_[g_t]*this->var_[r]*this->var_[s]*(v_new - E_K));
    this->cur_[IKur] = (1.0-this->block_coeff_[IKur]) * (this->prm_[g_kur]*this->var_[a_ur]*this->var_[i_ur]*(v_new - E_K));
    this->cur_[IK1] = (1.0-this->block_coeff_[IK1]) * ( (this->prm_[g_K1]*std::pow(this->var_[K_c], 0.445700)*(v_new - E_K))/(1.00000+std::exp(( 1.50000*((v_new - E_K)+3.60000)*this->prm_[F])/( this->prm_[R]*this->prm_[T]))) );
    
    double pip = 1.00000/(1.00000+std::exp((v_new+55.0000)/24.0000));
    this->cur_[IKr] = (1.0-this->block_coeff_[IKr]) * (this->prm_[g_Kr]*this->var_[pa]*pip*(v_new - E_K));
    this->cur_[IKs] = (1.0-this->block_coeff_[IKs]) * (this->prm_[g_Ks]*this->var_[n]*(v_new - E_K));
    double pow_Na_i_15 = std::pow(this->var_[Na_i], 1.50000);
    this->cur_[INaK] = (1.0-this->block_coeff_[INaK]) * ( ((((( this->prm_[i_NaK_max]*this->var_[K_c])/(this->var_[K_c]+this->prm_[K_NaK_K]))*pow_Na_i_15)/(pow_Na_i_15+this->prm_[pow_K_NaK_Na_15]))*(v_new+150.000))/(v_new+200.000) );
    
    double dK_i =  - (((this->cur_[Ito]+this->cur_[IKur]+this->cur_[IK1]+this->cur_[IKs]+this->cur_[IKr]) -  2.00000*this->cur_[INaK])+ stim_current*this->prm_[Cm])/( this->prm_[Vol_i]*this->prm_[F]);
    this->var_[K_i] = this->var_[K_i] + dt*dK_i;
    
    double dK_c = (this->prm_[K_b] - this->var_[K_c])/this->prm_[tau_K]+((this->cur_[Ito]+this->cur_[IKur]+this->cur_[IK1]+this->cur_[IKs]+this->cur_[IKr]) -  2.00000*this->cur_[INaK])/( this->prm_[Vol_c]*this->prm_[F]);
    this->var_[K_c] = this->var_[K_c] + dt*dK_c;
    
    double E_Na =  (( this->prm_[R]*this->prm_[T])/this->prm_[F])*std::log(this->var_[Na_c]/this->var_[Na_i]);
    this->cur_[INa] = (1.0-this->block_coeff_[INa]) * ( ((( this->prm_[P_Na]*this->var_[m]*this->var_[m]*this->var_[m]*( 0.900000*this->var_[h1]+ 0.100000*this->var_[h2])*this->var_[Na_c]*v_new*this->prm_[F]*this->prm_[F])/( this->prm_[R]*this->prm_[T]))*(std::exp(( (v_new - E_Na)*this->prm_[F])/( this->prm_[R]*this->prm_[T])) - 1.00000))/(std::exp(( v_new*this->prm_[F])/( this->prm_[R]*this->prm_[T])) - 1.00000) );
    this->cur_[INab] = (1.0-this->block_coeff_[INab]) * (this->prm_[g_B_Na]*(v_new - E_Na));
    this->cur_[INaCa] = (1.0-this->block_coeff_[INaCa]) * ( (this->prm_[K_NaCa]*( this->var_[Na_i]*this->var_[Na_i]*this->var_[Na_i]*this->var_[Ca_c]*std::exp(( this->prm_[F]*v_new*this->prm_[gamma_Na])/( this->prm_[R]*this->prm_[T])) -  this->var_[Na_c]*this->var_[Na_c]*this->var_[Na_c]*this->var_[Ca_i]*std::exp(( (this->prm_[gamma_Na] - 1.00000)*v_new*this->prm_[F])/( this->prm_[R]*this->prm_[T]))))/(1.00000+ this->prm_[d_NaCa]*( this->var_[Na_c]*this->var_[Na_c]*this->var_[Na_c]*this->var_[Ca_i]+ this->var_[Na_i]*this->var_[Na_i]*this->var_[Na_i]*this->var_[Ca_c])));
    double dNa_i =  - (this->cur_[INa]+this->cur_[INab]+ 3.00000*this->cur_[INaCa]+ 3.00000*this->cur_[INaK]+this->prm_[phi_Na_en])/( this->prm_[Vol_i]*this->prm_[F]);
    this->var_[Na_i] = this->var_[Na_i] + dt*dNa_i;
    
    double f_Ca = this->var_[Ca_d]/(this->var_[Ca_d]+this->prm_[k_Ca]);
    this->cur_[ICaL] = (1.0-this->block_coeff_[ICaL]) * (this->prm_[g_Ca_L]*this->var_[d_L]*( f_Ca*this->var_[f_L1] + (1.00000 - f_Ca)*this->var_[f_L2])*(v_new - this->prm_[E_Ca_app]) );
    double E_Ca =  (( this->prm_[R]*this->prm_[T])/( 2.00000*this->prm_[F]))*std::log(this->var_[Ca_c]/this->var_[Ca_i]);
    this->cur_[ICab] = (1.0-this->block_coeff_[ICab]) * (this->prm_[g_B_Ca]*(v_new - E_Ca));
    this->cur_[ICaP] = (1.0-this->block_coeff_[ICaP]) * ( ( this->prm_[i_CaP_max]*this->var_[Ca_i])/(this->var_[Ca_i]+this->prm_[k_CaP]) );
    double dCa_c = (this->prm_[Ca_b] - this->var_[Ca_c])/this->prm_[tau_Ca]+((this->cur_[ICaL]+this->cur_[ICab]+this->cur_[ICaP]) -  2.00000*this->cur_[INaCa])/( 2.00000*this->prm_[Vol_c]*this->prm_[F]);
    this->var_[Ca_c] = this->var_[Ca_c] + dt*dCa_c;
    
    double dNa_c = (this->prm_[Na_b] - this->var_[Na_c])/this->prm_[tau_Na]+(this->cur_[INa]+this->cur_[INab] + 3.00000*this->cur_[INaCa] + 3.00000*this->cur_[INaK]+this->prm_[phi_Na_en])/( this->prm_[Vol_c]*this->prm_[F]);
    this->var_[Na_c] = this->var_[Na_c] + dt*dNa_c;
    
    double i_di = ( (this->var_[Ca_d] - this->var_[Ca_i])*2.00000*this->prm_[Vol_d]*this->prm_[F])/this->prm_[tau_di];
    double dCa_d =  - (this->cur_[ICaL]+i_di)/( 2.00000*this->prm_[Vol_d]*this->prm_[F]);
    this->var_[Ca_d] = this->var_[Ca_d] + dt*dCa_d;
    
    this->cur_[IKACh] = (1.0-this->block_coeff_[IKACh]) * ((10.0000/(1.00000+( 9.13652*std::pow(1.00000, 0.477811))/std::pow(this->prm_[ACh], 0.477811)))*(0.0517000+0.451600/(1.00000+std::exp((v_new+59.5300)/17.1800)))*(v_new - E_K)*this->prm_[Cm] );
    this->cur_[MlcrCur::Iion] = (this->cur_[INa]+this->cur_[ICaL]+this->cur_[Ito]+this->cur_[IKur]+this->cur_[IK1]+this->cur_[IKr]+this->cur_[IKs]+this->cur_[INab]+this->cur_[ICab]+this->cur_[INaK]+this->cur_[ICaP]+this->cur_[INaCa]+this->cur_[IKACh])/this->prm_[Cm];
    this->var_[dvdt] = - (this->cur_[MlcrCur::Iion] - stim_current);
    
    double J_O_C =  200000.*this->var_[Ca_i]*(1.00000 - this->var_[O_C]) -  476.000*this->var_[O_C];
    double dO_C = J_O_C;
    this->var_[O_C] = this->var_[O_C] + dt*dO_C;
    
    double J_O_TC =  78400.0*this->var_[Ca_i]*(1.00000 - this->var_[O_TC]) -  392.000*this->var_[O_TC];
    double dO_TC = J_O_TC;
    this->var_[O_TC] = this->var_[O_TC] + dt*dO_TC;
    
    double J_O_TMgC =  200000.*this->var_[Ca_i]*((1.00000 - this->var_[O_TMgC]) - this->var_[O_TMgMg]) -  6.60000*this->var_[O_TMgC];
    double dO_TMgC = J_O_TMgC;
    this->var_[O_TMgC] = this->var_[O_TMgC] + dt*dO_TMgC;
    
    double J_O =  0.0800000*J_O_TC+ 0.160000*J_O_TMgC+ 0.0450000*J_O_C;
    double dO = J_O;
    this->var_[O] = this->var_[O] + dt*dO;
    
    double i_up = ( this->prm_[I_up_max]*(this->var_[Ca_i]/this->prm_[k_cyca] - ( this->prm_[k_xcs]*this->prm_[k_xcs]*this->var_[Ca_up])/this->prm_[k_srca]))/((this->var_[Ca_i]+this->prm_[k_cyca])/this->prm_[k_cyca]+( this->prm_[k_xcs]*(this->var_[Ca_up]+this->prm_[k_srca]))/this->prm_[k_srca]);
    double i_rel_f2 = this->var_[F2]/(this->var_[F2]+0.250000);
    double i_rel_factor =  i_rel_f2*i_rel_f2;
    double i_rel =  this->prm_[alpha_rel]*i_rel_factor*(this->var_[Ca_rel] - this->var_[Ca_i]);
    double dCa_i =  - ((this->cur_[ICab]+this->cur_[ICaP]+i_up) - (i_di+i_rel+ 2.00000*this->cur_[INaCa]))/( 2.00000*this->prm_[Vol_i]*this->prm_[F]) -  1.00000*J_O;
    this->var_[Ca_i] = this->var_[Ca_i] + dt*dCa_i;
    
    double i_tr = ( (this->var_[Ca_up] - this->var_[Ca_rel])*2.00000*this->prm_[Vol_rel]*this->prm_[F])/this->prm_[tau_tr];
    double dCa_up = (i_up - i_tr)/( 2.00000*this->prm_[Vol_up]*this->prm_[F]);
    this->var_[Ca_up] = this->var_[Ca_up] + dt*dCa_up;
    
    double J_O_Calse =  480.000*this->var_[Ca_rel]*(1.00000 - this->var_[O_Calse]) -  400.000*this->var_[O_Calse];
    double dO_Calse = J_O_Calse;
    this->var_[O_Calse] = this->var_[O_Calse] + dt*dO_Calse;
    
    double dCa_rel = (i_tr - i_rel)/( 2.00000*this->prm_[Vol_rel]*this->prm_[F]) -  31.0000*J_O_Calse;
    this->var_[Ca_rel] = this->var_[Ca_rel] + dt*dCa_rel; 
   
}



std::string Maleckar2009::PrintVariables() const 
{
    using namespace MlcrVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "Na_c = " << this->var_[Na_c] << "\n";
    oss << "Na_i = " << this->var_[Na_i] << "\n";
    oss << "m = " << this->var_[m] << "\n";
    oss << "h1 = " << this->var_[h1] << "\n";
    oss << "h2 = " << this->var_[h2] << "\n";
    oss << "Ca_d = " << this->var_[Ca_d] << "\n";
    oss << "d_L = " << this->var_[d_L] << "\n";
    oss << "f_L1 = " << this->var_[f_L1] << "\n";
    oss << "f_L2 = " << this->var_[f_L2] << "\n";
    oss << "K_c = " << this->var_[K_c] << "\n";
    oss << "K_i = " << this->var_[K_i] << "\n";
    oss << "r = " << this->var_[r] << "\n";
    oss << "s = " << this->var_[s] << "\n";
    oss << "a_ur = " << this->var_[a_ur] << "\n";
    oss << "i_ur = " << this->var_[i_ur] << "\n";
    oss << "n = " << this->var_[n] << "\n";
    oss << "pa = " << this->var_[pa] << "\n";
    oss << "Ca_c = " << this->var_[Ca_c] << "\n";
    oss << "Ca_i = " << this->var_[Ca_i] << "\n";
    oss << "O_C = " << this->var_[O_C] << "\n";
    oss << "O_TC = " << this->var_[O_TC] << "\n";
    oss << "O_TMgC = " << this->var_[O_TMgC] << "\n";
    oss << "O_TMgMg = " << this->var_[O_TMgMg] << "\n";
    oss << "O = " << this->var_[O] << "\n";
    oss << "Ca_rel = " << this->var_[Ca_rel] << "\n";
    oss << "Ca_up = " << this->var_[Ca_up] << "\n";
    oss << "O_Calse = " << this->var_[O_Calse] << "\n";
    oss << "F1 = " << this->var_[F1] << "\n";
    oss << "F2 = " << this->var_[F2];
    return oss.str();

}


std::string Maleckar2009::PrintParameters() const 
{
    using namespace MlcrPrm;    

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "R = " << this->prm_[R] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "F = " << this->prm_[F] << "\n";
    oss << "Cm = " << this->prm_[Cm] << "\n";
    oss << "P_Na = " << this->prm_[P_Na] << "\n";
    oss << "g_Ca_L = " << this->prm_[g_Ca_L] << "\n";
    oss << "E_Ca_app = " << this->prm_[E_Ca_app] << "\n";
    oss << "k_Ca = " << this->prm_[k_Ca] << "\n";
    oss << "g_t = " << this->prm_[g_t] << "\n";
    oss << "g_kur = " << this->prm_[g_kur] << "\n";
    oss << "g_K1 = " << this->prm_[g_K1] << "\n";
    oss << "g_Ks = " << this->prm_[g_Ks] << "\n";
    oss << "g_Kr = " << this->prm_[g_Kr] << "\n";
    oss << "g_B_Na = " << this->prm_[g_B_Na] << "\n";
    oss << "g_B_Ca = " << this->prm_[g_B_Ca] << "\n";
    oss << "K_NaK_K = " << this->prm_[K_NaK_K] << "\n";
    oss << "i_NaK_max = " << this->prm_[i_NaK_max] << "\n";
    oss << "pow_K_NaK_Na_15 = " << this->prm_[pow_K_NaK_Na_15] << "\n";
    oss << "i_CaP_max = " << this->prm_[i_CaP_max] << "\n";
    oss << "k_CaP = " << this->prm_[k_CaP] << "\n";
    oss << "K_NaCa = " << this->prm_[K_NaCa] << "\n";
    oss << "d_NaCa = " << this->prm_[d_NaCa] << "\n";
    oss << "gamma_Na = " << this->prm_[gamma_Na] << "\n";
    oss << "ACh = " << this->prm_[ACh] << "\n";
    oss << "phi_Na_en = " << this->prm_[phi_Na_en] << "\n";
    oss << "Vol_i = " << this->prm_[Vol_i] << "\n";
    oss << "Vol_d = " << this->prm_[Vol_d] << "\n";
    oss << "tau_di = " << this->prm_[tau_di] << "\n";
    oss << "Mg_i = " << this->prm_[Mg_i] << "\n";
    oss << "Vol_c = " << this->prm_[Vol_c] << "\n";
    oss << "tau_Na = " << this->prm_[tau_Na] << "\n";
    oss << "tau_K = " << this->prm_[tau_K] << "\n";
    oss << "tau_Ca = " << this->prm_[tau_Ca] << "\n";
    oss << "Na_b = " << this->prm_[Na_b] << "\n";
    oss << "Ca_b = " << this->prm_[Ca_b] << "\n";
    oss << "K_b = " << this->prm_[K_b] << "\n";
    oss << "I_up_max = " << this->prm_[I_up_max] << "\n";
    oss << "k_cyca = " << this->prm_[k_cyca] << "\n";
    oss << "k_srca = " << this->prm_[k_srca] << "\n";
    oss << "k_xcs = " << this->prm_[k_xcs] << "\n";
    oss << "alpha_rel = " << this->prm_[alpha_rel] << "\n";
    oss << "Vol_up = " << this->prm_[Vol_up] << "\n";
    oss << "Vol_rel = " << this->prm_[Vol_rel] << "\n";
    oss << "r_recov = " << this->prm_[r_recov] << "\n";
    oss << "tau_tr = " << this->prm_[tau_tr] << "\n";
    oss << "k_rel_i = " << this->prm_[k_rel_i] << "\n";
    oss << "k_rel_d = " << this->prm_[k_rel_d];   
    return oss.str();

}


std::string Maleckar2009::PrintCurrents() const
{
    using namespace MlcrCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "IKur = " << this->cur_[IKur] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "INab = " << this->cur_[INab] << "\n";
    oss << "ICab = " << this->cur_[ICab] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "ICaP = " << this->cur_[ICaP] << "\n";
    oss << "INaCa = " << this->cur_[INaCa] << "\n";
    oss << "IKACh = " << this->cur_[IKACh] << "\n";
    oss << "Iion = " << this->cur_[MlcrCur::Iion]; 
    return oss.str();

}


std::string Maleckar2009::PrintBlockCoeffs() const
{
    using namespace MlcrCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INa = " << this->block_coeff_[INa] << "\n";
    oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
    oss << "Ito = " << this->block_coeff_[Ito] << "\n";
    oss << "IKur = " << this->block_coeff_[IKur] << "\n";
    oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
    oss << "IKr = " << this->block_coeff_[IKr] << "\n";
    oss << "IKs = " << this->block_coeff_[IKs] << "\n";
    oss << "INab = " << this->block_coeff_[INab] << "\n";
    oss << "ICab = " << this->block_coeff_[ICab] << "\n";
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "ICaP = " << this->block_coeff_[ICaP] << "\n";
    oss << "INaCa = " << this->block_coeff_[INaCa] << "\n";
    oss << "IKACh = " << this->block_coeff_[IKACh];
    return oss.str();

}



} // End of namespace ELECTRA
