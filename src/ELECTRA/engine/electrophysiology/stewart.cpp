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


#include "ELECTRA/engine/electrophysiology/stewart.hpp"


namespace ELECTRA {

void Stewart::SetDataMapping()
{
    using namespace StrtVar;
    using namespace StrtPrm;
    using namespace StrtCur;

    // Set variables mapping.
    this->mapped_data_["v"]        = static_cast<std::size_t>(v);
    this->mapped_data_["K_i"]      = static_cast<std::size_t>(K_i);
    this->mapped_data_["Na_i"]     = static_cast<std::size_t>(Na_i);
    this->mapped_data_["Ca_i"]     = static_cast<std::size_t>(Ca_i);
    this->mapped_data_["y"]        = static_cast<std::size_t>(y);
    this->mapped_data_["Xr1"]      = static_cast<std::size_t>(Xr1);
    this->mapped_data_["Xr2"]      = static_cast<std::size_t>(Xr2);
    this->mapped_data_["Xs"]       = static_cast<std::size_t>(Xs);
    this->mapped_data_["m"]        = static_cast<std::size_t>(m);
    this->mapped_data_["h"]        = static_cast<std::size_t>(h);
    this->mapped_data_["j"]        = static_cast<std::size_t>(j);
    this->mapped_data_["Ca_ss"]    = static_cast<std::size_t>(Ca_ss);
    this->mapped_data_["d"]        = static_cast<std::size_t>(d);
    this->mapped_data_["f"]        = static_cast<std::size_t>(f);
    this->mapped_data_["f2"]       = static_cast<std::size_t>(f2);
    this->mapped_data_["fCass"]    = static_cast<std::size_t>(fCass);
    this->mapped_data_["s"]        = static_cast<std::size_t>(s);
    this->mapped_data_["r"]        = static_cast<std::size_t>(r);
    this->mapped_data_["Ca_SR"]    = static_cast<std::size_t>(Ca_SR);
    this->mapped_data_["R_prime"]  = static_cast<std::size_t>(R_prime);
    this->mapped_data_["V_c"]      = static_cast<std::size_t>(V_c);
    this->mapped_data_["P_kna"]    = static_cast<std::size_t>(P_kna);
    this->mapped_data_["K_o"]      = static_cast<std::size_t>(K_o);
    this->mapped_data_["Na_o"]     = static_cast<std::size_t>(Na_o);
    this->mapped_data_["Ca_o"]     = static_cast<std::size_t>(Ca_o);
    this->mapped_data_["g_f_Na"]   = static_cast<std::size_t>(g_f_Na);
    this->mapped_data_["g_f_K"]    = static_cast<std::size_t>(g_f_K);
    this->mapped_data_["g_f_K1"]   = static_cast<std::size_t>(g_f_K1);
    this->mapped_data_["g_f_Kr"]   = static_cast<std::size_t>(g_f_Kr);
    this->mapped_data_["g_f_Ks"]   = static_cast<std::size_t>(g_f_Ks);
    this->mapped_data_["g_Na"]     = static_cast<std::size_t>(g_Na);
    this->mapped_data_["g_bna"]    = static_cast<std::size_t>(g_bna);
    this->mapped_data_["g_CaL"]    = static_cast<std::size_t>(g_CaL);
    this->mapped_data_["g_bca"]    = static_cast<std::size_t>(g_bca);
    this->mapped_data_["g_to"]     = static_cast<std::size_t>(g_to);
    this->mapped_data_["g_sus"]    = static_cast<std::size_t>(g_sus);
    this->mapped_data_["p_NaK"]    = static_cast<std::size_t>(p_NaK);
    this->mapped_data_["K_mk"]     = static_cast<std::size_t>(K_mk);
    this->mapped_data_["K_mNa"]    = static_cast<std::size_t>(K_mNa);
    this->mapped_data_["K_NaCa"]   = static_cast<std::size_t>(K_NaCa);
    this->mapped_data_["K_sat"]    = static_cast<std::size_t>(K_sat);
    this->mapped_data_["alpha"]    = static_cast<std::size_t>(alpha);
    this->mapped_data_["gamma"]    = static_cast<std::size_t>(gamma);
    this->mapped_data_["Km_Ca"]    = static_cast<std::size_t>(Km_Ca);
    this->mapped_data_["Km_Nai"]   = static_cast<std::size_t>(Km_Nai);
    this->mapped_data_["g_pCa"]    = static_cast<std::size_t>(g_pCa);
    this->mapped_data_["K_pCa"]    = static_cast<std::size_t>(K_pCa);
    this->mapped_data_["g_pK"]     = static_cast<std::size_t>(g_pK);
    this->mapped_data_["k1_prime"] = static_cast<std::size_t>(k1_prime);
    this->mapped_data_["k2_prime"] = static_cast<std::size_t>(k2_prime);
    this->mapped_data_["k3"]       = static_cast<std::size_t>(k3);
    this->mapped_data_["k4"]       = static_cast<std::size_t>(k4);
    this->mapped_data_["EC"]       = static_cast<std::size_t>(EC);
    this->mapped_data_["max_sr"]   = static_cast<std::size_t>(max_sr);
    this->mapped_data_["min_sr"]   = static_cast<std::size_t>(min_sr);
    this->mapped_data_["V_rel"]    = static_cast<std::size_t>(V_rel);
    this->mapped_data_["V_xfer"]   = static_cast<std::size_t>(V_xfer);
    this->mapped_data_["K_up"]     = static_cast<std::size_t>(K_up);
    this->mapped_data_["V_leak"]   = static_cast<std::size_t>(V_leak);
    this->mapped_data_["Vmax_up"]  = static_cast<std::size_t>(Vmax_up);
    this->mapped_data_["Buf_c"]    = static_cast<std::size_t>(Buf_c);
    this->mapped_data_["K_buf_c"]  = static_cast<std::size_t>(K_buf_c);
    this->mapped_data_["Buf_sr"]   = static_cast<std::size_t>(Buf_sr);
    this->mapped_data_["K_buf_sr"] = static_cast<std::size_t>(K_buf_sr);
    this->mapped_data_["Buf_ss"]   = static_cast<std::size_t>(Buf_ss);
    this->mapped_data_["K_buf_ss"] = static_cast<std::size_t>(K_buf_ss);
    this->mapped_data_["V_sr"]     = static_cast<std::size_t>(V_sr);
    this->mapped_data_["V_ss"]     = static_cast<std::size_t>(V_ss);

    // Set parameters mapping.
    this->mapped_data_["R"]           = static_cast<std::size_t>(R);
    this->mapped_data_["T"]           = static_cast<std::size_t>(T);
    this->mapped_data_["Frdy"]        = static_cast<std::size_t>(Frdy);
    this->mapped_data_["Cm"]          = static_cast<std::size_t>(Cm);
    this->mapped_data_["y_inf"]       = static_cast<std::size_t>(y_inf);
    this->mapped_data_["xr1_inf"]     = static_cast<std::size_t>(xr1_inf);
    this->mapped_data_["xr2_inf"]     = static_cast<std::size_t>(xr2_inf);
    this->mapped_data_["xs_inf"]      = static_cast<std::size_t>(xs_inf);
    this->mapped_data_["m_inf"]       = static_cast<std::size_t>(m_inf);
    this->mapped_data_["h_inf"]       = static_cast<std::size_t>(h_inf);
    this->mapped_data_["j_inf"]       = static_cast<std::size_t>(j_inf);
    this->mapped_data_["d_inf"]       = static_cast<std::size_t>(d_inf);
    this->mapped_data_["f_inf"]       = static_cast<std::size_t>(f_inf);
    this->mapped_data_["f2_inf"]      = static_cast<std::size_t>(f2_inf);
    this->mapped_data_["fCass_inf"]   = static_cast<std::size_t>(fCass_inf);
    this->mapped_data_["s_inf"]       = static_cast<std::size_t>(s_inf);
    this->mapped_data_["r_inf"]       = static_cast<std::size_t>(r_inf);
    this->mapped_data_["E_Na"]        = static_cast<std::size_t>(E_Na);
    this->mapped_data_["alpha_y"]     = static_cast<std::size_t>(alpha_y);
    this->mapped_data_["alpha_xr1"]   = static_cast<std::size_t>(alpha_xr1);
    this->mapped_data_["alpha_xr2"]   = static_cast<std::size_t>(alpha_xr2);
    this->mapped_data_["alpha_xs"]    = static_cast<std::size_t>(alpha_xs);
    this->mapped_data_["alpha_m"]     = static_cast<std::size_t>(alpha_m);
    this->mapped_data_["alpha_h"]     = static_cast<std::size_t>(alpha_h);
    this->mapped_data_["alpha_j"]     = static_cast<std::size_t>(alpha_j);
    this->mapped_data_["alpha_d"]     = static_cast<std::size_t>(alpha_d);
    this->mapped_data_["tau_f"]       = static_cast<std::size_t>(tau_f);
    this->mapped_data_["tau_f2"]      = static_cast<std::size_t>(tau_f2);
    this->mapped_data_["tau_fCass"]   = static_cast<std::size_t>(tau_fCass);
    this->mapped_data_["tau_s"]       = static_cast<std::size_t>(tau_s);
    this->mapped_data_["tau_r"]       = static_cast<std::size_t>(tau_r);
    this->mapped_data_["E_K"]         = static_cast<std::size_t>(E_K);
    this->mapped_data_["beta_y"]      = static_cast<std::size_t>(beta_y);
    this->mapped_data_["beta_xr1"]    = static_cast<std::size_t>(beta_xr1);
    this->mapped_data_["beta_xr2"]    = static_cast<std::size_t>(beta_xr2);
    this->mapped_data_["beta_xs"]     = static_cast<std::size_t>(beta_xs);
    this->mapped_data_["beta_m"]      = static_cast<std::size_t>(beta_m);
    this->mapped_data_["beta_h"]      = static_cast<std::size_t>(beta_h);
    this->mapped_data_["beta_j"]      = static_cast<std::size_t>(beta_j);
    this->mapped_data_["beta_d"]      = static_cast<std::size_t>(beta_d);
    this->mapped_data_["E_Ks"]        = static_cast<std::size_t>(E_Ks);
    this->mapped_data_["tau_y"]       = static_cast<std::size_t>(tau_y);
    this->mapped_data_["tau_xr1"]     = static_cast<std::size_t>(tau_xr1);
    this->mapped_data_["tau_xr2"]     = static_cast<std::size_t>(tau_xr2);
    this->mapped_data_["tau_xs"]      = static_cast<std::size_t>(tau_xs);
    this->mapped_data_["tau_m"]       = static_cast<std::size_t>(tau_m);
    this->mapped_data_["tau_h"]       = static_cast<std::size_t>(tau_h);
    this->mapped_data_["tau_j"]       = static_cast<std::size_t>(tau_j);
    this->mapped_data_["gamma_d"]     = static_cast<std::size_t>(gamma_d);
    this->mapped_data_["E_Ca"]        = static_cast<std::size_t>(E_Ca);
    this->mapped_data_["tau_d"]       = static_cast<std::size_t>(tau_d);
    this->mapped_data_["xK1_inf"]     = static_cast<std::size_t>(xK1_inf);
    this->mapped_data_["a"]           = static_cast<std::size_t>(a);
    this->mapped_data_["kcasr"]       = static_cast<std::size_t>(kcasr);
    this->mapped_data_["Ca_i_bufc"]   = static_cast<std::size_t>(Ca_i_bufc),
    this->mapped_data_["k1"]          = static_cast<std::size_t>(k1);
    this->mapped_data_["k2"]          = static_cast<std::size_t>(k2);
    this->mapped_data_["O"]           = static_cast<std::size_t>(O);
    this->mapped_data_["Ca_sr_bufsr"] = static_cast<std::size_t>(Ca_sr_bufsr);
    this->mapped_data_["Ca_ss_bufss"] = static_cast<std::size_t>(Ca_ss_bufss);

    // Set currents mapping.
    this->mapped_data_["INa"]   = static_cast<std::size_t>(INa);
    this->mapped_data_["IfNa"]  = static_cast<std::size_t>(IfNa);
    this->mapped_data_["IfK"]   = static_cast<std::size_t>(IfK);
    this->mapped_data_["If"]    = static_cast<std::size_t>(If);
    this->mapped_data_["Irel"]  = static_cast<std::size_t>(Irel);
    this->mapped_data_["Isus"]  = static_cast<std::size_t>(Isus);
    this->mapped_data_["INaK"]  = static_cast<std::size_t>(INaK);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["IpCa"]  = static_cast<std::size_t>(IpCa);
    this->mapped_data_["IpK"]   = static_cast<std::size_t>(IpK);
    this->mapped_data_["Iup"]   = static_cast<std::size_t>(Iup);
    this->mapped_data_["Ileak"] = static_cast<std::size_t>(Ileak);
    this->mapped_data_["Ixfer"] = static_cast<std::size_t>(Ixfer);
    this->mapped_data_["IK1"]   = static_cast<std::size_t>(IK1);
    this->mapped_data_["IKr"]   = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"]   = static_cast<std::size_t>(IKs);
    this->mapped_data_["INa"]   = static_cast<std::size_t>(INa);
    this->mapped_data_["IbNa"]  = static_cast<std::size_t>(IbNa);
    this->mapped_data_["ICaL"]  = static_cast<std::size_t>(ICaL);
    this->mapped_data_["IbCa"]  = static_cast<std::size_t>(IbCa);
    this->mapped_data_["Ito"]   = static_cast<std::size_t>(Ito);

}


Stewart::Stewart()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Stewart;
    this->dt_stable_ = 0.02;
    this->var_.resize(69, 0.);
    this->prm_.resize(60, 0.);
    this->cur_.resize(21, 0.);
    this->block_coeff_.resize(20, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


Stewart::~Stewart()
{}


void Stewart::Initialize(CellType cell_type)
{
    using namespace StrtVar;
    using namespace StrtPrm;

    if (cell_type != CellType::purkinje) {
        std::string error_str = "Could not initialize the Stewart purkinje ap model. Expected cell type: ELECTRA::CellType::purkinje";
        throw std::invalid_argument(Logger::Error(error_str));
    }

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(69, 0.);
    this->prm_.clear();           this->prm_.resize(60, 0.);
    this->cur_.clear();           this->cur_.resize(21, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(20, 0.);

    // Set the model variables.
    this->var_[v] = -75.6120082079284;
    this->var_[dvdt] = 0.;
    this->var_[K_i] = 136.757624998959;
    this->var_[Na_i] = 8.8010751354651;
    this->var_[Ca_i] = 0.000148302222806612;
    this->var_[y] = 0.00757708477954872;
    this->var_[Xr1] = 0.39149983054504;
    this->var_[Xr2] = 0.373763519962402;
    this->var_[Xs] = 0.0397051345980805;
    this->var_[m] = 0.012407375005411;
    this->var_[h] = 0.356980249571247;
    this->var_[j] = 0.097975616674561;
    this->var_[Ca_ss] = 0.000555739066351023;
    this->var_[d] = 0.000121553774479047;
    this->var_[f] = 0.605783874119932;
    this->var_[f2] = 0.856306851055087;
    this->var_[fCass] = 0.985231302968553;
    this->var_[s] = 0.922552370758773;
    this->var_[r] = 0.000648943716577986;
    this->var_[Ca_SR] = 3.17317827043521;
    this->var_[R_prime] = 0.849862133544873;
    this->var_[V_c] = 0.016404;
    this->var_[P_kna] = 0.03;
    this->var_[K_o] = 5.4;
    this->var_[Na_o] = 140.;
    this->var_[Ca_o] = 2.;
    this->var_[g_f_Na] = 0.0145654;
    this->var_[g_f_K] = 0.0234346;
    this->var_[g_f_K1] = 0.065;
    this->var_[g_f_Kr] = 0.0918;
    this->var_[g_f_Ks] = 0.2352;
    this->var_[g_Na] = 130.5744;
    this->var_[g_bna] = 0.00029;
    this->var_[g_CaL] = 3.98e-05;
    this->var_[g_bca] = 0.000592;
    this->var_[g_to] = 0.08184;
    this->var_[g_sus] = 0.0227;
    this->var_[p_NaK] = 2.724;
    this->var_[K_mk] = 1.;
    this->var_[K_mNa] = 40.;
    this->var_[K_NaCa] = 1000.;
    this->var_[K_sat] = 0.1;
    this->var_[alpha] = 2.5;
    this->var_[gamma] = 0.35;
    this->var_[Km_Ca] = 1.38;
    this->var_[Km_Nai] = 87.5;
    this->var_[g_pCa] = 0.1238;
    this->var_[K_pCa] = 0.0005;
    this->var_[g_pK] = 0.0146;
    this->var_[k1_prime] = 0.15;
    this->var_[k2_prime] = 0.045;
    this->var_[k3] = 0.06;
    this->var_[k4] = 0.005;
    this->var_[EC] = 1.5;
    this->var_[max_sr] = 2.5;
    this->var_[min_sr] = 1.;
    this->var_[V_rel] = 0.102;
    this->var_[V_xfer] = 0.0038;
    this->var_[K_up] = 0.00025;
    this->var_[V_leak] = 0.00036;
    this->var_[Vmax_up] = 0.006375;
    this->var_[Buf_c] = 0.2;
    this->var_[K_buf_c] = 0.001;
    this->var_[Buf_sr] = 10.;
    this->var_[K_buf_sr] = 0.3;
    this->var_[Buf_ss] = 0.4;
    this->var_[K_buf_ss] = 0.00025;
    this->var_[V_sr] = 0.001094;
    this->var_[V_ss] = 5.468e-05;
 
    // Set the model parameters.
    this->prm_[R] = 8314.472;
    this->prm_[T] = 310.;
    this->prm_[Frdy] = 96485.3415;
    this->prm_[Cm] = 0.185;
    this->prm_[y_inf] = 0.324421570595336;
    this->prm_[xr1_inf] = 0.000834842379227446;
    this->prm_[xr2_inf] = 0.373749294458218;
    this->prm_[xs_inf] = 0.00640840263360323;
    this->prm_[m_inf] = 0.0124073892984643;
    this->prm_[h_inf] = 0.40115890572671;
    this->prm_[j_inf] = 0.40115890572671;
    this->prm_[d_inf] = 0.000121565359096558;
    this->prm_[f_inf] = 0.999645545227779;
    this->prm_[f2_inf] = 0.997981122282696;
    this->prm_[fCass_inf] = 0.999925880158706;
    this->prm_[s_inf] = 0.976783157422829;
    this->prm_[r_inf] = 0.000639124915845885;
    this->prm_[E_Na] = 73.9107934751248;
    this->prm_[alpha_y] = 1.13256059179852;
    this->prm_[alpha_xr1] = 20.131324929499;
    this->prm_[alpha_xr2] = 0.942570887821727;
    this->prm_[alpha_xs] = 1.69309441324577;
    this->prm_[alpha_m] = 0.0421924571636243;
    this->prm_[alpha_h] = 0.0298970890489957;
    this->prm_[alpha_j] = 0.00400443627029917;
    this->prm_[alpha_d] = 0.30897917053114;
    this->prm_[tau_f] = 198.197171197593;
    this->prm_[tau_f2] = 79.2037617647695;
    this->prm_[tau_fCass] = 81.9901173544941;
    this->prm_[tau_s] = 47.0283741219712;
    this->prm_[tau_r] = 12.465682714962;
    this->prm_[E_K] = -86.33383203891;
    this->prm_[beta_y] = 0.00893907673161732;
    this->prm_[beta_xr1] = 5.88844869212083;
    this->prm_[beta_xr2] = 1.11872957985673;
    this->prm_[beta_xs] = 0.999373124251383;
    this->prm_[beta_m] = 0.16517527478954;
    this->prm_[beta_h] = 0.00687450167743142;
    this->prm_[beta_j] = 0.000401675585180273;
    this->prm_[beta_d] = 1.39999896998161;
    this->prm_[E_Ks] = -71.015217402632;
    this->prm_[tau_y] = 3504.16220895677;
    this->prm_[tau_xr1] = 118.542273951768;
    this->prm_[tau_xr2] = 1.05448193331798;
    this->prm_[tau_xs] = 81.692033053418;
    this->prm_[tau_m] = 0.00696915070604752;
    this->prm_[tau_h] = 27.1949072706642;
    this->prm_[tau_j] = 226.957470168716;
    this->prm_[gamma_d] = 0.00186877550382267;
    this->prm_[E_Ca] = 127.015612247426;
    this->prm_[tau_d] = 0.434439295993189;
    this->prm_[xK1_inf] = 0.504300568716195;
    this->prm_[a] = 0.00864677966337645;
    this->prm_[kcasr] = 1.27396628532648;
    this->prm_[Ca_i_bufc] = 0.00654985485720445;
    this->prm_[k1] = 0.117742519348979;
    this->prm_[k2] = 0.0573284828396915;
    this->prm_[O] = 5.15118540890581e-07;
    this->prm_[Ca_sr_bufsr] = 0.800835674088394;
    this->prm_[Ca_ss_bufss] = 0.00645063485338993;

}


void Stewart::Compute(double v_new, double dt, double stim_current)
{
    using namespace StrtVar;
    using namespace StrtPrm;
    using namespace StrtCur;

    this->prm_[f_inf] = 1./(1.+std::exp((v_new+20.)/7.));
    this->prm_[tau_f] =  1102.5*std::exp(-std::pow(v_new+27., 2.)/225.)+200./(1.+std::exp((13. - v_new)/10.))+180./(1.+std::exp((v_new+30.)/10.))+20.;
    this->var_[f] = this->prm_[f_inf] - (this->prm_[f_inf] - this->var_[f])*std::exp(-dt/this->prm_[tau_f]);
    
    this->prm_[f2_inf] = 0.67/(1.+std::exp((v_new+35.)/7.))+0.33;
    this->prm_[tau_f2] =  562.*std::exp( - std::pow(v_new+27., 2.)/240.)+31./(1.+exp((25. - v_new)/10.))+80./(1.+exp((v_new+30.)/10.0000));
    this->var_[f2] = this->prm_[f2_inf] - (this->prm_[f2_inf] - this->var_[f2])*std::exp(-dt/this->prm_[tau_f2]);
    
    this->prm_[fCass_inf] = 0.6/(1.+std::pow(this->var_[Ca_ss]/0.05, 2.)) + 0.4;
    this->prm_[tau_fCass] = 80./(1.+std::pow(this->var_[Ca_ss]/0.05, 2.)) + 2.;
    this->var_[fCass] = this->prm_[fCass_inf] - (this->prm_[fCass_inf] - this->var_[fCass])*std::exp(-dt/this->prm_[tau_fCass]);
    
    this->prm_[s_inf] = 1./(1.+std::exp((v_new+27.)/13.));
    this->prm_[tau_s] =  85.*std::exp(-std::pow(v_new+25., 2.)/320.)+5./(1.+exp((v_new - 40.)/5.))+42.;
    this->var_[s] = this->prm_[s_inf] - (this->prm_[s_inf] - this->var_[s])*std::exp(-dt/this->prm_[tau_s]);
    
    this->prm_[r_inf] = 1./(1.+std::exp((20. - v_new)/13.));
    this->prm_[tau_r] =  10.45*std::exp(-std::pow(v_new+40., 2.)/1800.) + 7.3;
    this->var_[r] = this->prm_[r_inf] - (this->prm_[r_inf] - this->var_[r])*std::exp(-dt/this->prm_[tau_r]);
    
    this->prm_[y_inf] = 1./(1.+std::exp((v_new+80.6)/6.8));
    this->prm_[alpha_y] =  1.*std::exp(-2.9 - 0.04*v_new);
    this->prm_[beta_y] =  1.*std::exp(3.6 + 0.11*v_new);
    this->prm_[tau_y] = 4000./(this->prm_[alpha_y]+this->prm_[beta_y]);
    this->var_[y] = this->prm_[y_inf] - (this->prm_[y_inf] - this->var_[y])*std::exp(-dt/this->prm_[tau_y]);
    
    this->prm_[xr1_inf] = 1./(1.+std::exp((-26. - v_new)/7.));
    this->prm_[alpha_xr1] = 450./(1.+std::exp((-45. - v_new)/10.));
    this->prm_[beta_xr1] = 6./(1.+std::exp((v_new+30.)/11.5));
    this->prm_[tau_xr1] =  this->prm_[alpha_xr1] * this->prm_[beta_xr1];
    this->var_[Xr1] = this->prm_[xr1_inf] - (this->prm_[xr1_inf] - this->var_[Xr1])*std::exp(-dt/this->prm_[tau_xr1]);
    
    this->prm_[xr2_inf] = 1./(1.+std::exp((v_new+88.)/24.));
    this->prm_[alpha_xr2] = 3./(1.+std::exp((-60. - v_new)/20.));
    this->prm_[beta_xr2] = 1.12/(1.+std::exp((v_new - 60.)/20.));
    this->prm_[tau_xr2] =  this->prm_[alpha_xr2] * this->prm_[beta_xr2];
    this->var_[Xr2] = this->prm_[xr2_inf] - (this->prm_[xr2_inf] - this->var_[Xr2])*std::exp(-dt/this->prm_[tau_xr2]);
    
    this->prm_[xs_inf] = 1./(1.+std::exp((-5. - v_new)/14.));
    this->prm_[alpha_xs] = 1400./std::pow((1.+std::exp((5. - v_new)/6.)), 0.5);
    this->prm_[beta_xs] = 1./(1.+std::exp((v_new - 35.)/15.));
    this->prm_[tau_xs] =  this->prm_[alpha_xs] * this->prm_[beta_xs] + 80.;
    this->var_[Xs] = this->prm_[xs_inf] - (this->prm_[xs_inf] - this->var_[Xs])*std::exp(-dt/this->prm_[tau_xs]);
    
    this->prm_[m_inf] = 1./std::pow(1.+std::exp((-56.86 - v_new)/9.03), 2.);
    this->prm_[alpha_m] = 1./(1.+std::exp((-60. - v_new)/5.));
    this->prm_[beta_m] = 0.1/(1.+std::exp((v_new+35.)/5.))+0.1/(1.+std::exp((v_new - 50.)/200.));
    this->prm_[tau_m] = this->prm_[alpha_m] * this->prm_[beta_m];
    this->var_[m] = this->prm_[m_inf] - (this->prm_[m_inf] - this->var_[m])*std::exp(-dt/this->prm_[tau_m]);
    
    this->prm_[h_inf] = 1./std::pow(1.+std::exp((v_new+71.55)/7.43), 2.);

    if (v_new < -40) {
        this->prm_[alpha_h] = 0.057*std::exp(-(v_new+80.)/6.8);
        this->prm_[beta_h] = 2.7*std::exp(0.079*v_new) + 310000*std::exp(0.348500*v_new);
    } else {
        this->prm_[alpha_h] = 0.;
        this->prm_[beta_h] = 0.77 / (0.13*(1.+std::exp((v_new+10.66)/ - 11.1)));
    }

    this->prm_[tau_h] = 1./(this->prm_[alpha_h]+this->prm_[beta_h]);
    this->var_[h] = ALGORITHM::RushLarsen(this->prm_[h_inf], this->var_[h], dt, this->prm_[tau_h]);

    this->prm_[j_inf] = 1./std::pow(1.+std::exp((v_new+71.55)/7.43), 2.);
    
    if (v_new < -40.) {
        this->prm_[alpha_j] = (((-25428.*std::exp(0.2444*v_new) - 6.948e-6*std::exp(-0.04391*v_new))*(v_new+37.78)))/(1.+std::exp(0.311*(v_new+79.23)));
        this->prm_[beta_j] = (0.02424*std::exp(-0.01052*v_new))/(1.+std::exp(-0.1378*(v_new+40.14)));
    } else {
        this->prm_[alpha_j] = 0.;
        this->prm_[beta_j] = (0.6*std::exp(0.057*v_new))/(1.+std::exp(-0.1*(v_new+32.)));
    }

    this->prm_[tau_j] = 1./(this->prm_[alpha_j]+this->prm_[beta_j]);
    this->var_[j] = ALGORITHM::RushLarsen(this->prm_[j_inf], this->var_[j], dt, this->prm_[tau_j]);
    
    this->prm_[d_inf] = 1./(1.+std::exp((-8. - v_new)/7.5));
    this->prm_[alpha_d] = 1.4/(1.+std::exp((-35. - v_new)/13.))+0.25;
    this->prm_[beta_d] = 1.4/(1.+std::exp((v_new+5.)/5.));
    this->prm_[gamma_d] = 1./(1.+std::exp((50. - v_new)/20.));
    this->prm_[tau_d] = this->prm_[alpha_d] * this->prm_[beta_d] + this->prm_[gamma_d];
    this->var_[d] = ALGORITHM::RushLarsen(this->prm_[d_inf], this->var_[d], dt, this->prm_[tau_d]);
    
    this->prm_[E_Na] = ((this->prm_[R]*this->prm_[T])/this->prm_[Frdy]) * std::log(this->var_[Na_o]/this->var_[Na_i]);

    this->cur_[INaK] = ((((this->var_[p_NaK] * this->var_[K_o]) / (this->var_[K_o]+this->var_[K_mk]))*this->var_[Na_i]) / (this->var_[Na_i]+this->var_[K_mNa])) / (1.+0.1245*std::exp((-0.1*v_new*this->prm_[Frdy])/(this->prm_[R]*this->prm_[T])) + 0.0353*std::exp((-v_new*this->prm_[Frdy])/(this->prm_[R]*this->prm_[T])));
    this->cur_[INa] = this->var_[g_Na]*std::pow(this->var_[m], 3.)*this->var_[h]*this->var_[j]*(v_new - this->prm_[E_Na]);
    this->cur_[IbNa] = this->var_[g_bna]*(v_new - this->prm_[E_Na]);
    this->cur_[INaCa] = (this->var_[K_NaCa]*(std::exp((this->var_[gamma]*v_new*this->prm_[Frdy])/(this->prm_[R]*this->prm_[T]))*std::pow(this->var_[Na_i], 3.)*this->var_[Ca_o] - std::exp(((this->var_[gamma] - 1.)*v_new*this->prm_[Frdy])/(this->prm_[R]*this->prm_[T]))*std::pow(this->var_[Na_o], 3.)*this->var_[Ca_i]*this->var_[alpha]))/((std::pow(this->var_[Km_Nai], 3.)+std::pow(this->var_[Na_o], 3.))*(this->var_[Km_Ca]+this->var_[Ca_o])*(1.+this->var_[K_sat]*std::exp(((this->var_[gamma] - 1.)*v_new*this->prm_[Frdy])/(this->prm_[R]*this->prm_[T]))));
    this->cur_[IfNa] =  this->var_[y]*this->var_[g_f_Na]*(v_new - this->prm_[E_Na]);
    
    double dNa_i =  ((-1.*(this->cur_[INa] + this->cur_[IbNa] + this->cur_[IfNa] + 3.*this->cur_[INaK] + 3.*this->cur_[INaCa]))/(1.*this->var_[V_c]*this->prm_[Frdy]))*this->prm_[Cm];
    this->var_[Na_i] += dt*dNa_i;
    
    this->prm_[E_K] =  ((this->prm_[R]*this->prm_[T])/this->prm_[Frdy]) * std::log(this->var_[K_o]/this->var_[K_i]);
    this->prm_[xK1_inf] = 1./(1.+std::exp(0.1*(v_new+75.44)));
    this->cur_[IK1] = this->var_[g_f_K1]*this->prm_[xK1_inf]*((v_new - 8.) - this->prm_[E_K]);
    this->cur_[Ito] = this->var_[g_to]*this->var_[r]*this->var_[s]*(v_new - this->prm_[E_K]);
    
    this->prm_[a] = 1./(1.+std::exp((5. - v_new)/17.));
    this->cur_[Isus] = this->var_[g_sus]*this->prm_[a]*(v_new - this->prm_[E_K]);
    this->cur_[IKr] = this->var_[g_f_Kr]*std::pow((this->var_[K_o]/5.4), 0.5)*this->var_[Xr1]*this->var_[Xr2]*(v_new - this->prm_[E_K]);
    
    this->prm_[E_Ks] = ((this->prm_[R]*this->prm_[T])/this->prm_[Frdy]) * std::log((this->var_[K_o] + this->var_[P_kna]*this->var_[Na_o])/(this->var_[K_i] + this->var_[P_kna]*this->var_[Na_i]));
    this->cur_[IKs] = this->var_[g_f_Ks]*std::pow(this->var_[Xs], 2.)*(v_new - this->prm_[E_Ks]);
    this->cur_[ICaL] = (((this->var_[g_CaL]*this->var_[d]*this->var_[f]*this->var_[f2]*this->var_[fCass]*4.*(v_new - 15.)*std::pow(this->prm_[Frdy], 2.))/(this->prm_[R]*this->prm_[T]))*(0.25*this->var_[Ca_ss]*std::exp((2.*(v_new - 15.)*this->prm_[Frdy])/(this->prm_[R]*this->prm_[T])) - this->var_[Ca_o]))/(std::exp((2.*(v_new - 15.)*this->prm_[Frdy])/(this->prm_[R]*this->prm_[T])) - 1.);
    
    this->prm_[E_Ca] = ((0.5*this->prm_[R]*this->prm_[T])/this->prm_[Frdy]) * std::log(this->var_[Ca_o]/this->var_[Ca_i]);
    this->cur_[IbCa] = this->var_[g_bca]*(v_new - this->prm_[E_Ca]);
    this->cur_[IpK] = (this->var_[g_pK]*(v_new - this->prm_[E_K]))/(1.+std::exp((25. - v_new)/5.98));
    this->cur_[IpCa] = (this->var_[g_pCa]*this->var_[Ca_i])/(this->var_[Ca_i]+this->var_[K_pCa]);
    this->cur_[IfK] = this->var_[y]*this->var_[g_f_K]*(v_new - this->prm_[E_K]);
    this->cur_[If] = this->cur_[IfNa] + this->cur_[IfK];
    
    this->cur_[StrtCur::Iion] = this->cur_[IK1] + this->cur_[Ito] + this->cur_[Isus] + this->cur_[IKr] + this->cur_[IKs] + this->cur_[ICaL] + this->cur_[INaK] + this->cur_[INa] + this->cur_[IbNa] + this->cur_[INaCa] + this->cur_[IbCa] + this->cur_[IpK] + this->cur_[IpCa] + this->cur_[If];
    this->var_[dvdt] =  - (this->cur_[StrtCur::Iion] - stim_current);
    
    double dK_i =  ((-1.*((this->cur_[IK1] + this->cur_[Ito] + this->cur_[IfK] + this->cur_[Isus] + this->cur_[IKr] + this->cur_[IKs] + this->cur_[IpK]) - 2.*this->cur_[INaK]))/(1.*this->var_[V_c]*this->prm_[Frdy]))*this->prm_[Cm];
    this->var_[K_i] += dt*dK_i;
    
    this->cur_[Iup] = this->var_[Vmax_up]/(1.+std::pow(this->var_[K_up], 2.)/std::pow(this->var_[Ca_i], 2.));
    this->cur_[Ileak] = this->var_[V_leak]*(this->var_[Ca_SR] - this->var_[Ca_i]);
    this->cur_[Ixfer] = this->var_[V_xfer]*(this->var_[Ca_ss] - this->var_[Ca_i]);
    this->prm_[Ca_i_bufc] = 1./(1.+(this->var_[Buf_c]*this->var_[K_buf_c])/std::pow(this->var_[Ca_i]+this->var_[K_buf_c], 2.));
    double dCa_i = this->prm_[Ca_i_bufc]*((((this->cur_[Ileak] - this->cur_[Iup])*this->var_[V_sr])/this->var_[V_c]+this->cur_[Ixfer]) - (((this->cur_[IbCa]+this->cur_[IpCa]) - 2.*this->cur_[INaCa])*this->prm_[Cm])/(2.*this->var_[V_c]*this->prm_[Frdy]));
    this->var_[Ca_i] += dt*dCa_i;
    
    this->prm_[kcasr] = this->var_[max_sr] - (this->var_[max_sr] - this->var_[min_sr])/(1.+std::pow(this->var_[EC]/this->var_[Ca_SR], 2.));
    this->prm_[k2] = this->var_[k2_prime]*this->prm_[kcasr];
    double dR_prime = -this->prm_[k2]*this->var_[Ca_ss]*this->var_[R_prime] + this->var_[k4]*(1. - this->var_[R_prime]);
    this->var_[R_prime] += dt*dR_prime;
    
    this->prm_[k1] = this->var_[k1_prime] / this->prm_[kcasr];
    this->prm_[O] = (this->prm_[k1]*std::pow(this->var_[Ca_ss], 2.)*this->var_[R_prime])/(this->var_[k3] + this->prm_[k1]*std::pow(this->var_[Ca_ss], 2.));
    this->cur_[Irel] = this->var_[V_rel]*this->prm_[O]*(this->var_[Ca_SR] - this->var_[Ca_ss]);
    
    this->prm_[Ca_sr_bufsr] = 1./(1.+(this->var_[Buf_sr]*this->var_[K_buf_sr])/std::pow(this->var_[Ca_SR]+this->var_[K_buf_sr], 2.));
    double dCa_SR = this->prm_[Ca_sr_bufsr]*(this->cur_[Iup] - (this->cur_[Irel] + this->cur_[Ileak]));
    this->var_[Ca_SR] += dt*dCa_SR;
    
    this->prm_[Ca_ss_bufss] = 1./(1.+(this->var_[Buf_ss]*this->var_[K_buf_ss])/std::pow(this->var_[Ca_ss]+this->var_[K_buf_ss], 2.));
    double dCa_ss = this->prm_[Ca_ss_bufss]*(((-1.*this->cur_[ICaL]*this->prm_[Cm])/(2.*this->var_[V_ss]*this->prm_[Frdy]) + (this->cur_[Irel]*this->var_[V_sr])/this->var_[V_ss]) - (this->cur_[Ixfer]*this->var_[V_c])/this->var_[V_ss]);
    this->var_[Ca_ss] += dt*dCa_ss;
    
}


std::string Stewart::PrintVariables() const 
{
    using namespace StrtVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "K_i = " << this->var_[K_i] << "\n";
    oss << "Na_i = " << this->var_[Na_i] << "\n";
    oss << "Ca_i = " << this->var_[Ca_i] << "\n";
    oss << "y = " << this->var_[y] << "\n";
    oss << "Xr1 = " << this->var_[Xr1] << "\n";
    oss << "Xr2 = " << this->var_[Xr2] << "\n";
    oss << "Xs = " << this->var_[Xs] << "\n";
    oss << "m = " << this->var_[m] << "\n";
    oss << "h = " << this->var_[h] << "\n";
    oss << "j = " << this->var_[j] << "\n";
    oss << "Ca_ss = " << this->var_[Ca_ss] << "\n";
    oss << "d = " << this->var_[d] << "\n";
    oss << "f = " << this->var_[f] << "\n";
    oss << "f2 = " << this->var_[f2] << "\n";
    oss << "fCass = " << this->var_[fCass] << "\n";
    oss << "s = " << this->var_[s] << "\n";
    oss << "r = " << this->var_[r] << "\n";
    oss << "Ca_SR = " << this->var_[Ca_SR] << "\n";
    oss << "R_prime = " << this->var_[R_prime] << "\n";
    oss << "V_c = " << this->var_[V_c] << "\n";
    oss << "P_kna = " << this->var_[P_kna] << "\n";
    oss << "K_o = " << this->var_[K_o] << "\n";
    oss << "Na_o = " << this->var_[Na_o] << "\n";
    oss << "Ca_o = " << this->var_[Ca_o] << "\n";
    oss << "g_f_Na = " << this->var_[g_f_Na] << "\n";
    oss << "g_f_K = " << this->var_[g_f_K] << "\n";
    oss << "g_f_K1 = " << this->var_[g_f_K1] << "\n";
    oss << "g_f_Kr = " << this->var_[g_f_Kr] << "\n";
    oss << "g_f_Ks = " << this->var_[g_f_Ks] << "\n";
    oss << "g_Na = " << this->var_[g_Na] << "\n";
    oss << "g_bna = " << this->var_[g_bna] << "\n";
    oss << "g_CaL = " << this->var_[g_CaL] << "\n";
    oss << "g_bca = " << this->var_[g_bca] << "\n";
    oss << "g_to = " << this->var_[g_to] << "\n";
    oss << "g_sus = " << this->var_[g_sus] << "\n";
    oss << "p_NaK = " << this->var_[p_NaK] << "\n";
    oss << "K_mk = " << this->var_[K_mk] << "\n";
    oss << "K_mNa = " << this->var_[K_mNa] << "\n";
    oss << "K_NaCa = " << this->var_[K_NaCa] << "\n";
    oss << "K_sat = " << this->var_[K_sat] << "\n";
    oss << "alpha = " << this->var_[alpha] << "\n";
    oss << "gamma = " << this->var_[gamma] << "\n";
    oss << "Km_Ca = " << this->var_[Km_Ca] << "\n";
    oss << "Km_Nai = " << this->var_[Km_Nai] << "\n";
    oss << "g_pCa = " << this->var_[g_pCa] << "\n";
    oss << "K_pCa = " << this->var_[K_pCa] << "\n";
    oss << "g_pK = " << this->var_[g_pK] << "\n";
    oss << "k1_prime = " << this->var_[k1_prime] << "\n";
    oss << "k2_prime = " << this->var_[k2_prime] << "\n";
    oss << "k3 = " << this->var_[k3] << "\n";
    oss << "k4 = " << this->var_[k4] << "\n";
    oss << "EC = " << this->var_[EC] << "\n";
    oss << "max_sr = " << this->var_[max_sr] << "\n";
    oss << "min_sr = " << this->var_[min_sr] << "\n";
    oss << "V_rel = " << this->var_[V_rel] << "\n";
    oss << "V_xfer = " << this->var_[V_xfer] << "\n";
    oss << "K_up = " << this->var_[K_up] << "\n";
    oss << "V_leak = " << this->var_[V_leak] << "\n";
    oss << "Vmax_up = " << this->var_[Vmax_up] << "\n";
    oss << "Buf_c = " << this->var_[Buf_c] << "\n";
    oss << "K_buf_c = " << this->var_[K_buf_c] << "\n";
    oss << "Buf_sr = " << this->var_[Buf_sr] << "\n";
    oss << "K_buf_sr = " << this->var_[K_buf_sr] << "\n";
    oss << "Buf_ss = " << this->var_[Buf_ss] << "\n";
    oss << "K_buf_ss = " << this->var_[K_buf_ss] << "\n";
    oss << "V_sr = " << this->var_[V_sr] << "\n";
    oss << "V_ss = " << this->var_[V_ss];
    return oss.str();

}


std::string Stewart::PrintParameters() const 
{
    using namespace StrtPrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "R = " << this->prm_[R] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "Frdy = " << this->prm_[Frdy] << "\n";
    oss << "Cm = " << this->prm_[Cm] << "\n";
    oss << "y_inf = " << this->prm_[y_inf] << "\n";
    oss << "xr1_inf = " << this->prm_[xr1_inf] << "\n";
    oss << "xr2_inf = " << this->prm_[xr2_inf] << "\n";
    oss << "xs_inf = " << this->prm_[xs_inf] << "\n";
    oss << "m_inf = " << this->prm_[m_inf] << "\n";
    oss << "h_inf = " << this->prm_[h_inf] << "\n";
    oss << "j_inf = " << this->prm_[j_inf] << "\n";
    oss << "d_inf = " << this->prm_[d_inf] << "\n";
    oss << "f_inf = " << this->prm_[f_inf] << "\n";
    oss << "f2_inf = " << this->prm_[f2_inf] << "\n";
    oss << "fCass_inf = " << this->prm_[fCass_inf] << "\n";
    oss << "s_inf = " << this->prm_[s_inf] << "\n";
    oss << "r_inf = " << this->prm_[r_inf] << "\n";
    oss << "E_Na = " << this->prm_[E_Na] << "\n";
    oss << "alpha_y = " << this->prm_[alpha_y] << "\n";
    oss << "alpha_xr1 = " << this->prm_[alpha_xr1] << "\n";
    oss << "alpha_xr2 = " << this->prm_[alpha_xr2] << "\n";
    oss << "alpha_xs = " << this->prm_[alpha_xs] << "\n";
    oss << "alpha_m = " << this->prm_[alpha_m] << "\n";
    oss << "alpha_h = " << this->prm_[alpha_h] << "\n";
    oss << "alpha_j = " << this->prm_[alpha_j] << "\n";
    oss << "alpha_d = " << this->prm_[alpha_d] << "\n";
    oss << "tau_f = " << this->prm_[tau_f] << "\n";
    oss << "tau_f2 = " << this->prm_[tau_f2] << "\n";
    oss << "tau_fCass = " << this->prm_[tau_fCass] << "\n";
    oss << "tau_s = " << this->prm_[tau_s] << "\n";
    oss << "tau_r = " << this->prm_[tau_r] << "\n";
    oss << "E_K = " << this->prm_[E_K] << "\n";
    oss << "beta_y = " << this->prm_[beta_y] << "\n";
    oss << "beta_xr1 = " << this->prm_[beta_xr1] << "\n";
    oss << "beta_xr2 = " << this->prm_[beta_xr2] << "\n";
    oss << "beta_xs = " << this->prm_[beta_xs] << "\n";
    oss << "beta_m = " << this->prm_[beta_m] << "\n";
    oss << "beta_h = " << this->prm_[beta_h] << "\n";
    oss << "beta_j = " << this->prm_[beta_j] << "\n";
    oss << "beta_d = " << this->prm_[beta_d] << "\n";
    oss << "E_Ks = " << this->prm_[E_Ks] << "\n";
    oss << "tau_y = " << this->prm_[tau_y] << "\n";
    oss << "tau_xr1 = " << this->prm_[tau_xr1] << "\n";
    oss << "tau_xr2 = " << this->prm_[tau_xr2] << "\n";
    oss << "tau_xs = " << this->prm_[tau_xs] << "\n";
    oss << "tau_m = " << this->prm_[tau_m] << "\n";
    oss << "tau_h = " << this->prm_[tau_h] << "\n";
    oss << "tau_j = " << this->prm_[tau_j] << "\n";
    oss << "gamma_d = " << this->prm_[gamma_d] << "\n";
    oss << "E_Ca = " << this->prm_[E_Ca] << "\n";
    oss << "tau_d = " << this->prm_[tau_d] << "\n";
    oss << "xK1_inf = " << this->prm_[xK1_inf] << "\n";
    oss << "a = " << this->prm_[a] << "\n";
    oss << "kcasr = " << this->prm_[kcasr] << "\n";
    oss << "Ca_i_bufc = " << this->prm_[Ca_i_bufc] << "\n";
    oss << "k1 = " << this->prm_[k1] << "\n";
    oss << "k2 = " << this->prm_[k2] << "\n";
    oss << "O = " << this->prm_[O] << "\n";
    oss << "Ca_sr_bufsr = " << this->prm_[Ca_sr_bufsr] << "\n";
    oss << "Ca_ss_bufss = " << this->prm_[Ca_ss_bufss] << "\n";
    return oss.str();

}


std::string Stewart::PrintCurrents() const
{
    using namespace StrtCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "IfNa = " << this->cur_[IfNa] << "\n";
    oss << "IfK = " << this->cur_[IfK] << "\n";
    oss << "If = " << this->cur_[If] << "\n";
    oss << "Irel = " << this->cur_[Irel] << "\n";
    oss << "Isus = " << this->cur_[Isus] << "\n";
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "INaCa = " << this->cur_[INaCa] << "\n";
    oss << "IpCa = " << this->cur_[IpCa] << "\n";
    oss << "IpK = " << this->cur_[IpK] << "\n";
    oss << "Iup = " << this->cur_[Iup] << "\n";
    oss << "Ileak = " << this->cur_[Ileak] << "\n";
    oss << "Ixfer = " << this->cur_[Ixfer] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "IbNa = " << this->cur_[IbNa] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "IbCa = " << this->cur_[IbCa] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "Iion = " << this->cur_[StrtCur::Iion];
    return oss.str();

}


std::string Stewart::PrintBlockCoeffs() const
{
    using namespace StrtCur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "IfNa = " << this->block_coeff_[IfNa] << "\n";
    oss << "IfK = " << this->block_coeff_[IfK] << "\n";
    oss << "If = " << this->block_coeff_[If] << "\n";
    oss << "Irel = " << this->block_coeff_[Irel] << "\n";
    oss << "Isus = " << this->block_coeff_[Isus] << "\n";
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "INaCa = " << this->block_coeff_[INaCa] << "\n";
    oss << "IpCa = " << this->block_coeff_[IpCa] << "\n";
    oss << "IpK = " << this->block_coeff_[IpK] << "\n";
    oss << "Iup = " << this->block_coeff_[Iup] << "\n";
    oss << "Ileak = " << this->block_coeff_[Ileak] << "\n";
    oss << "Ixfer = " << this->block_coeff_[Ixfer] << "\n";
    oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
    oss << "IKr = " << this->block_coeff_[IKr] << "\n";
    oss << "IKs = " << this->block_coeff_[IKs] << "\n";
    oss << "INa = " << this->block_coeff_[INa] << "\n";
    oss << "IbNa = " << this->block_coeff_[IbNa] << "\n";
    oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
    oss << "IbCa = " << this->block_coeff_[IbCa] << "\n";
    oss << "Ito = " << this->block_coeff_[Ito];
    return oss.str();

}

} // End of namespace ELECTRA
