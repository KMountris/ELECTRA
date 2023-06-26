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


#include "ELECTRA/engine/electrophysiology/tentusscher2006.hpp"


namespace ELECTRA {

void TenTusscher2006::SetDataMapping()
{
    using namespace Tnt06Var;
    using namespace Tnt06Prm;
    using namespace Tnt06Cur;

    // Set variables mapping.
    this->mapped_data_["v"] = static_cast<std::size_t>(v);
    this->mapped_data_["ki"] = static_cast<std::size_t>(ki);
    this->mapped_data_["nai"] = static_cast<std::size_t>(nai);
    this->mapped_data_["cai"] = static_cast<std::size_t>(cai);
    this->mapped_data_["Xr1"] = static_cast<std::size_t>(Xr1);
    this->mapped_data_["Xr2"] = static_cast<std::size_t>(Xr2);
    this->mapped_data_["Xs"] = static_cast<std::size_t>(Xs);
    this->mapped_data_["m"] = static_cast<std::size_t>(m);
    this->mapped_data_["h"] = static_cast<std::size_t>(h);
    this->mapped_data_["j"] = static_cast<std::size_t>(j);
    this->mapped_data_["ca_ss"] = static_cast<std::size_t>(ca_ss);
    this->mapped_data_["d"] = static_cast<std::size_t>(d);
    this->mapped_data_["f"] = static_cast<std::size_t>(f);
    this->mapped_data_["f2"] = static_cast<std::size_t>(f2);
    this->mapped_data_["fCass"] = static_cast<std::size_t>(fCass);
    this->mapped_data_["s"] = static_cast<std::size_t>(s);
    this->mapped_data_["r"] = static_cast<std::size_t>(r);
    this->mapped_data_["ca_SR"] = static_cast<std::size_t>(ca_SR);
    this->mapped_data_["Rprime"] = static_cast<std::size_t>(Rprime);
  
    // Set parameters mapping.
    this->mapped_data_["R"] = static_cast<std::size_t>(R);
    this->mapped_data_["T"] = static_cast<std::size_t>(T);
    this->mapped_data_["F"] = static_cast<std::size_t>(F);
    this->mapped_data_["Cm"] = static_cast<std::size_t>(Cm);
    this->mapped_data_["Vc"] = static_cast<std::size_t>(Vc);
    this->mapped_data_["Pkna"] = static_cast<std::size_t>(Pkna);
    this->mapped_data_["ko"] = static_cast<std::size_t>(ko);
    this->mapped_data_["nao"] = static_cast<std::size_t>(nao);
    this->mapped_data_["cao"] = static_cast<std::size_t>(cao);
    this->mapped_data_["gK1"] = static_cast<std::size_t>(gK1);
    this->mapped_data_["gKr"] = static_cast<std::size_t>(gKr);
    this->mapped_data_["gKs"] = static_cast<std::size_t>(gKs);
    this->mapped_data_["gNa"] = static_cast<std::size_t>(gNa);
    this->mapped_data_["gNab"] = static_cast<std::size_t>(gNab);
    this->mapped_data_["gCaL"] = static_cast<std::size_t>(gCaL);
    this->mapped_data_["gCab"] = static_cast<std::size_t>(gCab);
    this->mapped_data_["gto"] = static_cast<std::size_t>(gto);
    this->mapped_data_["PNaK"] = static_cast<std::size_t>(PNaK);
    this->mapped_data_["Kmk"] = static_cast<std::size_t>(Kmk);
    this->mapped_data_["Kmna"] = static_cast<std::size_t>(Kmna);
    this->mapped_data_["Knaca"] = static_cast<std::size_t>(Knaca);
    this->mapped_data_["Ksat"] = static_cast<std::size_t>(Ksat);
    this->mapped_data_["alpha"] = static_cast<std::size_t>(alpha);
    this->mapped_data_["gamma"] = static_cast<std::size_t>(gamma);
    this->mapped_data_["Kmca"] = static_cast<std::size_t>(Kmca);
    this->mapped_data_["Kmnai"] = static_cast<std::size_t>(Kmnai);
    this->mapped_data_["gCap"] = static_cast<std::size_t>(gCap);
    this->mapped_data_["Kcap"] = static_cast<std::size_t>(Kcap);
    this->mapped_data_["gKp"] = static_cast<std::size_t>(gKp);
    this->mapped_data_["k1_prime"] = static_cast<std::size_t>(k1_prime);
    this->mapped_data_["k2_prime"] = static_cast<std::size_t>(k2_prime);
    this->mapped_data_["k3"] = static_cast<std::size_t>(k3);
    this->mapped_data_["k4"] = static_cast<std::size_t>(k4);
    this->mapped_data_["EC"] = static_cast<std::size_t>(EC);
    this->mapped_data_["max_sr"] = static_cast<std::size_t>(max_sr);
    this->mapped_data_["min_sr"] = static_cast<std::size_t>(min_sr);
    this->mapped_data_["Vrel"] = static_cast<std::size_t>(Vrel);
    this->mapped_data_["Vxfer"] = static_cast<std::size_t>(Vxfer);
    this->mapped_data_["Kup"] = static_cast<std::size_t>(Kup);
    this->mapped_data_["Vleak"] = static_cast<std::size_t>(Vleak);
    this->mapped_data_["Vmax_up"] = static_cast<std::size_t>(Vmax_up);
    this->mapped_data_["Buf_c"] = static_cast<std::size_t>(Buf_c);
    this->mapped_data_["K_buf_c"] = static_cast<std::size_t>(K_buf_c);
    this->mapped_data_["Buf_sr"] = static_cast<std::size_t>(Buf_sr);
    this->mapped_data_["K_buf_sr"] = static_cast<std::size_t>(K_buf_sr);
    this->mapped_data_["Buf_ss"] = static_cast<std::size_t>(Buf_ss);
    this->mapped_data_["K_buf_ss"] = static_cast<std::size_t>(K_buf_ss);
    this->mapped_data_["Vsr"] = static_cast<std::size_t>(Vsr);
    this->mapped_data_["Vss"] = static_cast<std::size_t>(Vss);

    // Set currents mapping.
    this->mapped_data_["INaK"] = static_cast<std::size_t>(INaK);
    this->mapped_data_["INa"] = static_cast<std::size_t>(INa);
    this->mapped_data_["INab"] = static_cast<std::size_t>(INab);
    this->mapped_data_["INaCa"] = static_cast<std::size_t>(INaCa);
    this->mapped_data_["IK1"] = static_cast<std::size_t>(IK1);
    this->mapped_data_["Ito"] = static_cast<std::size_t>(Ito);
    this->mapped_data_["IKr"] = static_cast<std::size_t>(IKr);
    this->mapped_data_["IKs"] = static_cast<std::size_t>(IKs);
    this->mapped_data_["ICaL"] = static_cast<std::size_t>(ICaL);
    this->mapped_data_["ICab"] = static_cast<std::size_t>(ICab);
    this->mapped_data_["IKp"] = static_cast<std::size_t>(IKp);
    this->mapped_data_["ICap"] = static_cast<std::size_t>(ICap);
 
}


TenTusscher2006::TenTusscher2006()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::TenTusscher2006;
    this->dt_stable_ = 0.02;
    this->var_.resize(20, 0.);
    this->prm_.resize(50, 0.);
    this->cur_.resize(13, 0.);
    this->block_coeff_.resize(12, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


TenTusscher2006::~TenTusscher2006()
{}


void TenTusscher2006::Initialize(CellType cell_type)
{
    using namespace Tnt06Var;
    using namespace Tnt06Prm;

    //Initialize the model's data.
    this->var_.clear();           this->var_.resize(20, 0.);
    this->prm_.clear();           this->prm_.resize(50, 0.);
    this->cur_.clear();           this->cur_.resize(13, 0.);
    this->block_coeff_.clear();   this->block_coeff_.resize(12, 0.);

    // Set data according to the cell type.
    switch (cell_type)
    {
        case CellType::endo :
            // Set the model variables for endocardium cell type.
            this->var_[v] = -86.709;
            this->var_[dvdt] = 0.;
            this->var_[ki] = 138.4;
            this->var_[nai] = 10.355;
            this->var_[cai] = 0.00013;
            this->var_[Xr1] = 0.00448;
            this->var_[Xr2] = 0.476;
            this->var_[Xs] = 0.0087;
            this->var_[m] = 0.00155;
            this->var_[h] = 0.7573;
            this->var_[j] = 0.7225; 
            this->var_[ca_ss] = 0.00036;
            this->var_[d] = 3.164e-5;
            this->var_[f] = 0.8009;
            this->var_[f2] = 0.9778;
            this->var_[fCass] = 0.9953;
            this->var_[s] = 0.3212;
            this->var_[r] = 2.235e-8;
            this->var_[ca_SR] = 3.715;
            this->var_[Rprime] = 0.9068;

            // Set the model parameters for endocardium cell type.
            this->prm_[gKs] = 0.392;
            this->prm_[gto] = 0.073;
            this->prm_[delta] = 0.;
            break;

        case CellType::mid :
            // Set the model variables for midmyocardium cell type.
            this->var_[v] = -85.423;
            this->var_[dvdt] = 0.;
            this->var_[ki] = 138.52;
            this->var_[nai] = 10.132;
            this->var_[cai] = 0.000153;
            this->var_[Xr1] = 0.0165;
            this->var_[Xr2] = 0.473;
            this->var_[Xs] = 0.0174;
            this->var_[m] = 0.00165;
            this->var_[h] = 0.749;
            this->var_[j] = 0.6788;
            this->var_[ca_ss] = 0.00042;
            this->var_[d] = 3.288e-5;
            this->var_[f] = 0.7026;
            this->var_[f2] = 0.9526;
            this->var_[fCass] = 0.9942;
            this->var_[s] = 0.999998;
            this->var_[r] = 2.347e-8;
            this->var_[ca_SR] = 4.272;
            this->var_[Rprime] = 0.8978;

            // Set the model parameters for midmyocardium cell type.
            this->prm_[gKs] = 0.098;
            this->prm_[gto] = 0.294;
            this->prm_[delta] = 1.;
            break;

        case CellType::epi : 
            // Set the model variables for epicardium cell type.
            this->var_[v] = -85.23;
            this->var_[dvdt] = 0.;
            this->var_[ki] = 136.89;      
            this->var_[nai] = 8.604;      
            this->var_[cai] = 0.000126;   
            this->var_[Xr1] = 0.00621;     
            this->var_[Xr2] = 0.4712;      
            this->var_[Xs] = 0.0095;       
            this->var_[m] = 0.00172;       
            this->var_[h] = 0.7444;        
            this->var_[j] = 0.7045;        
            this->var_[ca_ss] = 0.00036;   
            this->var_[d] = 3.373e-5;      
            this->var_[f] = 0.7888;        
            this->var_[f2] = 0.9755;       
            this->var_[fCass] = 0.9953;      
            this->var_[s] = 0.999998;      
            this->var_[r] = 2.42e-8;       
            this->var_[ca_SR] = 3.64;      
            this->var_[Rprime] = 0.9073;  
        
            // Set the model parameters for epicardium cell type.
            this->prm_[gKs] = 0.392;
            this->prm_[gto] = 0.294;
            this->prm_[delta] = 1.;
            break;
    
        default:
            throw std::invalid_argument(Logger::Error("Could not initialize TenTusscher2006 electrophysiology model. Expected: CellType::endo | CellType::mid | CellType::epi"));
            break;
    }

    // Set the cell model parameters that are independent from the cell type.
    this->prm_[R] = 8314.472;       
    this->prm_[T] = 310.;            
    this->prm_[F] = 96485.3415;     
    this->prm_[Cm] = 0.185;         
    this->prm_[Vc] = 0.016404;     
    this->prm_[Pkna] = 0.03;       
    this->prm_[ko] = 5.4;          
    this->prm_[nao] = 140.;         
    this->prm_[cao] = 2.;           
    this->prm_[gK1] = 5.405;       
    this->prm_[gKr] = 0.153;       
    this->prm_[gNa] = 14.838;      
    this->prm_[gNab] = 0.00029;    
    this->prm_[gCaL] = 0.0000398;  
    this->prm_[gCab] = 0.000592;   
    this->prm_[PNaK] = 2.724;      
    this->prm_[Kmk] = 1.;           
    this->prm_[Kmna] = 40.;         
    this->prm_[Knaca] = 1000.;      
    this->prm_[Ksat] = 0.1;        
    this->prm_[alpha] = 2.5;        
    this->prm_[gamma] = 0.35;       
    this->prm_[Kmca] = 1.38;       
    this->prm_[Kmnai] = 87.5;      
    this->prm_[gCap] = 0.1238;     
    this->prm_[Kcap] = 0.0005;     
    this->prm_[gKp] = 0.0146;      
    this->prm_[k1_prime] = 0.15;    
    this->prm_[k2_prime] = 0.045;   
    this->prm_[k3] = 0.06;          
    this->prm_[k4] = 0.005;         
    this->prm_[EC] = 1.5;           
    this->prm_[max_sr] = 2.5;       
    this->prm_[min_sr] = 1.;         
    this->prm_[Vrel] = 0.102;      
    this->prm_[Vxfer] = 0.0038;    
    this->prm_[Kup] = 0.00025;     
    this->prm_[Vleak] = 0.00036;   
    this->prm_[Vmax_up] = 0.006375; 
    this->prm_[Buf_c] = 0.2;        
    this->prm_[K_buf_c] = 0.001;    
    this->prm_[Buf_sr] = 10.;        
    this->prm_[K_buf_sr] = 0.3;     
    this->prm_[Buf_ss] = 0.4;       
    this->prm_[K_buf_ss] = 0.00025; 
    this->prm_[Vsr] = 0.001094;    
    this->prm_[Vss] = 0.00005468;  
}


void TenTusscher2006::Compute(double v_new, double dt, double stim_current)
{
    using namespace Tnt06Var;
    using namespace Tnt06Prm;
    using namespace Tnt06Cur;

    double f_inf = 1. / (1. + std::exp((v_new+20.)/7.));
    double tau_f = 1102.5 * std::exp( -std::pow(v_new+27., 2.) / 225.) + 200. / (1. + std::exp((13.0-v_new)/10.)) + 180. / (1.+std::exp((v_new+30.)/10.)) + 20.;
    this->var_[f] = ALGORITHM::RushLarsen(f_inf, this->var_[f], dt, tau_f);

    double f2_inf = 0.67 / (1.+std::exp((v_new+35.)/7.)) + 0.33;
    double tau_f2 = 562.*std::exp( -std::pow(v_new+27., 2.) / 240.) + 31. / (1. + std::exp((25.0-v_new)/10.)) + 80. / (1.+std::exp((v_new+30.)/10.));
    this->var_[f2] = ALGORITHM::RushLarsen(f2_inf, this->var_[f2], dt, tau_f2);

    double fCass_inf = 0.6 / (1.+std::pow(this->var_[ca_ss]/0.05, 2.)) + 0.4;
    double tau_fCass = 80. / (1.+std::pow(this->var_[ca_ss]/0.05, 2.)) + 2.;
    this->var_[fCass] = ALGORITHM::RushLarsen(fCass_inf, this->var_[fCass], dt, tau_fCass);

    double s_inf = 1. / (1.+std::exp((v_new+20.)/5.));
    double tau_s = 0.;
    if (static_cast<int>(this->prm_[delta]) == 0) 
        tau_s =  1000.*std::exp( -std::pow(v_new+67., 2.)/1000.) + 8.;
    else
        tau_s =  85.*std::exp( -std::pow(v_new+45., 2.)/320.) + 5./(1.+std::exp((v_new - 20.)/5.)) + 3.;
    this->var_[s] = ALGORITHM::RushLarsen(s_inf, this->var_[s], dt, tau_s);

    double r_inf = 1. / (1.+std::exp((20.0-v_new)/6.));
    double tau_r = 9.5*std::exp( -std::pow(v_new+40., 2.)/1800.) + 0.8;
    this->var_[r] = ALGORITHM::RushLarsen(r_inf, this->var_[r], dt, tau_r);

    double xr1_inf = 1. / (1.+std::exp(( -26.0-v_new)/7.));
    double alpha_xr1 = 450. / (1.+std::exp(( -45.0-v_new)/10.));
    double beta_xr1 = 6. / (1.+std::exp((v_new+30.)/11.5));
    double tau_xr1 = alpha_xr1 * beta_xr1;
    this->var_[Xr1] = ALGORITHM::RushLarsen(xr1_inf, this->var_[Xr1], dt, tau_xr1);

    double xr2_inf = 1. / (1.+std::exp((v_new+88.)/24.));
    double alpha_xr2 = 3. / (1.+std::exp((-60.0-v_new)/20.));
    double beta_xr2 = 1.12 / (1.+std::exp((v_new-60.0000)/20.));
    double tau_xr2 = alpha_xr2 * beta_xr2;
    this->var_[Xr2] = ALGORITHM::RushLarsen(xr2_inf, this->var_[Xr2], dt, tau_xr2);

    double xs_inf = 1. / (1.+std::exp((-5.0-v_new)/14.));
    double alpha_xs = 1400. / std::pow((1.+std::exp((5.0-v_new)/6.)), 0.5);
    double beta_xs = 1. / (1.+std::exp((v_new-35.)/15.));
    double tau_xs = alpha_xs * beta_xs + 80.;
    this->var_[Xs] = ALGORITHM::RushLarsen(xs_inf, this->var_[Xs], dt, tau_xs);
 
    double m_inf = 1. / std::pow(1.+std::exp((-56.86-v_new)/9.03), 2.);
    double alpha_m = 1. / (1.+std::exp((-60.0-v_new)/5.));
    double beta_m = 0.1 / (1.+std::exp((v_new+35.)/5.)) + 0.1/(1.+std::exp((v_new-50.)/200.));
    double tau_m = alpha_m * beta_m;
    this->var_[m] = ALGORITHM::RushLarsen(m_inf, this->var_[m], dt, tau_m);

    double h_inf = 1. / std::pow(1.+std::exp((v_new+71.55)/7.43), 2.);
    double alpha_h = 0.;
    double beta_h = 0.77 / (0.13 * (1.+std::exp((v_new+10.66) / -11.1)));
    if (v_new< - 40.) {
        alpha_h = 0.057 * std::exp(-(v_new+80.)/6.8);
        beta_h = 2.7 * std::exp(0.079*v_new) + 310000.*std::exp(0.3485*v_new);
    }
    double tau_h = 1. / (alpha_h+beta_h);
    this->var_[h] = ALGORITHM::RushLarsen(h_inf, this->var_[h], dt, tau_h);
    
    double j_inf = 1. / std::pow(1.+std::exp((v_new+71.55)/7.43), 2.);
    double alpha_j = 0.;
    double beta_j = (0.6*std::exp(0.057*v_new)) / (1.+std::exp(-0.1*(v_new+32.)));
    if (v_new < -40.) {
        alpha_j = (( (-25428.*std::exp(0.2444*v_new) - 6.948e-06*std::exp(-0.04391*v_new)) * (v_new+37.78))) / (1.+std::exp(0.311*(v_new+79.23)));
        beta_j = (0.02424*std::exp(-0.01052*v_new)) / (1.+std::exp(-0.1378*(v_new+40.14)));
    } 
    double tau_j = 1. / (alpha_j+beta_j);
    this->var_[j] = ALGORITHM::RushLarsen(j_inf, this->var_[j], dt, tau_j);
    
    double d_inf = 1. / (1.+std::exp((-8.0-v_new) / 7.5));
    double alpha_d = 1.4 / (1.+std::exp((-35.0-v_new) / 13.)) + 0.25;
    double beta_d = 1.4 / (1.+std::exp((v_new+5.) / 5.));
    double gamma_d = 1. / (1.+std::exp((50.0-v_new)/20.));
    double tau_d = alpha_d*beta_d + gamma_d;
    this->var_[d] = ALGORITHM::RushLarsen(d_inf, this->var_[d], dt, tau_d);
    
    this->cur_[INaK] = (( ((this->prm_[PNaK]*this->prm_[ko]) / (this->prm_[ko]+this->prm_[Kmk]))*this->var_[nai]) / (this->var_[nai]+this->prm_[Kmna])) / (1.+0.1245*std::exp((-0.1*v_new*this->prm_[F])/(this->prm_[R]*this->prm_[T])) + 0.0353*std::exp((-v_new*this->prm_[F])/(this->prm_[R]*this->prm_[T])));
    double E_Na = ((this->prm_[R]*this->prm_[T])/this->prm_[F]) * std::log(this->prm_[nao]/this->var_[nai]);
    this->cur_[INa] = this->prm_[gNa]*this->var_[m]*this->var_[m]*this->var_[m]*this->var_[h]*this->var_[j]*(v_new - E_Na);
    this->cur_[INab] =  this->prm_[gNab] * (v_new - E_Na);
    this->cur_[INaCa] = ( this->prm_[Knaca] * (std::exp((this->prm_[gamma]*v_new*this->prm_[F]) / (this->prm_[R]*this->prm_[T]))*std::pow(this->var_[nai], 3.) * this->prm_[cao] - std::exp(((this->prm_[gamma] - 1.)*v_new*this->prm_[F])/(this->prm_[R]*this->prm_[T]))*std::pow(this->prm_[nao], 3.)*this->var_[cai]*this->prm_[alpha]))/((std::pow(this->prm_[Kmnai], 3.)+std::pow(this->prm_[nao], 3.))*(this->prm_[Kmca]+this->prm_[cao])*(1.+this->prm_[Ksat]*std::exp(((this->prm_[gamma] - 1.)*v_new*this->prm_[F])/(this->prm_[R]*this->prm_[T]))));
    double dNa_i = (( -(this->cur_[INa]+this->cur_[INab] + 3.*this->cur_[INaK] + 3.*this->cur_[INaCa])) / (this->prm_[Vc]*this->prm_[F]))*this->prm_[Cm];
    this->var_[nai] = ALGORITHM::ForwardEuler(this->var_[nai], dt, dNa_i);

    double E_K = ((this->prm_[R]*this->prm_[T])/this->prm_[F]) * std::log(this->prm_[ko]/this->var_[ki]);
    double alpha_K1 = 0.1 / (1.+std::exp(0.06*((v_new-E_K) - 200.)));
    double beta_K1 = (3.*std::exp(0.0002*((v_new-E_K)+100.)) + std::exp(0.1*((v_new - E_K) - 10.))) / (1.+std::exp(-0.5*(v_new - E_K)));
    double xK1_inf = alpha_K1 / (alpha_K1+beta_K1);
    this->cur_[IK1] = this->prm_[gK1] * xK1_inf*std::pow((this->prm_[ko]/5.4), 0.5) * (v_new-E_K);
    this->cur_[Ito] = this->prm_[gto] * this->var_[r]*this->var_[s]*(v_new - E_K);
    this->cur_[IKr] = this->prm_[gKr] * std::pow((this->prm_[ko]/5.4), 0.5) * this->var_[Xr1]*this->var_[Xr2]*(v_new - E_K);

    double E_Ks = ((this->prm_[R]*this->prm_[T])/this->prm_[F])*std::log((this->prm_[ko] + this->prm_[Pkna]*this->prm_[nao])/(this->var_[ki] + this->prm_[Pkna]*this->var_[nai]));
    this->cur_[IKs] = this->prm_[gKs]*std::pow(this->var_[Xs], 2.) * (v_new - E_Ks);
    this->cur_[ICaL] = (((this->prm_[gCaL]*this->var_[d]*this->var_[f]*this->var_[f2]*this->var_[fCass]*4.*(v_new-15.) * std::pow(this->prm_[F], 2.)) / (this->prm_[R]*this->prm_[T]))*(0.25*this->var_[ca_ss]*std::exp((2.*(v_new-15.)*this->prm_[F])/(this->prm_[R]*this->prm_[T])) - this->prm_[cao])) / (std::exp((2.*(v_new-15.)*this->prm_[F]) / (this->prm_[R]*this->prm_[T])) - 1.);

    double E_Ca = ((0.5*this->prm_[R]*this->prm_[T])/this->prm_[F])*std::log(this->prm_[cao]/this->var_[cai]);
    this->cur_[ICab] =  this->prm_[gCab]*(v_new - E_Ca);
    this->cur_[IKp] = (this->prm_[gKp]*(v_new - E_K)) / (1.+std::exp((25.0-v_new) / 5.98));
    this->cur_[ICap] = (this->prm_[gCap]*this->var_[cai]) / (this->var_[cai]+this->prm_[Kcap]); 

    // Compute the total Iion current.
    this->cur_[Tnt06Cur::Iion] = this->cur_[IK1] + this->cur_[Ito] + this->cur_[IKr] + this->cur_[IKs] + this->cur_[ICaL] + 
                            this->cur_[INaK] + this->cur_[INa] + this->cur_[INab] + this->cur_[INaCa] + this->cur_[ICab] + 
                            this->cur_[IKp] +this->cur_[ICap];

    // Set the new value to the membrante potential time derivative.
    this->var_[dvdt] = - (this->cur_[Tnt06Cur::Iion] - stim_current);    

    double dK_i = (( -((this->cur_[IK1]+this->cur_[Ito]+this->cur_[IKr]+this->cur_[IKs]+this->cur_[IKp]-stim_current) - 2.*this->cur_[INaK])) / (this->prm_[Vc]*this->prm_[F]))*this->prm_[Cm];
    this->var_[ki] = ALGORITHM::ForwardEuler(this->var_[ki], dt, dK_i);

    double i_up = this->prm_[Vmax_up] / (1.+std::pow(this->prm_[Kup], 2.) / std::pow(this->var_[cai], 2.));
    double i_leak = this->prm_[Vleak] * (this->var_[ca_SR] - this->var_[cai]);
    double i_xfer = this->prm_[Vxfer] * (this->var_[ca_ss] - this->var_[cai]);
    double Ca_i_bufc = 1. / (1.+(this->prm_[Buf_c]*this->prm_[K_buf_c]) / std::pow(this->var_[cai]+this->prm_[K_buf_c], 2.));
    double dCa_i = Ca_i_bufc*((((i_leak - i_up)*this->prm_[Vsr])/this->prm_[Vc]+i_xfer) - (((this->cur_[ICab]+this->cur_[ICap]) - 2.*this->cur_[INaCa])*this->prm_[Cm]) / (2.*this->prm_[Vc]*this->prm_[F]));
    this->var_[cai] = ALGORITHM::ForwardEuler(this->var_[cai], dt, dCa_i);

    double kcasr = this->prm_[max_sr] - (this->prm_[max_sr] - this->prm_[min_sr]) / (1.+std::pow(this->prm_[EC]/this->var_[ca_SR], 2.));
    double k2 =  this->prm_[k2_prime]*kcasr;
    double dR_prime = -k2*this->var_[ca_ss]*this->var_[Rprime] + this->prm_[k4]*(1.0-this->var_[Rprime]);
    this->var_[Rprime] = ALGORITHM::ForwardEuler(this->var_[Rprime], dt, dR_prime);

    double k1 = this->prm_[k1_prime] / kcasr;
    double O = (k1*std::pow(this->var_[ca_ss], 2.)*this->var_[Rprime]) / (this->prm_[k3] + k1*std::pow(this->var_[ca_ss], 2.));
    double i_rel =  this->prm_[Vrel] * O * (this->var_[ca_SR] - this->var_[ca_ss]);
    double Ca_sr_bufsr = 1. / (1.+(this->prm_[Buf_sr]*this->prm_[K_buf_sr]) / std::pow(this->var_[ca_SR]+this->prm_[K_buf_sr], 2.));
    double dCa_SR = Ca_sr_bufsr * (i_up - (i_rel+i_leak));
    this->var_[ca_SR] = ALGORITHM::ForwardEuler(this->var_[ca_SR], dt, dCa_SR);

    double Ca_ss_bufss = 1. / (1.+(this->prm_[Buf_ss]*this->prm_[K_buf_ss]) / std::pow(this->var_[ca_ss]+this->prm_[K_buf_ss], 2.));
    double dCa_ss = Ca_ss_bufss*(((-this->cur_[ICaL]*this->prm_[Cm]) / (2.*this->prm_[Vss]*this->prm_[F]) + (i_rel*this->prm_[Vsr]) / this->prm_[Vss]) - (i_xfer*this->prm_[Vc])/this->prm_[Vss]);
    this->var_[ca_ss] = ALGORITHM::ForwardEuler(this->var_[ca_ss], dt, dCa_ss);   
}


std::string TenTusscher2006::PrintVariables() const 
{
    using namespace Tnt06Var;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "ki = " << this->var_[ki] << "\n";
    oss << "nai = " << this->var_[nai] << "\n";
    oss << "cai = " << this->var_[cai] << "\n";
    oss << "Xr1 = " << this->var_[Xr1] << "\n";
    oss << "Xr2 = " << this->var_[Xr2] << "\n";
    oss << "Xs = " << this->var_[Xs] << "\n";
    oss << "m = " << this->var_[m] << "\n";
    oss << "h = " << this->var_[h] << "\n";
    oss << "j = " << this->var_[j] << "\n";
    oss << "ca_ss = " << this->var_[ca_ss] << "\n";
    oss << "d = " << this->var_[d] << "\n";
    oss << "f = " << this->var_[f] << "\n";
    oss << "f2 = " << this->var_[f2] << "\n";
    oss << "fCass = " << this->var_[fCass] << "\n";
    oss << "s = " << this->var_[s] << "\n";
    oss << "r = " << this->var_[r] << "\n";
    oss << "ca_SR = " << this->var_[ca_SR] << "\n";
    oss << "Rprime = " << this->var_[Rprime];
    return oss.str();
}


std::string TenTusscher2006::PrintParameters() const 
{
    using namespace Tnt06Prm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "R = " << this->prm_[R] << "\n";
    oss << "T = " << this->prm_[T] << "\n";
    oss << "F = " << this->prm_[F] << "\n";
    oss << "Cm = " << this->prm_[Cm] << "\n";
    oss << "Vc = " << this->prm_[Vc] << "\n";
    oss << "Pkna = " << this->prm_[Pkna] << "\n";
    oss << "ko = " << this->prm_[ko] << "\n";
    oss << "nao = " << this->prm_[nao] << "\n";
    oss << "cao = " << this->prm_[cao] << "\n";
    oss << "gK1 = " << this->prm_[gK1] << "\n";
    oss << "gKr = " << this->prm_[gKr] << "\n";
    oss << "gKs = " << this->prm_[gKs] << "\n";
    oss << "gNa = " << this->prm_[gNa] << "\n";
    oss << "gNab = " << this->prm_[gNab] << "\n";
    oss << "gCaL = " << this->prm_[gCaL] << "\n";
    oss << "gCab = " << this->prm_[gCab] << "\n";
    oss << "gto = " << this->prm_[gto] << "\n";
    oss << "PNaK = " << this->prm_[PNaK] << "\n";
    oss << "Kmk = " << this->prm_[Kmk] << "\n";
    oss << "Kmna = " << this->prm_[Kmna] << "\n";
    oss << "Knaca = " << this->prm_[Knaca] << "\n";
    oss << "Ksat = " << this->prm_[Ksat] << "\n";
    oss << "alpha = " << this->prm_[alpha] << "\n";
    oss << "gamma = " << this->prm_[gamma] << "\n";
    oss << "Kmca = " << this->prm_[Kmca] << "\n";
    oss << "Kmnai = " << this->prm_[Kmnai] << "\n";
    oss << "gCap = " << this->prm_[gCap] << "\n";
    oss << "Kcap = " << this->prm_[Kcap] << "\n";
    oss << "gKp = " << this->prm_[gKp] << "\n";
    oss << "k1_prime = " << this->prm_[k1_prime] << "\n";
    oss << "k2_prime = " << this->prm_[k2_prime] << "\n";
    oss << "k3 = " << this->prm_[k3] << "\n";
    oss << "k4 = " << this->prm_[k4] << "\n";
    oss << "EC = " << this->prm_[EC] << "\n";
    oss << "max_sr = " << this->prm_[max_sr] << "\n";
    oss << "min_sr = " << this->prm_[min_sr] << "\n";
    oss << "Vrel = " << this->prm_[Vrel] << "\n";
    oss << "Vxfer = " << this->prm_[Vxfer] << "\n";
    oss << "Kup = " << this->prm_[Kup] << "\n";
    oss << "Vleak = " << this->prm_[Vleak] << "\n";
    oss << "Vmax_up = " << this->prm_[Vmax_up] << "\n";
    oss << "Buf_c = " << this->prm_[Buf_c] << "\n";
    oss << "K_buf_c = " << this->prm_[K_buf_c] << "\n";
    oss << "Buf_sr = " << this->prm_[Buf_sr] << "\n";
    oss << "K_buf_sr = " << this->prm_[K_buf_sr] << "\n";
    oss << "Buf_ss = " << this->prm_[Buf_ss] << "\n";
    oss << "K_buf_ss = " << this->prm_[K_buf_ss] << "\n";
    oss << "Vsr = " << this->prm_[Vsr] << "\n";
    oss << "Vss = " << this->prm_[Vss] << "\n";
    oss << "delta = " << this->prm_[delta];
    return oss.str();
}


std::string TenTusscher2006::PrintCurrents() const
{
    using namespace Tnt06Cur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INaK = " << this->cur_[INaK] << "\n";
    oss << "INa = " << this->cur_[INa] << "\n";
    oss << "INab = " << this->cur_[INab] << "\n";
    oss << "INaCa = " << this->cur_[INaCa] << "\n";
    oss << "IK1 = " << this->cur_[IK1] << "\n";
    oss << "Ito = " << this->cur_[Ito] << "\n";
    oss << "IKr = " << this->cur_[IKr] << "\n";
    oss << "IKs = " << this->cur_[IKs] << "\n";
    oss << "ICaL = " << this->cur_[ICaL] << "\n";
    oss << "ICab = " << this->cur_[ICab] << "\n";
    oss << "IKp = " << this->cur_[IKp] << "\n";
    oss << "ICap = " << this->cur_[ICap] << "\n";
    oss << "Iion = " << this->cur_[Tnt06Cur::Iion];
    return oss.str();

}


std::string TenTusscher2006::PrintBlockCoeffs() const
{
    using namespace Tnt06Cur;

    // Create output string stream to pass the currents and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "INaK = " << this->block_coeff_[INaK] << "\n";
    oss << "INa = " << this->block_coeff_[INa] << "\n";
    oss << "INab = " << this->block_coeff_[INab] << "\n";
    oss << "INaCa = " << this->block_coeff_[INaCa] << "\n";
    oss << "IK1 = " << this->block_coeff_[IK1] << "\n";
    oss << "Ito = " << this->block_coeff_[Ito] << "\n";
    oss << "IKr = " << this->block_coeff_[IKr] << "\n";
    oss << "IKs = " << this->block_coeff_[IKs] << "\n";
    oss << "ICaL = " << this->block_coeff_[ICaL] << "\n";
    oss << "ICab = " << this->block_coeff_[ICab] << "\n";
    oss << "IKp = " << this->block_coeff_[IKp] << "\n";
    oss << "ICap = " << this->block_coeff_[ICap];
    return oss.str();
}

} // End of namespace ELECTRA
