/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/electrophysiology/bueno.hpp"


namespace ELECTRA {


void Bueno::SetDataMapping()
{
    using namespace BueVar;
    using namespace BuePrm;
    using namespace BueCur;

    // Set variables mapping.
    this->mapped_data_["v"]           = static_cast<std::size_t>(v);
    this->mapped_data_["Vm"]          = static_cast<std::size_t>(Vm);
    this->mapped_data_["g_u"]         = static_cast<std::size_t>(g_u);
    this->mapped_data_["g_v"]         = static_cast<std::size_t>(g_v);
    this->mapped_data_["g_w"]         = static_cast<std::size_t>(g_w);
    this->mapped_data_["g_s"]         = static_cast<std::size_t>(g_s);
    this->mapped_data_["p"]           = static_cast<std::size_t>(p);
    this->mapped_data_["tau_s"]       = static_cast<std::size_t>(tau_s);
    this->mapped_data_["m"]           = static_cast<std::size_t>(m);
    this->mapped_data_["v_inf"]       = static_cast<std::size_t>(v_inf);
    this->mapped_data_["q"]           = static_cast<std::size_t>(q);
    this->mapped_data_["tau_v_minus"] = static_cast<std::size_t>(tau_v_minus);
    this->mapped_data_["r"]           = static_cast<std::size_t>(r);
    this->mapped_data_["w_inf"]       = static_cast<std::size_t>(w_inf);
    this->mapped_data_["tau_w_minus"] = static_cast<std::size_t>(tau_w_minus);
    this->mapped_data_["tau_o"]       = static_cast<std::size_t>(tau_o);
    this->mapped_data_["tau_so"]      = static_cast<std::size_t>(tau_so);
    

    // Set parameters mapping.
    this->mapped_data_["V0"]           = static_cast<std::size_t>(V0);
    this->mapped_data_["Vfi"]          = static_cast<std::size_t>(Vfi);
    this->mapped_data_["u_m"]          = static_cast<std::size_t>(u_m);
    this->mapped_data_["u_p"]          = static_cast<std::size_t>(u_p);
    this->mapped_data_["tau_v_plus"]   = static_cast<std::size_t>(tau_v_plus);
    this->mapped_data_["tau_s1"]       = static_cast<std::size_t>(tau_s1);
    this->mapped_data_["k_s"]          = static_cast<std::size_t>(k_s);
    this->mapped_data_["u_s"]          = static_cast<std::size_t>(u_s);
    this->mapped_data_["tau_s2"]       = static_cast<std::size_t>(tau_s2);
    this->mapped_data_["tau_v2_minus"] = static_cast<std::size_t>(tau_v2_minus);
    this->mapped_data_["u_q"]          = static_cast<std::size_t>(u_q);
    this->mapped_data_["tau_o1"]       = static_cast<std::size_t>(tau_o1);
    this->mapped_data_["u_r"]          = static_cast<std::size_t>(u_r);
    this->mapped_data_["tau_fi"]       = static_cast<std::size_t>(tau_fi);
    this->mapped_data_["u_u"]          = static_cast<std::size_t>(u_u);
    this->mapped_data_["tau_so1"]      = static_cast<std::size_t>(tau_so1);
    this->mapped_data_["tau_so2"]      = static_cast<std::size_t>(tau_so2);
    this->mapped_data_["k_so"]         = static_cast<std::size_t>(k_so);
    this->mapped_data_["u_so"]         = static_cast<std::size_t>(u_so);
    this->mapped_data_["tau_si"]       = static_cast<std::size_t>(tau_si);
    this->mapped_data_["tau_winf"]     = static_cast<std::size_t>(tau_winf);
    this->mapped_data_["wstar_inf"]    = static_cast<std::size_t>(wstar_inf);
    this->mapped_data_["tau_w1_minus"] = static_cast<std::size_t>(tau_w1_minus);
    this->mapped_data_["tau_w2_minus"] = static_cast<std::size_t>(tau_w2_minus);
    this->mapped_data_["k_w_minus"]    = static_cast<std::size_t>(k_w_minus);
    this->mapped_data_["u_w_minus"]    = static_cast<std::size_t>(u_w_minus);
    this->mapped_data_["tau_w_plus"]   = static_cast<std::size_t>(tau_w_plus);
    this->mapped_data_["tau_v1_minus"] = static_cast<std::size_t>(tau_v1_minus);
    this->mapped_data_["tau_o2"]       = static_cast<std::size_t>(tau_o2);    

    // Set currents mapping.
    this->mapped_data_["Ifi"]  = static_cast<std::size_t>(Ifi);
    this->mapped_data_["Iso"]  = static_cast<std::size_t>(Iso);
    this->mapped_data_["Isi"]  = static_cast<std::size_t>(Isi);
    
}


Bueno::Bueno()
{
    //Initialize the model's data.
    this->model_type_ = EpModelType::Bueno;
    this->dt_stable_ = 0.02;
    this->var_.resize(18, 0.);
    this->prm_.resize(29, 0.);
    this->cur_.resize(4, 0.);
    this->block_coeff_.resize(3, 0.);

    // Set mapped data.
    this->SetDataMapping();
}


Bueno::~Bueno()
{}


void Bueno::Initialize(CellType cell_type)
{
    using namespace BueVar;
    using namespace BuePrm;

    //Initialize the model's data.
    this->var_.clear();    this->var_.resize(18, 0.);
    this->prm_.clear();  this->prm_.resize(29, 0.);
    this->cur_.clear();    this->cur_.resize(4, 0.);
    this->block_coeff_.clear();  this->block_coeff_.resize(3, 0.);


    // Set the cell model parameters.
    // Same for each cell type.
    this->prm_[V0]  = -84.;
    this->prm_[Vfi] = 2.7;
    this->prm_[u_m] = 0.3;
    this->prm_[u_p] = 0.13;
    this->prm_[tau_v_plus] = 1.45;
    this->prm_[tau_s1]     = 2.7342;
    this->prm_[k_s] = 2.0994;
    this->prm_[u_s] = 0.9087;

    // Different for each cell type.
    switch (cell_type)
    {
        case CellType::endo :
            // Set the cell model parameters for endocardium cell type.
            this->prm_[tau_s2] = 2.;
            this->prm_[tau_v2_minus] = 10.;
            this->prm_[u_q] = 0.024;
            this->prm_[tau_o1] = 470.;
            this->prm_[u_r] = 0.006;
            this->prm_[tau_fi] = 0.104;
            this->prm_[u_u] = 1.56;
            this->prm_[tau_so1] = 40.;
            this->prm_[tau_so2] = 1.2;
            this->prm_[k_so] = 2.;
            this->prm_[u_so] = 0.65;
            this->prm_[tau_si] = 2.9013;
            this->prm_[tau_winf] = 0.0273;
            this->prm_[wstar_inf] = 0.78;
            this->prm_[tau_w1_minus] = 6.;
            this->prm_[tau_w2_minus] = 140;
            this->prm_[k_w_minus] = 200.;
            this->prm_[u_w_minus] = 0.016;
            this->prm_[tau_w_plus] = 280.;
            this->prm_[tau_v1_minus] = 75.;
            this->prm_[tau_o2] = 6.; 

            break;

        case CellType::mid :

            // Set the cell model parameters for midmyocardium cell type.
            this->prm_[tau_s2] = 4.;
            this->prm_[tau_v2_minus] = 1.45;
            this->prm_[u_q] = 0.1;
            this->prm_[tau_o1] = 410.;
            this->prm_[u_r] = 0.005;
            this->prm_[tau_fi] = 0.078;
            this->prm_[u_u] = 1.61;
            this->prm_[tau_so1] = 91.;
            this->prm_[tau_so2] = 0.8;
            this->prm_[k_so] = 2.1;
            this->prm_[u_so] = 0.6;
            this->prm_[tau_si] = 3.3849;
            this->prm_[tau_winf] = 0.01;
            this->prm_[wstar_inf] = 0.5;
            this->prm_[tau_w1_minus] = 70.;
            this->prm_[tau_w2_minus] = 8.;
            this->prm_[k_w_minus] = 200.;
            this->prm_[u_w_minus] = 0.016;
            this->prm_[tau_w_plus] = 280.;
            this->prm_[tau_v1_minus] = 80.;
            this->prm_[tau_o2] = 7.;

            break;

        case CellType::epi : 

            // Set the cell model parameters for epicardium cell type.
            this->prm_[tau_s2] = 16.;
            this->prm_[tau_v2_minus] = 1150.;
            this->prm_[u_q] = 0.006;
            this->prm_[tau_o1] = 400.;
            this->prm_[u_r] = 0.006;
            this->prm_[tau_fi] = 0.11;
            this->prm_[u_u] = 1.55;
            this->prm_[tau_so1] = 30.02;
            this->prm_[tau_so2] = 0.996;
            this->prm_[k_so] = 2.04;
            this->prm_[u_so] = 0.65;
            this->prm_[tau_si] = 1.8875;
            this->prm_[tau_winf] = 0.07;
            this->prm_[wstar_inf] = 0.94;
            this->prm_[tau_w1_minus] = 60.;
            this->prm_[tau_w2_minus] = 15.;
            this->prm_[k_w_minus] = 65.;
            this->prm_[u_w_minus] = 0.03;
            this->prm_[tau_w_plus] = 200.;
            this->prm_[tau_v1_minus] = 60.;
            this->prm_[tau_o2] = 6.;

            break;
    
        default:
            throw std::invalid_argument(Logger::Error("Could not initialize Bueno ap model. Not supported cell type. Expected: [CellType::endo | CellType::mid | CellType::epi]"));
            break;
    }

    // Set the cell variables.
    this->var_[v] = this->prm_[V0];
    this->var_[g_u] = 0.;
    this->var_[g_v] = 1.;
    this->var_[g_w] = 1.;
    this->var_[g_s] = 0.;
    this->var_[p] = (this->var_[g_u] < this->prm_[u_p]) ? 0. : 1.; 
    this->var_[tau_s] =  (1. - this->var_[p])*this->prm_[tau_s1] + this->var_[p]*this->prm_[tau_s2];
    this->var_[m] = (this->var_[g_u] < this->prm_[u_m]) ? 0. : 1.;
    this->var_[v_inf] = (this->var_[g_u] < this->prm_[u_q]) ? 1. : 0.;
    this->var_[q] = (this->var_[g_u] < this->prm_[u_q]) ? 0. : 1.;
    this->var_[tau_v_minus] =  this->var_[q]*this->prm_[tau_v2_minus] + (1. - this->var_[q])*this->prm_[tau_v1_minus];
    this->var_[r] = (this->var_[g_u] < this->prm_[u_r]) ? 0. : 1.;
    this->var_[w_inf] = (1. - this->var_[r]) * (1. - this->var_[g_u]/this->prm_[tau_winf]) + this->var_[r]*this->prm_[wstar_inf];
    this->var_[tau_w_minus] = this->prm_[tau_w1_minus] + 0.5*((this->prm_[tau_w2_minus] - this->prm_[tau_w1_minus]) * (1. + std::tanh(this->prm_[k_w_minus] * (this->var_[g_u] - this->prm_[u_w_minus]))));
    this->var_[tau_o] =  (1. - this->var_[r]) * this->prm_[tau_o1] + this->var_[r] * this->prm_[tau_o2];
    this->var_[tau_so] = this->prm_[tau_so1] + 0.5*((this->prm_[tau_so2] - this->prm_[tau_so1]) * (1. + std::tanh(this->prm_[k_so]*(this->var_[g_u] - this->prm_[u_so]))));
    this->var_[Vm] = this->prm_[V0] + this->var_[g_u] * (this->prm_[Vfi] - this->prm_[V0]);

}


void Bueno::Compute(double v_new, double dt, double stim_current)
{
    using namespace BueVar;
    using namespace BuePrm;
    using namespace BueCur;

    // Convert new gate u value to dimensionless.
    this->var_[g_u] = (v_new-this->prm_[V0]) / (this->prm_[Vfi] - this->prm_[V0]); 

    // Compute slow inward current gate s derivative.
    this->var_[p] = (this->var_[g_u] < this->prm_[u_p]) ? 0. : 1.;
    this->var_[tau_s] = (1. -this->var_[p]) * this->prm_[tau_s1] + this->var_[p]*this->prm_[tau_s2];
    double ds = (0.5*(1. + std::tanh(this->prm_[k_s]*(this->var_[g_u] - this->prm_[u_s]))) - this->var_[g_s]) / this->var_[tau_s];
    this->var_[g_s] = ALGORITHM::ForwardEuler(this->var_[g_s], dt, ds);

    // Compute fast inward current gate v derivative.
    this->var_[m] = (this->var_[g_u] < this->prm_[u_m]) ? 0. : 1.;
    this->var_[v_inf] = (this->var_[g_u] < this->prm_[u_q]) ? 1. : 0.;
    this->var_[q] = (this->var_[g_u] < this->prm_[u_q]) ? 0. : 1.;
    this->var_[tau_v_minus] = this->var_[q]*this->prm_[tau_v2_minus] + (1. - this->var_[q])*this->prm_[tau_v1_minus];
    double dv = ((1. - this->var_[m])*(this->var_[v_inf] - this->var_[g_v])) / this->var_[tau_v_minus] - (this->var_[m]*this->var_[g_v])/this->prm_[tau_v_plus];
    this->var_[g_v] = ALGORITHM::ForwardEuler(this->var_[g_v], dt, dv);

    // Compute slow inward current gate w derivative.
    this->var_[r] = (this->var_[g_u] < this->prm_[u_r]) ? 0. : 1.;
    this->var_[w_inf] = (1. - this->var_[r]) * (1. - this->var_[g_u]/this->prm_[tau_winf]) + this->var_[r]*this->prm_[wstar_inf];
    this->var_[tau_w_minus] = this->prm_[tau_w1_minus] + 0.5*((this->prm_[tau_w2_minus] - this->prm_[tau_w1_minus]) * (1. + std::tanh(this->prm_[k_w_minus]*(this->var_[g_u] - this->prm_[u_w_minus]))) );
    double dw = ((1. - this->var_[r]) * (this->var_[w_inf] - this->var_[g_w])) / this->var_[tau_w_minus] - (this->var_[r] * this->var_[g_w]) / this->prm_[tau_w_plus];
    this->var_[g_w] = ALGORITHM::ForwardEuler(this->var_[g_w], dt, dw);

    // Compute fast inward current.
    this->cur_[Ifi] = (1.0-this->block_coeff_[Ifi]) * ((- this->var_[m] * this->var_[g_v] * (this->var_[g_u] - this->prm_[u_m]) * (this->prm_[u_u] - this->var_[g_u])) / this->prm_[tau_fi]);

    // Compute slow outward current.
    this->var_[tau_o]  = (1. - this->var_[r])*this->prm_[tau_o1] + this->var_[r]*this->prm_[tau_o2];
    this->var_[tau_so] = this->prm_[tau_so1] + 0.5*((this->prm_[tau_so2] - this->prm_[tau_so1]) * (1. + std::tanh(this->prm_[k_so]*(this->var_[g_u] - this->prm_[u_so]))));
    this->cur_[Iso] = (1. - this->block_coeff_[Iso]) * ((this->var_[g_u]*(1. - this->var_[p]))/this->var_[tau_o] + this->var_[p]/this->var_[tau_so]);

    // Compute slow inward current.
    this->cur_[Isi] = (1. - this->block_coeff_[Isi]) * ((-this->var_[p]*this->var_[g_w]*this->var_[g_s]) / this->prm_[tau_si]);

    // Compute total ionic current.
    this->cur_[BueCur::Iion] = -(this->cur_[Ifi] + this->cur_[Iso] + this->cur_[Isi]);

    // Compute membrane potential gate u.
    this->var_[g_u] = ALGORITHM::ForwardEuler(this->var_[g_u], dt, (this->cur_[BueCur::Iion] + stim_current));
    
    this->var_[Vm] = this->prm_[V0] + this->var_[g_u]*(this->prm_[Vfi] - this->prm_[V0]);  

    this->var_[dvdt] = (this->var_[Vm] - v_new) / dt;
}


std::string Bueno::PrintVariables() const 
{
    using namespace BueVar;

    // Create output string stream to pass the variables and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "v = " << this->var_[v] << "\n";      
    oss << "dvdt = " << this->var_[dvdt] << "\n";
    oss << "Vm = " << this->var_[Vm] << "\n";
    oss << "g_u = " << this->var_[g_u] << "\n";
    oss << "g_v = " << this->var_[g_v] << "\n";
    oss << "g_w = " << this->var_[g_w] << "\n";
    oss << "g_s = " << this->var_[g_s] << "\n";
    oss << "p = " << this->var_[p] << "\n";
    oss << "tau_s = " << this->var_[tau_s] << "\n";
    oss << "m = " << this->var_[m] << "\n";
    oss << "v_inf = " << this->var_[v_inf] << "\n";
    oss << "q = " << this->var_[q] << "\n";
    oss << "tau_v_minus = " << this->var_[tau_v_minus] << "\n";
    oss << "r = " << this->var_[r] << "\n";
    oss << "w_inf = " << this->var_[w_inf] << "\n";
    oss << "tau_w_minus = " << this->var_[tau_w_minus] << "\n";
    oss << "tau_o = " << this->var_[tau_o] << "\n";
    oss << "tau_so = " << this->var_[tau_so];
    return oss.str();

}


std::string Bueno::PrintParameters() const 
{
    using namespace BuePrm;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "V0 = " << this->prm_[V0] << "\n";
    oss << "Vfi = " << this->prm_[Vfi] << "\n";
    oss << "u_m = " << this->prm_[u_m] << "\n";
    oss << "u_p = " << this->prm_[u_p] << "\n";
    oss << "tau_v_plus = " << this->prm_[tau_v_plus] << "\n";
    oss << "tau_s1 = " << this->prm_[tau_s1] << "\n";
    oss << "k_s = " << this->prm_[k_s] << "\n";
    oss << "u_s = " << this->prm_[u_s] << "\n";
    oss << "tau_s2 = " << this->prm_[tau_s2] << "\n";
    oss << "tau_v2_minus = " << this->prm_[tau_v2_minus] << "\n";
    oss << "u_q = " << this->prm_[u_q] << "\n";
    oss << "tau_o1 = " << this->prm_[tau_o1] << "\n";
    oss << "u_r = " << this->prm_[u_r] << "\n";
    oss << "tau_fi = " << this->prm_[tau_fi] << "\n";
    oss << "u_u = " << this->prm_[u_u] << "\n";
    oss << "tau_so1 = " << this->prm_[tau_so1] << "\n";
    oss << "tau_so2 = " << this->prm_[tau_so2] << "\n";
    oss << "k_so = " << this->prm_[k_so] << "\n";
    oss << "u_so = " << this->prm_[u_so] << "\n";
    oss << "tau_si = " << this->prm_[tau_si] << "\n";
    oss << "tau_winf = " << this->prm_[tau_winf] << "\n";
    oss << "wstar_inf = " << this->prm_[wstar_inf] << "\n";
    oss << "tau_w1_minus = " << this->prm_[tau_w1_minus] << "\n";
    oss << "tau_w2_minus = " << this->prm_[tau_w2_minus] << "\n";
    oss << "k_w_minus = " << this->prm_[k_w_minus] << "\n";
    oss << "u_w_minus = " << this->prm_[u_w_minus] << "\n";
    oss << "tau_w_plus = " << this->prm_[tau_w_plus] << "\n";
    oss << "tau_v1_minus = " << this->prm_[tau_v1_minus] << "\n";
    oss << "tau_o2 = " << this->prm_[tau_o2];
    return oss.str();

}


std::string Bueno::PrintCurrents() const
{
    using namespace BueCur;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "Ifi = " << this->cur_[Ifi] << "\n";
    oss << "Iso = " << this->cur_[Iso] << "\n";
    oss << "Isi = " << this->cur_[Isi] << "\n";
    oss << "Iion = " << this->cur_[BueCur::Iion];
    return oss.str();

}


std::string Bueno::PrintBlockCoeffs() const
{
    using namespace BueCur;

    // Create output string stream to pass the parameters and their values.
    std::ostringstream oss;
    oss.precision(15);
    oss << "Ifi = " << this->block_coeff_[Ifi] << "\n";
    oss << "Iso = " << this->block_coeff_[Iso] << "\n";
    oss << "Isi = " << this->block_coeff_[Isi];
    return oss.str();
}


} // End of namespace ELECTRA
